////////////////////////////////////////////////////////////////////////
/// \file   ROIDeconvolution.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/CalData/DeconTools/ROIDeconvolution.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "uboone/CalData/DeconTools/BaselineStandard.h"
#include "uboone/CalData/DeconTools/BaselineMostProbAve.h"

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{
    
//----------------------------------------------------------------------
// Constructor.
ROIDeconvolution::ROIDeconvolution(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIDeconvolution::~ROIDeconvolution()
{
}
    
void ROIDeconvolution::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fFFTSize     = pset.get< size_t >("FFTSize"                );
    
    //wire-by-wire calibration
    fDodQdxCalib = pset.get< bool >("DodQdxCalib", false);
    
    if (fDodQdxCalib)
    {
        fdQdxCalibFileName = pset.get< std::string >("dQdxCalibFileName");
        std::string fullname;
        cet::search_path sp("FW_SEARCH_PATH");
        sp.find_file(fdQdxCalibFileName, fullname);
        
        if (fullname.empty())
        {
            std::cout << "Input file " << fdQdxCalibFileName << " not found" << std::endl;
            throw cet::exception("File not found");
        }
        else
            std::cout << "Applying wire-by-wire calibration using file " << fdQdxCalibFileName << std::endl;
        
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        
        while (std::getline(inFile,line))
        {
            unsigned int channel;
            float        constant;
            std::stringstream linestream(line);
            linestream >> channel >> constant;
            fdQdxCalib[channel] = constant;
            if (channel%1000==0) std::cout<<"Channel "<<channel<<" correction factor "<<fdQdxCalib[channel]<<std::endl;
        }
    }
    
    // Recover the baseline tool

    fhicl::ParameterSet pb = pset.get<fhicl::ParameterSet>("Baseline");
    std::string pb_type = pb.get<std::string>("tool_type");
    if(pb_type == std::string("BaselineStandard")) {
      fBaseline  = std::unique_ptr<uboone_tool::IBaseline>(new uboone_tool::BaselineStandard(pb));
    }
    else if(pb_type == std::string("BaselineMostProbAve")) {
      fBaseline  = std::unique_ptr<uboone_tool::IBaseline>(new uboone_tool::BaselineMostProbAve(pb));
    }
    else {
      throw art::Exception(art::errors::Configuration) << "Unknown baseline tool type" << pb_type;
    }
    
    // Get signal shaping service.
    fSignalShaping = art::ServiceHandle<util::SignalShapingServiceMicroBooNE>();
    fFFT           = art::ServiceHandle<util::LArFFT>();
    
    return;
}
void ROIDeconvolution::Deconvolve(const IROIFinder::Waveform&        waveform,
                                  raw::ChannelID_t                   channel,
                                  IROIFinder::CandidateROIVec const& roiVec,
                                  recob::Wire::RegionsOfInterest_t&  ROIVec) const
{
    double deconNorm = fSignalShaping->GetDeconNorm();
    
    // And now process them
    for(auto const& roi : roiVec)
    {
        // First up: copy out the relevent ADC bins into the ROI holder
        size_t roiLen = roi.second - roi.first;
        
        // We want the deconvolution buffer size to be a power of 2 in length
        // to facilitate the FFT
        size_t deconSize = fFFTSize;
        
        while(1)
        {
            if (roiLen > deconSize) deconSize *= 2;
            else break;
        }
        
        // In theory, most ROI's are around the same size so this should mostly be a noop
        fSignalShaping->SetDecon(deconSize, channel);
        
        deconSize = fFFT->FFTSize();
        
        std::vector<float> holder(deconSize);
        
        // Pad with zeroes if the deconvolution buffer is larger than the input waveform
        if (deconSize > waveform.size()) holder.resize(deconSize, 0.);
        
        // Watch for the case where the input ROI is long enough to want an deconvolution buffer that is
        // larger than the input waveform.
        size_t maxActualSize = std::min(deconSize, waveform.size());
        
        // Extend the ROI to accommodate the extra bins for the FFT
        // The idea is to try to center the desired ROI in the buffer used by deconvolution
        size_t halfLeftOver = (maxActualSize - roiLen) / 2;           // Number bins either side of ROI
        int    roiStartInt  = halfLeftOver;                           // Start in the buffer of the ROI
        int    roiStopInt   = halfLeftOver + roiLen;                  // Stop in the buffer of the ROI
        int    firstOffset  = roi.first - halfLeftOver;               // Offset into the ADC vector of buffer start
        int    secondOffset = roi.second + halfLeftOver + roiLen % 2; // Offset into the ADC vector of buffer end
        
        // Check for the two edge conditions - starting before the ADC vector or running off the end
        // In either case we shift the actual roi within the FFT buffer
        // First is the case where we would be starting before the ADC vector
        if (firstOffset < 0)
        {
            roiStartInt  += firstOffset;  // remember that firstOffset is negative
            roiStopInt   += firstOffset;
            secondOffset -= firstOffset;
            firstOffset   = 0;
        }
        // Second is the case where we would overshoot the end
        else if (size_t(secondOffset) > waveform.size())
        {
            size_t overshoot = secondOffset - waveform.size();
            
            roiStartInt  += overshoot;
            roiStopInt   += overshoot;
            firstOffset  -= overshoot;
            secondOffset  = waveform.size();
        }
        
        size_t roiStart(roiStartInt);
        size_t roiStop(roiStopInt);
        
        // Fill the buffer and do the deconvolution
        std::copy(waveform.begin()+firstOffset, waveform.begin()+secondOffset, holder.begin());
        
        // Deconvolute the raw signal
        fSignalShaping->Deconvolute(channel,holder);
        
        // "normalize" the vector
        std::transform(holder.begin(),holder.end(),holder.begin(),[deconNorm](float& deconVal){return deconVal/deconNorm;});
        
        // Now we do the baseline determination and correct the ROI
        float base = fBaseline->GetBaseline(holder, channel, roiStart, roiLen);
        
        std::transform(holder.begin(),holder.end(),holder.begin(),[base](float& adcVal){return adcVal - base;});
        
        // Get rid of the leading and trailing "extra" bins needed to keep the FFT happy
        if (roiStart > 0) std::copy(holder.begin() + roiStart, holder.begin() + roiStop, holder.begin());
        
        // Resize the holder to the ROI length
        holder.resize(roiLen);
        
        // apply wire-by-wire calibration
        if (fDodQdxCalib)
        {
            if(fdQdxCalib.find(channel) != fdQdxCalib.end())
            {
                float constant = fdQdxCalib.at(channel);
                
                for (size_t iholder = 0; iholder < holder.size(); ++iholder) holder[iholder] *= constant;
            }
        }

        // add the range into ROIVec
        ROIVec.add_range(roi.first, std::move(holder));
    } // loop over candidate roi's
    
    return;
}
    

    
void ROIDeconvolution::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "ROIDeconvolutionPlane_" + std::to_string(fPlane);
 
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
 
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fROIDeconvolutionVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "ROIDeconvolutionPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "ROIDeconvolution;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fROIDeconvolutionVec.at(bin).Re());
    }
*/
    
    return;
}
    
}
