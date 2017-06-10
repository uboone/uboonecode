////////////////////////////////////////////////////////////////////////
/// \file   MCC7Deconvolution.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/CalData/DeconTools/MCC7Deconvolution.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "uboone/CalData/DeconTools/BaselineStandard.h"
#include "uboone/CalData/DeconTools/BaselineMostProbAve.h"

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{

//----------------------------------------------------------------------
// Constructor.
MCC7Deconvolution::MCC7Deconvolution(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
MCC7Deconvolution::~MCC7Deconvolution()
{
}
    
void MCC7Deconvolution::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
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
    
void MCC7Deconvolution::Deconvolve(IROIFinder::Waveform const&       waveform,
                                   raw::ChannelID_t                  channel,
                                   IROIFinder::CandidateROIVec const& roiVec,
                                   recob::Wire::RegionsOfInterest_t& ROIVec) const
{
    // The goal of this function is to reproduce "exactly" the operation of the deconvolution process in MCC7
    // hence the copying over of some of the code that has been pushed into external tools.
    
    // The size of the input waveform **should** be the raw buffer size
    size_t dataSize = waveform.size();
    
    // Make sure the deconvolution size is set correctly (this will probably be a noop after first call)
    fSignalShaping->SetDecon(dataSize, channel);
    
    size_t transformSize = fFFT->FFTSize();
    
    // now make a buffer to contain the waveform which will be of the right size
    std::vector<float> rawAdcLessPedVec;
    
    rawAdcLessPedVec.resize(transformSize,0.);
    
    size_t binOffset = transformSize > dataSize ? (transformSize - dataSize) / 2 : 0;
    double deconNorm = fSignalShaping->GetDeconNorm();
    
    // Copy the input (assumed pedestal subtracted) waveforms into our zero padded deconvolution buffer
    std::copy(waveform.begin(),waveform.end(),rawAdcLessPedVec.begin()+binOffset);
    
    // Strategy is to run deconvolution on the entire channel and then pick out the ROI's we found above
    fSignalShaping->Deconvolute(channel,rawAdcLessPedVec);
    
    std::vector<float> holder;
    
    for(size_t roiIdx = 0; roiIdx < roiVec.size(); roiIdx++)
    {
        const auto roi = roiVec[roiIdx];
        
        // First up: copy out the relevent ADC bins into the ROI holder
        size_t roiLen = roi.second - roi.first;
        
        holder.resize(roiLen);
        
        std::copy(rawAdcLessPedVec.begin()+binOffset+roi.first, rawAdcLessPedVec.begin()+binOffset+roi.second, holder.begin());
        std::transform(holder.begin(),holder.end(),holder.begin(),[deconNorm](float& deconVal){return deconVal/deconNorm;});
        
        // Now we do the baseline determination (and I'm left wondering if there is a better way using the entire waveform?)
        float base = fBaseline->GetBaseline(holder, channel, 0, roiLen);
        
        std::transform(holder.begin(),holder.end(),holder.begin(),[base](float& adcVal){return adcVal - base;});
        
        // apply wire-by-wire calibration
        if (fDodQdxCalib){
            if(fdQdxCalib.find(channel) != fdQdxCalib.end()){
                float constant = fdQdxCalib.at(channel);
                //std::cout<<channel<<" "<<constant<<std::endl;
                for (size_t iholder = 0; iholder < holder.size(); ++iholder){
                    holder[iholder] *= constant;
                }
            }
        }

        // add the range into ROIVec
        ROIVec.add_range(roi.first, std::move(holder));
    }
    
    return;
}

    
void MCC7Deconvolution::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "MCC7DeconvolutionPlane_" + std::to_string(fPlane);
 
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
 
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fMCC7DeconvolutionVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "MCC7DeconvolutionPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "MCC7Deconvolution;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fMCC7DeconvolutionVec.at(bin).Re());
    }
*/
    
    return;
}
    
}
