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

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{

//----------------------------------------------------------------------
// Constructor.
MCC7Deconvolution::MCC7Deconvolution(const fhicl::ParameterSet& pset) :
                      fROIPropertiesAlg(pset.get<fhicl::ParameterSet>("ROIPropertiesAlg"))
{
    configure(pset);
}
    
MCC7Deconvolution::~MCC7Deconvolution()
{
}
    
void MCC7Deconvolution::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fDoBaselineSub                       = pset.get< bool >("DoBaselineSub"                      );
    fMinROIAverageTickThreshold          = pset.get<float> ("MinROIAverageTickThreshold",    -0.5);
    fDoBaselineSub_WaveformPropertiesAlg = pset.get< bool >("DoBaselineSub_WaveformPropertiesAlg");
    
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
        bool  baseSet(false);
        float base(0.);
        if(fDoBaselineSub) // && fPreROIPad[thePlane] > 0 ) <-- this part doesn't make sense?
        {
            //1. Check Baseline match?
            // If not, include next ROI(if none, go to the end of signal)
            // If yes, proceed
            size_t binsToAve(20);
            float  basePre  = std::accumulate(holder.begin(),holder.begin()+binsToAve,0.) / float(binsToAve);
            float  basePost = std::accumulate(holder.end()-binsToAve,holder.end(),0.) / float(binsToAve);
            
            // emulate method for refining baseline from original version of CalWireROI
            float deconNoise = 1.26491 * fSignalShaping->GetDeconNoise(channel);    // 4./sqrt(10) * noise
            
            // If the estimated baseline from the front of the roi does not agree well with that from the end
            // of the roi then we'll extend the roi hoping for good agreement
            if (!(fabs(basePre - basePost) < deconNoise))
            {
                int   nTries(0);
                
                // get start of roi and find the maximum we can extend to
                std::vector<float>::iterator rawAdcRoiStartItr = rawAdcLessPedVec.begin() + binOffset + roi.first;
                std::vector<float>::iterator rawAdcMaxItr      = rawAdcLessPedVec.end()   - binOffset;
                
                // if this is not the last roi then limit max range to start of next roi
                if (roiIdx < roiVec.size() - 1)
                    rawAdcMaxItr = rawAdcLessPedVec.begin() + binOffset + roiVec[roiIdx+1].first;
                
                // Basically, allow roi to be extended until we get good agreement unless it seems pointless
                while (!(fabs(basePre - basePost) < deconNoise) && nTries++ < 3)
                {
                    size_t nBinsToAdd(100);
                    int    nBinsLeft      = std::distance(rawAdcRoiStartItr+roiLen,rawAdcMaxItr) > 0
                    ? std::distance(rawAdcRoiStartItr+roiLen,rawAdcMaxItr) : 0;
                    size_t roiLenAddition = std::min(nBinsToAdd, size_t(nBinsLeft));
                    
                    if (roiLenAddition > 0)
                    {
                        std::vector<float> additionVec(roiLenAddition);
                        
                        std::copy(rawAdcRoiStartItr + roiLen, rawAdcRoiStartItr + roiLen + roiLenAddition, additionVec.begin());
                        
                        holder.resize(holder.size() + roiLenAddition);
                        std::transform(additionVec.begin(),additionVec.end(),holder.begin() + roiLen,[deconNorm](float& deconVal){return deconVal/deconNorm;});
                        
                        basePost = std::accumulate(holder.end()-binsToAve,holder.end(),0.) / float(binsToAve);
                        
                        roiLen = holder.size();
                    }
                    else break;
                }
            }
            
            baseSet = true;
            base    = SubtractBaseline(holder, basePre,basePost,roi.first,roiLen,dataSize);
        } // fDoBaselineSub ...
        else if(fDoBaselineSub_WaveformPropertiesAlg)
        {
            baseSet = true;
            base    = fROIPropertiesAlg.GetWaveformPedestal(holder);
        }
        
        if (baseSet) std::transform(holder.begin(),holder.end(),holder.begin(),[base](float& adcVal){return adcVal - base;});
        
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
        
        //wes 23.12.2016 --- sum up the roi, and if it's very negative get rid of it
        float average_val = std::accumulate(holder.begin(),holder.end(),0.0) / holder.size();
        float min = *std::min_element(holder.begin(),holder.end());
        float max = *std::max_element(holder.begin(),holder.end());
        if(average_val>fMinROIAverageTickThreshold || std::abs(min)<std::abs(max)){
            // add the range into ROIVec
            ROIVec.add_range(roi.first, std::move(holder));
        }
    }
    
    return;
}
    
    
float MCC7Deconvolution::SubtractBaseline(std::vector<float>& holder,
                                          float               basePre,
                                          float               basePost,
                                          size_t              roiStart,
                                          size_t              roiLen,
                                          size_t              dataSize) const
{
    float base=0;
    
    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20){
        base = basePost;
        // can not trust the later part
    }else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20){
        base = basePre;
        // can trust both
    }else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20){
        if (fabs(basePre-basePost)<3){
            base = (basePre+basePost)/2.;
        }else{
            if (basePre < basePost){
                base = basePre;
            }else{
                base = basePost;
            }
        }
        // can not use both
    }else{
        float min = 0,max=0;
        for (unsigned int bin = 0; bin < roiLen; bin++){
            if (holder[bin] > max) max = holder[bin];
            if (holder[bin] < min) min = holder[bin];
        }
        int nbin = max - min;
        if (nbin!=0){
            TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
            for (unsigned int bin = 0; bin < roiLen; bin++){
                h1->Fill(holder[bin]);
            }
            float ped = h1->GetMaximum();
            float ave=0,ncount = 0;
            for (unsigned int bin = 0; bin < roiLen; bin++){
                if (fabs(holder[bin]-ped)<2){
                    ave +=holder[bin];
                    ncount ++;
                }
            }
            if (ncount==0) ncount=1;
            ave = ave/ncount;
            h1->Delete();
            base = ave;
        }
    }
    
    return base;
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
