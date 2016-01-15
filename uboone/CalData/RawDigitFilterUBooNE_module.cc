////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitFilterUBooNE
// Module Type: producer
// File:        RawDigitFilterUBooNE_module.cc
//
//              The intent of this module is to filter out "bad" channels
//              in an input RawDigit data stream. In the current implementation,
//              "bad" is defined as the truncated rms for a channel crossing
//              a user controlled threshold
//
// Configuration parameters:
//
// DigitModuleLabel      - the source of the RawDigit collection
// TruncMeanFraction     - the fraction of waveform bins to discard when
//                         computing the means and rms
// RMSRejectionCutHi     - vector of maximum allowed rms values to keep channel
// RMSRejectionCutLow    - vector of lowest allowed rms values to keep channel
// RMSSelectionCut       - vector of rms values below which to not correct
// TheChoseWire          - Wire chosen for "example" hists
// MaxPedestalDiff       - Baseline difference to pedestal to flag
// SmoothCorrelatedNoise - Turns on the correlated noise suppression
// NumWiresToGroup       - When removing correlated noise, # wires to group
// FillHistograms        - Turn on histogram filling for diagnostics
// RunFFTInputWires      - FFT analyze the input RawDigits if true - diagnostics
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
// TruncateTicks:        - Provide mechanism to truncate a readout window to a smaller size
// WindowSize:           - The desired size of the output window
// NumTicksToDropFront:  - The number ticks dropped off the front of the original window
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on August 17, 2015
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "Utilities/SimpleTimeService.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"

#include "NoiseFilterAlgs/RawDigitNoiseFilterDefs.h"
#include "NoiseFilterAlgs/RawDigitBinAverageAlg.h"
#include "NoiseFilterAlgs/RawDigitCharacterizationAlg.h"
#include "NoiseFilterAlgs/RawDigitCorrelatedCorrectionAlg.h"
#include "NoiseFilterAlgs/RawDigitFilterAlg.h"

class RawDigitFilterUBooNE : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit RawDigitFilterUBooNE(fhicl::ParameterSet const & pset);
    virtual ~RawDigitFilterUBooNE();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:
    
    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, raw::ChannelID_t&, caldata::RawDigitVector&, float, float);

    // Fcl parameters.
    std::string          fDigitModuleLabel;      ///< The full collection of hits
    float                fTruncMeanFraction;     ///< Fraction for truncated mean
    std::vector<float>   fRmsRejectionCutHi;     ///< Maximum rms for input channels, reject if larger
    std::vector<float>   fRmsRejectionCutLow;    ///< Minimum rms to consider channel "alive"
    std::vector<float>   fRmsSelectionCut;       ///< Don't use/apply correction to wires below this
    std::vector<short>   fMinMaxSelectionCut;    ///< Plane by plane cuts for spread cut
    unsigned int         fTheChosenWire;         ///< For example hist
    double               fMaxPedestalDiff;       ///< Max pedestal diff to db to warn
    bool                 fApplyBinAverage;       ///< Do bin averaging to get rid of high frequency noise
    std::vector<size_t>  fBinsToAverage;         ///< # bins to average by view
    bool                 fApplyTopHatFilter;     ///< Apply the top hat filter
    size_t               fStructuringElement;    ///< Structuring element to use with Top Hat filter
    bool                 fSmoothCorrelatedNoise; ///< Should we smooth the noise?
    bool                 fApplyCorSmoothing;     ///< Attempt to smooth the correlated noise correction?
    bool                 fApplyFFTCorrection;    ///< Use an FFT to get the correlated noise correction
    bool                 fFillFFTHistograms;     ///< Fill associated FFT histograms
    std::vector<size_t>  fFFTHistsWireGroup;     ///< Wire Group to pick on
    std::vector<size_t>  fFFTNumHists;           ///< Number of hists total per view
    std::vector<double>  fFFTHistsStartTick;     ///< Starting tick for histograms
    std::vector<double>  fFFTMinPowerThreshold;  ///< Threshold for trimming FFT power spectrum
    std::vector<size_t>  fNumWiresToGroup;       ///< If smoothing, the number of wires to look at
    bool                 fFillHistograms;        ///< if true then will fill diagnostic hists
    bool                 fRunFFTInput;           ///< Should we run FFT's on input wires?
    bool                 fRunFFTCorrected;       ///< Should we run FFT's on corrected wires?
    bool                 fTruncateTicks;         ///< If true then drop channels off ends of wires
    unsigned int         fWindowSize;            ///< # ticks to keep in window
    unsigned int         fNumTicksToDropFront;   ///< # ticks to drop from front of waveform

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Correction algorithms
    caldata::RawDigitBinAverageAlg               fBinAverageAlg;
    caldata::RawDigitCharacterizationAlg         fCharacterizationAlg;
    caldata::RawDigitCorrelatedCorrectionAlg     fCorCorrectAlg;
    caldata::RawDigitFilterAlg                   fFilterAlg;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>            fGeometry;             ///< pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;   ///< Detector properties service
    const lariov::IDetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};

DEFINE_ART_MODULE(RawDigitFilterUBooNE)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitFilterUBooNE::RawDigitFilterUBooNE(fhicl::ParameterSet const & pset) :
                      fNumEvent(0),
                      fBinAverageAlg(pset),
                      fCharacterizationAlg(pset),
                      fCorCorrectAlg(pset),
                      fFilterAlg(pset),
                      fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())

{
    reconfigure(pset);
    produces<std::vector<raw::RawDigit> >();

    // Report.
    mf::LogInfo("RawDigitFilterUBooNE") << "RawDigitFilterUBooNE configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFilterUBooNE::~RawDigitFilterUBooNE()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFilterUBooNE::reconfigure(fhicl::ParameterSet const & pset)
{
    fDigitModuleLabel      = pset.get<std::string>        ("DigitModuleLabel",                                        "daq");
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fRmsRejectionCutHi     = pset.get<std::vector<float>> ("RMSRejectionCutHi",     std::vector<float>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<float>> ("RMSRejectionCutLow",    std::vector<float>() = {0.70,0.70,0.70});
    fRmsSelectionCut       = pset.get<std::vector<float>> ("RMSSelectionCut",       std::vector<float>() = {1.40,1.40,1.00});
    fMinMaxSelectionCut    = pset.get<std::vector<short>> ("MinMaxSelectionCut",        std::vector<short>() = {13, 13, 11});
    fTheChosenWire         = pset.get<unsigned int>       ("TheChosenWire",                                            1200);
    fMaxPedestalDiff       = pset.get<double>             ("MaxPedestalDiff",                                           10.);
    fApplyBinAverage       = pset.get<bool>               ("ApplyBinAverage",                                          true);
    fBinsToAverage         = pset.get<std::vector<size_t>>("NumBinsToAverage",              std::vector<size_t>() = {2,2,2});
    fApplyTopHatFilter     = pset.get<bool>               ("ApplyTopHatFilter",                                        true);
    fStructuringElement    = pset.get<size_t>             ("StructuringElement",                                         30);
    fSmoothCorrelatedNoise = pset.get<bool>               ("SmoothCorrelatedNoise",                                    true);
    fApplyCorSmoothing     = pset.get<bool>               ("ApplyCorSmoothing",                                        true);
    fApplyFFTCorrection    = pset.get<bool>               ("ApplyFFTCorrection",                                       true);
    fFillFFTHistograms     = pset.get<bool>               ("FillFFTHistograms",                                        true);
    fFFTHistsWireGroup     = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fFFTNumHists           = pset.get<std::vector<size_t>>("FFTNumWaveHistograms",       std::vector<size_t>() = {10,48,48});
    fFFTHistsStartTick     = pset.get<std::vector<double>>("FFTWaveHistsStartTick", std::vector<double>() = {96.,96.,7670.});
    fFFTMinPowerThreshold  = pset.get<std::vector<double>>("FFTPowerThreshold",     std::vector<double>() = {100.,75.,500.});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                          false);
    fRunFFTInput           = pset.get<bool>               ("RunFFTInputWires",                                        false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                    false);
    fTruncateTicks         = pset.get<bool>               ("TruncateTicks",                                           false);
    fWindowSize            = pset.get<size_t>             ("WindowSize",                                               6400);
    fNumTicksToDropFront   = pset.get<size_t>             ("NumTicksToDropFront",                                      2250);
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitFilterUBooNE::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    fCharacterizationAlg.initializeHists(tfs);
    fCorCorrectAlg.initializeHists(tfs);
    fFilterAlg.initializeHists(tfs);
    
    return;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void RawDigitFilterUBooNE::produce(art::Event & event)
{
    ++fNumEvent;
    
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<raw::RawDigit> > filteredRawDigit(new std::vector<raw::RawDigit>);
    
    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    event.getByLabel(fDigitModuleLabel, digitVecHandle);
    
    // Require a valid handle
    if (digitVecHandle.isValid())
    {
        unsigned int maxChannels    = fGeometry->Nchannels();
        unsigned int maxTimeSamples = fDetectorProperties->NumberTimeSamples();
        
        // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
        // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
        std::vector<const raw::RawDigit*> rawDigitVec;
        
        // Ugliness to fill the pointer vector...
        for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.push_back(&digitVecHandle->at(idx)); //art::Ptr<raw::RawDigit>(digitVecHandle, idx).get());
        
        // Sort (use a lambda to sort by channel id)
        std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
    
        // Ok, to do the correlated noise removal we are going to need a rather impressive data structure...
        // Because we need to unpack each wire's data, we will need to "explode" it out into a data structure
        // here... with the good news that we'll release the memory at the end of the module so should not
        // impact downstream processing (I hope!).
        // What we are going to do is make a vector over views of vectors over wires of vectors over time samples
        //std::vector<RawDigitVector> rawDataWireTimeVec;
        std::vector<caldata::RawDigitVector> rawDataWireTimeVec;
        std::vector<float>                   truncMeanWireVec;
        std::vector<float>                   truncRmsWireVec;
        std::vector<short>                   meanWireVec;
        std::vector<short>                   medianWireVec;
        std::vector<short>                   modeWireVec;
        std::vector<float>                   skewnessWireVec;
        std::vector<float>                   fullRmsWireVec;
        std::vector<short>                   minMaxWireVec;
        std::vector<float>                   neighborRatioWireVec;
        std::vector<float>                   pedCorWireVec;
        std::vector<raw::ChannelID_t>        channelWireVec;
        caldata::GroupToDigitIdxPairMap      groupToDigitIdxPairMap;
        
        // Declare a temporary digit holder and resize it if downsizing the waveform
        caldata::RawDigitVector tempVec;
    
        // Commence looping over raw digits
        for(const auto& rawDigit : rawDigitVec)
        {
            raw::ChannelID_t channel = rawDigit->Channel();
        
            bool goodChan(true);
        
            // The below try-catch block may no longer be necessary
            // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids;
            try {
                wids = fGeometry->ChannelToWire(channel);
            }
            catch(...)
            {
                //std::cout << "===>> Found illegal channel with id: " << channel << std::endl;
                goodChan = false;
            }
        
            if (channel >= maxChannels || !goodChan) continue;

        
            // Recover plane and wire in the plane
            unsigned int view = wids[0].Plane;
            unsigned int wire = wids[0].Wire;
        
            unsigned int dataSize = rawDigit->Samples();
            unsigned int wireIdx  = wire % fNumWiresToGroup[view];
            
            // Cross check that our storage arrays are the correct size
            // (note there is a possible boundary issue here that we are going to ignore...)
            if (rawDataWireTimeVec.size() != fNumWiresToGroup[view])
            {
                // For each view we need to presize the vector to the number of wires
                rawDataWireTimeVec.resize(fNumWiresToGroup[view]);
                truncMeanWireVec.resize(fNumWiresToGroup[view]);
                truncRmsWireVec.resize(fNumWiresToGroup[view]);
                meanWireVec.resize(fNumWiresToGroup[view]);
                medianWireVec.resize(fNumWiresToGroup[view]);
                modeWireVec.resize(fNumWiresToGroup[view]);
                skewnessWireVec.resize(fNumWiresToGroup[view]);
                fullRmsWireVec.resize(fNumWiresToGroup[view]);
                minMaxWireVec.resize(fNumWiresToGroup[view]);
                neighborRatioWireVec.resize(fNumWiresToGroup[view]);
                pedCorWireVec.resize(fNumWiresToGroup[view]);
                channelWireVec.resize(fNumWiresToGroup[view]);
                groupToDigitIdxPairMap.clear();
            }
        
            // vector holding uncompressed adc values
            std::vector<short>& rawadc = rawDataWireTimeVec[wireIdx];
        
            channelWireVec[wireIdx] = channel;
            
            // If we are trying to truncate the incoming RawDigit collection then we need to do so when we extract from the input raw digits
            // This causes a small division here...
            if (fTruncateTicks)
            {
                maxTimeSamples = fWindowSize;
                
                if (rawadc.size()  != maxTimeSamples) rawadc.resize(maxTimeSamples);
                if (tempVec.size() != dataSize)       tempVec.resize(dataSize);
                
                // And now uncompress
                raw::Uncompress(rawDigit->ADCs(), tempVec, rawDigit->Compression());
                
                std::copy(tempVec.begin() + fNumTicksToDropFront, tempVec.begin() + fNumTicksToDropFront + fWindowSize, rawadc.begin());
            }
            else
            {
                maxTimeSamples = dataSize;
                
                if (rawadc.size() != dataSize) rawadc.resize(maxTimeSamples);
                
                // And now uncompress
                raw::Uncompress(rawDigit->ADCs(), rawadc, rawDigit->Compression());
            }
            
            // Apply the high frequency filter
            if (fApplyBinAverage) fBinAverageAlg.doTwoBinAverage(rawadc);
            
            // Apply the top hat filter
//            if (fApplyTopHatFilter) doTopHatFilter(rawadc);
//            if (fApplyTopHatFilter) doAdaptiveFilter(rawadc);
            
            // Get the kitchen sink
            fCharacterizationAlg.getWaveformParams(rawadc,
                                                   channel,
                                                   view,
                                                   wire,
                                                   truncMeanWireVec[wireIdx],
                                                   truncRmsWireVec[wireIdx],
                                                   meanWireVec[wireIdx],
                                                   medianWireVec[wireIdx],
                                                   modeWireVec[wireIdx],
                                                   skewnessWireVec[wireIdx],
                                                   fullRmsWireVec[wireIdx],
                                                   minMaxWireVec[wireIdx],
                                                   neighborRatioWireVec[wireIdx],
                                                   pedCorWireVec[wireIdx]);
            
            // If we are not performing noise corrections then we are done with this wire
            // Store it and move on
            if (!fSmoothCorrelatedNoise)
            {
                // Filter out the very high noise wires
                if (truncRmsWireVec[wireIdx] < fRmsRejectionCutHi[view])
                    saveRawDigits(filteredRawDigit, channel, rawadc, truncMeanWireVec[wireIdx], truncRmsWireVec[wireIdx]);
                else
                {
                    // Eventually we'll interface to some sort of channel status communication mechanism.
                    // For now use the log file
                    mf::LogInfo("RawDigitFilterUBooNE") <<  "--> Rejecting channel for large rms, channel: " << channel << ", rmsVal: " << pedCorWireVec[wireIdx] << ", truncMean: " << truncMeanWireVec[wireIdx] << ", pedestal: " << pedCorWireVec[wireIdx] << std::endl;
                }
                
                continue;
            }
            
            // Add this wire to the map and try to do some classification here
            fCharacterizationAlg.classifyRawDigitVec(rawadc, view, wire, truncRmsWireVec[wireIdx], minMaxWireVec[wireIdx], meanWireVec[wireIdx], skewnessWireVec[wireIdx], neighborRatioWireVec[wireIdx], groupToDigitIdxPairMap);

            // Are we at the correct boundary for dealing with the noise?
            if (!((wireIdx + 1) % fNumWiresToGroup[view]))
            {
                int baseWireIdx = wire - wire % fNumWiresToGroup[view];

                // Now go through the groups to remove correlated noise in those groups
                for(auto& groupToDigitIdxPair : groupToDigitIdxPairMap)
                {
                    fCorCorrectAlg.removeCorrelatedNoise(groupToDigitIdxPair.second,
                                                         view,
                                                         truncMeanWireVec,
                                                         truncRmsWireVec,
                                                         minMaxWireVec,
                                                         meanWireVec,
                                                         skewnessWireVec,
                                                         neighborRatioWireVec,
                                                         pedCorWireVec);
                }
                
                // One more pass through to store the good channels
                for (size_t locWireIdx = 0; locWireIdx < fNumWiresToGroup[view]; locWireIdx++)
                {
                    // Try baseline correction?
//                    if (fApplyTopHatFilter && view < 2 && skewnessWireVec[locWireIdx] != 0.) doTopHatFilter(rawDataWireTimeVec[locWireIdx]);
                    if (fApplyTopHatFilter && view != 2 && skewnessWireVec[locWireIdx] > 0.)
                    {
                        //doAdaptiveFilter(rawDataWireTimeVec[locWireIdx]);
                        fFilterAlg.doTopHatFilter(rawDataWireTimeVec[locWireIdx], baseWireIdx + locWireIdx);
                    }
                    
                    // recalculate rms for the output
                    double rmsVal   = 0.;
                    double pedestal = truncMeanWireVec[locWireIdx];
                    double pedCor   = pedCorWireVec[locWireIdx];
                    
                    caldata::RawDigitVector& rawDataVec = rawDataWireTimeVec[locWireIdx];
                    
                    for(const auto& adcVal : rawDataVec)
                    {
                        double adcLessPed = adcVal - pedestal + pedCor;
                        rmsVal += adcLessPed * adcLessPed;
                    }
                        
                    rmsVal = std::sqrt(std::max(0.,rmsVal / double(rawDataVec.size())));
                    
                    // The ultra high noise channels are simply zapped
                    if (rmsVal < fRmsRejectionCutHi[view]) // && ImAGoodWire(view,baseWireIdx + locWireIdx))
                    {
                        saveRawDigits(filteredRawDigit, channelWireVec[locWireIdx], rawDataVec, pedestal, rmsVal);
                    }
                    else
                    {
                        mf::LogInfo("RawDigitFilterUBooNE") <<  "--> Rejecting channel for large rms, channel: " << channel << ", rmsVal: " << truncRmsWireVec[wireIdx] << ", truncMean: " << truncMeanWireVec[wireIdx] << ", pedestal: " << pedCorWireVec[wireIdx] << std::endl;
                    }
                }
                
                groupToDigitIdxPairMap.clear();
            }
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(filteredRawDigit));
}

void RawDigitFilterUBooNE::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                         raw::ChannelID_t&                             channel,
                                         caldata::RawDigitVector&                      rawDigitVec,
                                         float                                         pedestal,
                                         float                                         rms)
{
    filteredRawDigit->emplace_back(raw::RawDigit(channel, rawDigitVec.size(), rawDigitVec, raw::kNone));
    filteredRawDigit->back().SetPedestal(pedestal,rms);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void RawDigitFilterUBooNE::endJob()
{
    mf::LogInfo("RawDigitFilterUBooNE") << "Looked at " << fNumEvent << " events" << std::endl;
}
