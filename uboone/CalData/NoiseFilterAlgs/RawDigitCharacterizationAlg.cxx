
#include "RawDigitCharacterizationAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cmath>
#include <algorithm>

namespace caldata
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitCharacterizationAlg::RawDigitCharacterizationAlg(fhicl::ParameterSet const & pset) :
                      fHistsInitialized(false),
                      fFirstEvent(true),
                      fChannelGroups(pset),
                      fPedestalRetrievalAlg(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider())

{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitCharacterizationAlg") << "RawDigitCharacterizationAlg configured\n";
}
    
//----------------------------------------------------------------------------
/// Destructor.
RawDigitCharacterizationAlg::~RawDigitCharacterizationAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitCharacterizationAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fRmsRejectionCutHi     = pset.get<std::vector<float>> ("RMSRejectionCutHi",     std::vector<float>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<float>> ("RMSRejectionCutLow",    std::vector<float>() = {0.70,0.70,0.70});
    fRmsSelectionCut       = pset.get<std::vector<float>> ("RMSSelectionCut",       std::vector<float>() = {1.40,1.40,1.00});
    fMinMaxSelectionCut    = pset.get<std::vector<short>> ("MinMaxSelectionCut",        std::vector<short>() = {13, 13, 11});
    fTheChosenWire         = pset.get<unsigned int>       ("TheChosenWire",                                            1200);
    fMaxPedestalDiff       = pset.get<double>             ("MaxPedestalDiff",                                           10.);
    fHistsWireGroup        = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                          false);
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitCharacterizationAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fFillHistograms)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fAdcCntHist[0]    = tfs->make<TH1D>("CntUPlane", ";#adc",  200, 9000., 10000.);
        fAdcCntHist[1]    = tfs->make<TH1D>("CntVPlane", ";#adc",  200, 9000., 10000.);
        fAdcCntHist[2]    = tfs->make<TH1D>("CntWPlane", ";#adc",  200, 9000., 10000.);
        fAveValHist[0]    = tfs->make<TH1D>("AveUPlane", ";Ave",   120,  -30.,    30.);
        fAveValHist[1]    = tfs->make<TH1D>("AveVPlane", ";Ave",   120,  -30.,    30.);
        fAveValHist[2]    = tfs->make<TH1D>("AveWPlane", ";Ave",   120,  -30.,    30.);
        fRmsValHist[0]    = tfs->make<TH1D>("RmsUPlane", ";RMS",   200,    0.,    50.);
        fRmsValHist[1]    = tfs->make<TH1D>("RmsVPlane", ";RMS",   200,    0.,    50.);
        fRmsValHist[2]    = tfs->make<TH1D>("RmsWPlane", ";RMS",   200,    0.,    50.);
        fPedValHist[0]    = tfs->make<TH1D>("PedUPlane", ";Ped",   200,  1950,  2150.);
        fPedValHist[1]    = tfs->make<TH1D>("PedVPlane", ";Ped",   200,  1950,  2150.);
        fPedValHist[2]    = tfs->make<TH1D>("PedWPlane", ";Ped",   200,   350,   550.);
    
        fRmsValProf[0]    = tfs->make<TProfile>("RmsUPlaneProf",    ";Wire #",  2400, 0., 2400., 0., 100.);
        fRmsValProf[1]    = tfs->make<TProfile>("RmsVPlaneProf",    ";Wire #",  2400, 0., 2400., 0., 100.);
        fRmsValProf[2]    = tfs->make<TProfile>("RmsWPlaneProf",    ";Wire #",  3456, 0., 3456., 0., 100.);
    
        fMinMaxValProf[0] = tfs->make<TProfile>("MinMaxUPlaneProf", ";Wire #",  2400, 0., 2400., 0., 200.);
        fMinMaxValProf[1] = tfs->make<TProfile>("MinMaxVPlaneProf", ";Wire #",  2400, 0., 2400., 0., 200.);
        fMinMaxValProf[2] = tfs->make<TProfile>("MinMaxWPlaneProf", ";Wire #",  3456, 0., 3456., 0., 200.);

        fPedValProf[0]    = tfs->make<TProfile>("PedUPlaneProf",    ";Wire #",  2400, 0., 2400., 1500., 2500.);
        fPedValProf[1]    = tfs->make<TProfile>("PedVPlaneProf",    ";Wire #",  2400, 0., 2400., 1500., 2500.);
        fPedValProf[2]    = tfs->make<TProfile>("PedWPlaneProf",    ";Wire #",  3456, 0., 3456.,    0., 1000.);
    
        fAverageHist[0]   = tfs->make<TH1D>("AverageU", ";Bin", 1000, 1500., 2500.);
        fAverageHist[1]   = tfs->make<TH1D>("AverageV", ";Bin", 1000, 1500., 2500.);
        fAverageHist[2]   = tfs->make<TH1D>("AverageW", ";Bin", 1000,    0., 1000.);
    
        fMinMaxProfiles.resize(3);
        fSkewnessProfiles.resize(3);
        fModeRatioProfiles.resize(3);
        
        for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
        {
            std::string minMaxName = "MinMax_" + std::to_string(planeIdx);
        
            fMinMaxProfiles[planeIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Min/Max Profiles;Wire", fNumWiresToGroup[planeIdx], 0., fNumWiresToGroup[planeIdx], 0., 200.);
        
            minMaxName = "Skewness_" + std::to_string(planeIdx);
        
            fSkewnessProfiles[planeIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Skewness Profiles;Wire", fNumWiresToGroup[planeIdx], 0., fNumWiresToGroup[planeIdx], -4., 4.);
        
            minMaxName = "ModeRatio_" + std::to_string(planeIdx);
        
            fModeRatioProfiles[planeIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Mode Ratio;Wire", fNumWiresToGroup[planeIdx], 0., fNumWiresToGroup[planeIdx], 0., 1.2);
        }
    
        fHistsInitialized = true;
    }
    
    return;
}

// Basic waveform mean, rms and pedestal offset
void RawDigitCharacterizationAlg::getWaveformParams(const RawDigitVector& rawWaveform,
                                                    unsigned int          channel,
                                                    unsigned int          view,
                                                    unsigned int          wire,
                                                    float&                truncMean,
                                                    float&                truncRms,
                                                    short&                mean,
                                                    short&                median,
                                                    short&                mode,
                                                    float&                skewness,
                                                    float&                rms,
                                                    short&                minMax,
                                                    float&                neighborRatio,
                                                    float&                pedCorVal) const
{
    // Recover the truncated mean and rms, plus the ped correction
    float minMaxFloat;
    
    getMeanRmsMinMaxAndPedCor(rawWaveform, channel, truncMean, truncRms, minMaxFloat, pedCorVal);
    
    // Convert to a short
    minMax = std::max(float(-4095),std::min(float(4095.),minMaxFloat));
    
    // We also want mean, median, rms, etc., for all ticks on the waveform
    std::vector<short> localTimeVec = rawWaveform;
    
    std::sort(localTimeVec.begin(),localTimeVec.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float realMean(float(std::accumulate(localTimeVec.begin(),localTimeVec.end(),0))/float(localTimeVec.size()));
    
    median = localTimeVec[localTimeVec.size()/2];
    mean   = std::round(realMean);
    
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(localTimeVec.size());
    
    std::transform(localTimeVec.begin(),localTimeVec.end(),adcLessPedVec.begin(),std::bind2nd(std::minus<short>(),mean));
    
    rms      = std::sqrt(std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.) / float(adcLessPedVec.size()));
    skewness = 3. * float(mean - median) / rms;
    
    // Final task is to get the mode and neighbor ratio
    // To do this we need to set up a vector of counts by ADC value...
    std::map<short,int> countMap;
    short               adcValMP(-4095);
    int                 countMP(0);
    
    for(const auto& adc : rawWaveform)
    {
        if (++countMap[adc] > countMP)
        {
            adcValMP = adc;
            countMP  = countMap[adc];
        }
    }
    
    short neighborSum(0);
    short leftNeighbor(countMP);
    short rightNeighbor(countMP);
    
    mode = adcValMP;
    
    if (countMap.find(adcValMP-1) != countMap.end())
    {
        leftNeighbor  = countMap.find(adcValMP-1)->second;
        neighborSum  += leftNeighbor;
    }
    
    if (countMap.find(adcValMP+1) != countMap.end())
    {
        rightNeighbor  = countMap.find(adcValMP+1)->second;
        neighborSum   += rightNeighbor;
    }
    
    neighborRatio = float(neighborSum) / float(2*countMP);

    neighborRatio = float(std::min(leftNeighbor,rightNeighbor)) / float(countMP);
    
    // Fill some histograms here
    if (fHistsInitialized)
    {
//        fAdcCntHist[view]->Fill(curBinCnt, 1.);
        fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,double(pedCorVal))), 1.);
        fRmsValHist[view]->Fill(std::min(49.9, double(truncRms)), 1.);
        fRmsValProf[view]->Fill(wire, double(truncRms), 1.);
        fMinMaxValProf[view]->Fill(wire, double(minMax), 1.);
        fPedValProf[view]->Fill(wire, truncMean, 1.);
        fPedValHist[view]->Fill(truncMean, 1.);
    }
    
    
    if (wire / fNumWiresToGroup[view] == fHistsWireGroup[view])
    {
        float  leastNeighborRatio = float(std::min(leftNeighbor,rightNeighbor)) / float(countMP);
        size_t wireIdx            = wire % fNumWiresToGroup[view];
        
        if (skewness > 0. && leastNeighborRatio < 0.7)
        {
            short threshold(6);
            
            RawDigitVector::const_iterator stopChirpItr = std::find_if(rawWaveform.begin(),rawWaveform.end(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
        }
        
        if (fHistsInitialized)
        {
            fMinMaxProfiles[view]->Fill(double(wireIdx+0.5), double(minMax), 1.);
            fSkewnessProfiles[view]->Fill(double(wireIdx+0.5), double(skewness), 1.);
            fModeRatioProfiles[view]->Fill(double(wireIdx+0.5), double(leastNeighborRatio), 1.);
        }
    }
    
    return;
}
void RawDigitCharacterizationAlg::getTruncatedRMS(const RawDigitVector& rawWaveform,
                                                  float&                pedestal,
                                                  float&                truncRms) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(rawWaveform.size());
    
    // Fill the vector
    std::transform(rawWaveform.begin(),rawWaveform.end(),adcLessPedVec.begin(),std::bind2nd(std::minus<short>(),pedestal));
    
    // sort in ascending order so we can truncate the sume
    std::sort(adcLessPedVec.begin(), adcLessPedVec.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    int minNumBins = (1. - fTruncMeanFraction) * rawWaveform.size();
    
    // Get the truncated sum
    truncRms = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.begin() + minNumBins, adcLessPedVec.begin(), 0.);
    truncRms = std::sqrt(std::max(0.,truncRms / double(minNumBins)));
    
    return;
}

void RawDigitCharacterizationAlg::getMeanRmsMinMaxAndPedCor(const RawDigitVector& rawWaveform,
                                                            unsigned int          channel,
                                                            float&                truncMean,
                                                            float&                rmsVal,
                                                            float&                minMax,
                                                            float&                pedCorVal) const
{
    // First simply get the mean and rms...
    getMeanRmsAndMinMax(rawWaveform, truncMean, rmsVal, minMax, fTruncMeanFraction);
    
    // Recover the database version of the pedestal
    float pedestal = fPedestalRetrievalAlg.PedMean(channel);
    
    pedCorVal = truncMean - pedestal;
    
    // Fill some histograms here
    if (fHistsInitialized)
    {
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        
        // Recover plane and wire in the plane
        unsigned int view = wids[0].Plane;
        unsigned int wire = wids[0].Wire;
        
        std::pair<RawDigitVector::const_iterator,RawDigitVector::const_iterator> minMaxPair = std::minmax_element(rawWaveform.begin(),rawWaveform.end());
        
        short maxVal = *minMaxPair.second;
        short minVal = *minMaxPair.first;
        short minMax = std::min(maxVal - minVal,199);
        
//        fAdcCntHist[view]->Fill(curBinCnt, 1.);
        fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,double(truncMean - pedestal))), 1.);
        fRmsValHist[view]->Fill(std::min(49.9, double(rmsVal)), 1.);
        fRmsValProf[view]->Fill(wire, double(rmsVal), 1.);
        fMinMaxValProf[view]->Fill(wire, double(minMax), 1.);
        fPedValProf[view]->Fill(wire, truncMean, 1.);
        fPedValHist[view]->Fill(truncMean, 1.);
    }
    
    // Output a message is there is significant different to the pedestal
    if (abs(truncMean - pedestal) > fMaxPedestalDiff)
    {
        mf::LogInfo("RawDigitCharacterizationAlg") << ">>> Pedestal mismatch, channel: " << channel << ", new value: " << truncMean << ", original: " << pedestal << ", rms: " << rmsVal << std::endl;
    }
    
    return;
}

void RawDigitCharacterizationAlg::getMeanRmsAndMinMax(const RawDigitVector& rawWaveform,
                                                      float&                aveVal,
                                                      float&                rmsVal,
                                                      float&                minMax,
                                                      float                 fracBins) const
{
    // first step is to copy the input waveform into a local work vector - convert to floats
    std::vector<float> locWaveform(rawWaveform.size());
    
    std::copy(rawWaveform.begin(),rawWaveform.end(),locWaveform.begin());
    
    // Note that we are not meant to assume that the incoming waveform has been zero suppressed
    // So first task is to get an initial mean and rms...
    // And it can be tricky getting an initial mean, so to start let's really try to get the most probable value from
    // the input vector
    std::map<int,int> frequencyMap;
    int               mpCount(0);
    int               mpVal(0);
    
    for(const auto& val : rawWaveform)
    {
        frequencyMap[val]++;
        
        if (frequencyMap.at(val) > mpCount)
        {
            mpCount = frequencyMap.at(val);
            mpVal   = val;
        }
    }
    
    if (frequencyMap.size() > 1) minMax = frequencyMap.rbegin()->first - frequencyMap.begin()->first;
    else                         minMax = 0;
    
    // take a weighted average of two neighbor bins
    int meanCnt = mpCount;
    int meanSum = mpVal * mpCount;
    
    for(int idx = -3; idx < 4; idx++)
    {
        if (!idx) continue;
        
        std::map<int,int>::iterator neighborItr = frequencyMap.find(mpVal+idx);
    
        if (neighborItr != frequencyMap.end() && 5 * neighborItr->second > mpCount)
        {
            meanSum += neighborItr->first * neighborItr->second;
            meanCnt += neighborItr->second;
        }
    }
    
    float meanWaveform = float(meanSum) / float(meanCnt);
    
    // Use this to do an initial zero suppression
    std::transform(locWaveform.begin(),locWaveform.end(),locWaveform.begin(),std::bind2nd(std::minus<float>(),meanWaveform));
    
    // sort in ascending order so we can truncate before counting any possible signal
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    // And now the local rms calculation...
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.end(), locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size())));
    
    // We use this local rms calculation to try to limit the range of the waveform we use to calculate the truncated mean and rms
    float threshold = std::min(4.,1.5 * localRMS);
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int((1.-fracBins) * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Ok, now we calculate the mean for this range
    float aveSum = std::accumulate(locWaveform.begin(), locWaveform.begin() + minNumBins, 0.);
    aveVal       = aveSum / minNumBins;
    
    std::transform(locWaveform.begin(),locWaveform.begin() + minNumBins,locWaveform.begin(), std::bind2nd(std::minus<float>(),aveVal));
    
    rmsVal = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    rmsVal = std::sqrt(std::max(float(0.),rmsVal / float(minNumBins)));
    
    // Finally adjust the baseline
    aveVal += meanWaveform;
    
    return;
}

bool RawDigitCharacterizationAlg::classifyRawDigitVec(RawDigitVector&         rawWaveform,
                                                      unsigned int            planeIdx,
                                                      unsigned int            wire,
                                                      float                   truncRms,
                                                      short                   minMax,
                                                      short                   mean,
                                                      float                   skewness,
                                                      float                   neighborRatio,
                                                      GroupToDigitIdxPairMap& groupToDigitIdxPairMap) const
{
    // This simply classifies the input waveform:
    // a) determines if it should be added to the list of waveforms to process
    // b) if to be analyzed, places in the group of wires to process
    bool classified(false);
    
    // Dereference the selection/rejection cut
    float selectionCut = fMinMaxSelectionCut[planeIdx];
    float rejectionCut = fRmsRejectionCutHi[planeIdx];
    
    // Selection to process
    if (minMax > selectionCut && truncRms < rejectionCut)
    {
        size_t group = fChannelGroups.channelGroup(planeIdx,wire);
        
        if (groupToDigitIdxPairMap.find(group) == groupToDigitIdxPairMap.end())
            groupToDigitIdxPairMap.insert(std::pair<size_t,RawDigitAdcIdxPair>(group,RawDigitAdcIdxPair()));
        
        groupToDigitIdxPairMap.at(group).first.insert(WireToRawDigitVecPair(wire,rawWaveform));
        groupToDigitIdxPairMap.at(group).second.insert(std::pair<size_t,RawDigitVectorIdxPair>(wire,RawDigitVectorIdxPair(0,rawWaveform.size())));
        
        // Look for chirping wire sections. Confine this to only the V plane
        if (planeIdx == 1)
        {
            // Do wire shape corrections to look for chirping wires and other oddities to avoid
            // Recover our objects...
            WireToAdcIdxMap& wireToAdcIdxMap = groupToDigitIdxPairMap.at(group).second;

            // Set a threshold
            short threshold(6);
            
            // If going from quiescent to on again, then the min/max will be large
            //if (skewnessWireVec[wireIdx] > 0. && minMaxWireVec[wireIdx] > 50 && truncRmsWireVec[wireIdx] > 2.)
            if (skewness > 0. && neighborRatio < 0.7 && minMax > 50)
            {
                RawDigitVector::iterator stopChirpItr = std::find_if(rawWaveform.begin(),rawWaveform.end(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
                
                size_t threshIndex = std::distance(rawWaveform.begin(),stopChirpItr);
                
                if (threshIndex > 60) wireToAdcIdxMap[wire].first = threshIndex;
            }
            // Check in the reverse direction?
            else if (minMax > 20 && neighborRatio < 0.7)
            {
                threshold = 3;
                
                RawDigitVector::reverse_iterator startChirpItr = std::find_if(rawWaveform.rbegin(),rawWaveform.rend(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
                
                size_t threshIndex = std::distance(rawWaveform.rbegin(),startChirpItr);
                
                if (threshIndex > 60) wireToAdcIdxMap[wire].second = rawWaveform.size() - threshIndex;
            }
        }
        
        classified = true;
    }
    
    return classified;
}

template<class T> T RawDigitCharacterizationAlg::getMedian(std::vector<T>& valuesVec, T defaultValue) const
{
    T medianValue(defaultValue);
    
    if (!valuesVec.empty())
    {
        std::sort(valuesVec.begin(),valuesVec.end());
        
        size_t medianIdx = valuesVec.size() / 2;
        
        medianValue = valuesVec[medianIdx];
        
        if (valuesVec.size() > 1 && medianIdx % 2) medianValue = (medianValue + valuesVec[medianIdx+1]) / 2;
    }
    
    return std::max(medianValue,defaultValue);
}
    
}
