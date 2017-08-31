
#include <cmath>
#include <algorithm>
#include <vector>

#include "RawDigitCorrelatedCorrectionAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cmath>
#include <algorithm>

#include "TVirtualFFT.h"

namespace caldata
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitCorrelatedCorrectionAlg::RawDigitCorrelatedCorrectionAlg(fhicl::ParameterSet const & pset) :
    fFFTAlg(pset)
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitCorrelatedCorrectionAlg") << "RawDigitCorrelatedCorrectionAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitCorrelatedCorrectionAlg::~RawDigitCorrelatedCorrectionAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitCorrelatedCorrectionAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fApplyCorSmoothing     = pset.get<bool>               ("ApplyCorSmoothing",                                        true);
    fApplyFFTCorrection    = pset.get<bool>               ("ApplyFFTCorrection",                                       true);
    fFillFFTHistograms     = pset.get<bool>               ("FillFFTHistograms",                                       false);
    fFFTHistsWireGroup     = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fFFTNumHists           = pset.get<std::vector<size_t>>("FFTNumWaveHistograms",       std::vector<size_t>() = {10,48,48});
    fFFTHistsStartTick     = pset.get<std::vector<double>>("FFTWaveHistsStartTick", std::vector<double>() = {96.,96.,7670.});
    fFFTMinPowerThreshold  = pset.get<std::vector<double>>("FFTPowerThreshold",     std::vector<double>() = {100.,75.,500.});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                          false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                    false);
    fNumRmsToSmoothVec     = pset.get<std::vector<float>> ("NumRmsToSmooth",          std::vector<float>() = {3.6, 3.6, 4.});
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitCorrelatedCorrectionAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    // Following to determine min/max frequencies
    double sampleRate  = fDetectorProperties->SamplingRate();
    double readOutSize = fDetectorProperties->ReadOutWindowSize();
    double maxFreq     = 1000000. / (2. * sampleRate);
    double minFreq     = 1000000. / (2. * sampleRate * readOutSize);
    int    numSamples  = (readOutSize / 2 + 1) / 4;
    
    // quick aside
//    std::cout << "++> plane 0 offset: " << fDetectorProperties->GetXTicksOffset(0,0,0) << std::endl;
//    std::cout << "++> plane 1 offset: " << fDetectorProperties->GetXTicksOffset(1,0,0) << std::endl;
//    std::cout << "++> plane 2 offset: " << fDetectorProperties->GetXTicksOffset(2,0,0) << std::endl;
    
    fFFTHist[0]       = tfs->make<TProfile>("FFTPlaneU",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[1]       = tfs->make<TProfile>("FFTPlaneV",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[2]       = tfs->make<TProfile>("FFTPlaneW",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistLow[0]    = tfs->make<TProfile>("FFTLowPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[1]    = tfs->make<TProfile>("FFTLowPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[2]    = tfs->make<TProfile>("FFTLowPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistCor[0]    = tfs->make<TProfile>("FFTCorPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[1]    = tfs->make<TProfile>("FFTCorPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[2]    = tfs->make<TProfile>("FFTCorPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);

    fCorValHist[0]    = tfs->make<TProfile>("FFTCorValU",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fCorValHist[1]    = tfs->make<TProfile>("FFTCorValV",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fCorValHist[2]    = tfs->make<TProfile>("FFTCorValW",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    
    fFFTCorValHist[0] = tfs->make<TProfile>("FFTHistValU",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    fFFTCorValHist[1] = tfs->make<TProfile>("FFTHistValV",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    fFFTCorValHist[2] = tfs->make<TProfile>("FFTHistValW",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    
    fFFTCorHist[0]    = tfs->make<TProfile>("FFTCorHistU",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fFFTCorHist[1]    = tfs->make<TProfile>("FFTCorHistV",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fFFTCorHist[2]    = tfs->make<TProfile>("FFTCorHistW",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    
    if (fFillFFTHistograms)
    {
        fInputWaveHists.resize(3);
        fOutputWaveHists.resize(3);
        fStdCorWaveHists.resize(3);
        fWaveformProfHists.resize(3);
        fWaveCorHists.resize(3);
    
        for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
        {
            fInputWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fOutputWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fStdCorWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fWaveformProfHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);

            for(size_t histIdx = 0; histIdx < fNumWiresToGroup[viewIdx]; histIdx++)
            {
                std::string histName = "InputWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fInputWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Raw Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
            
                histName = "OutputWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fOutputWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "FFT Corrected Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
            
                histName = "StdCorWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fStdCorWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Standard Corrected Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
                
                histName = "InputWaveform_prof_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
                
                fWaveformProfHists[viewIdx][histIdx] = tfs->make<TH1D>(histName.c_str(), "Waveform Profile;dADC;Count", 400, -200., 200.);
            }
            
            fWaveCorHists[viewIdx].resize(fFFTNumHists[viewIdx]);
            
            for(size_t histIdx = 0; histIdx < fFFTNumHists[viewIdx]; histIdx++)
            {
                std::string histName = "WaveCor_" + std::to_string(viewIdx) + "_tick_" + std::to_string(fFFTHistsStartTick[viewIdx]+histIdx);
                
                fWaveCorHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Tick;ADC", fNumWiresToGroup[viewIdx], 0., fNumWiresToGroup[viewIdx], -500., 500.);
            }
        }
    }
    
    fFFTvsMBProf[0]   = tfs->make<TProfile2D>("FFTvsMBPlaneU", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[1]   = tfs->make<TProfile2D>("FFTvsMBPlaneV", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[2]   = tfs->make<TProfile2D>("FFTvsMBPlaneW", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 3456/16, 0., 3456/16);

    fCorrectionHistVec.resize(3);
    fCorrFixedHistVec.resize(3);
    fErosionHistVec.resize(3);
    fDilationHistVec.resize(3);
    fAverageHistVec.resize(3);
    fDiffHistVec.resize(3);

    fCorMaxHists = 50;
    
    for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
    {
        for(int wireKey = fNumWiresToGroup.at(planeIdx) - 1; wireKey < int(fCorMaxHists * fNumWiresToGroup.at(planeIdx)); wireKey += fNumWiresToGroup.at(planeIdx))
        {
            std::string histName = std::to_string(planeIdx) + "_" + std::to_string(wireKey);
            
            fCorrectionHistVec.at(planeIdx)[wireKey] = tfs->make<TH1D>(("Correction_" + histName).c_str(), ";Tick", 6400, 0, 6400);
            fCorrFixedHistVec.at(planeIdx)[wireKey]  = tfs->make<TH1D>(("CorrFixed_"  + histName).c_str(), ";Tick", 6400, 0, 6400);
            fErosionHistVec.at(planeIdx)[wireKey]    = tfs->make<TH1D>(("Erosion_"    + histName).c_str(), ";Tick", 6400, 0, 6400);
            fDilationHistVec.at(planeIdx)[wireKey]   = tfs->make<TH1D>(("Dilation_"   + histName).c_str(), ";Tick", 6400, 0, 6400);
            fAverageHistVec.at(planeIdx)[wireKey]    = tfs->make<TH1D>(("Average_"    + histName).c_str(), ";Tick", 6400, 0, 6400);
            fDiffHistVec.at(planeIdx)[wireKey]       = tfs->make<TH1D>(("Difference_" + histName).c_str(), ";Tick", 6400, 0, 6400);
        }
    }
 
}

void RawDigitCorrelatedCorrectionAlg::smoothCorrectionVec(std::vector<float>& corValVec, raw::ChannelID_t channel) const
{
    // Try working (again) with erosion/dilation vectors
    std::vector<float> erosionVec;
    std::vector<float> dilationVec;
    std::vector<float> averageVec;
    std::vector<float> differenceVec;
    
    getErosionDilationAverageDifference(corValVec,
                                        erosionVec,
                                        dilationVec,
                                        averageVec,
                                        differenceVec);
    
    // Get the mean and rms of the difference vector
    float diffMean;
    float diffRMS;
    
    getTruncatedMeanRMS(differenceVec, diffMean, diffRMS);
    
    // set a low bar...
    float diffThreshold = diffMean + 5. * diffRMS;

    // Recover plane/wire from channel
    std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
    
    unsigned int plane = wids[0].Plane;
    unsigned int wire  = wids[0].Wire;

    // Fill histograms?
    std::map<int,TH1D*>::const_iterator corHistMapItr = fCorrectionHistVec.at(plane).find(wire);
    
    if (corHistMapItr != fCorrectionHistVec.at(plane).end())
    {
        std::cout << "--> Plane: " << plane << ", wire: " << wire << ", diff mean: " << diffMean << ", rms: " << diffRMS << ", threshold: " << diffThreshold << std::endl;
        
        TH1D* corHist        = fCorrectionHistVec.at(plane).at(wire);
        TH1D* erosionHist    = fErosionHistVec.at(plane).at(wire);
        TH1D* dilationHist   = fDilationHistVec.at(plane).at(wire);
        TH1D* averageHist    = fAverageHistVec.at(plane).at(wire);
        TH1D* differenceHist = fDiffHistVec.at(plane).at(wire);
        
        for(size_t tickIdx = 0; tickIdx < corValVec.size(); tickIdx++)
        {
            corHist->Fill(tickIdx,corValVec.at(tickIdx));
            erosionHist->Fill(tickIdx,erosionVec.at(tickIdx));
            dilationHist->Fill(tickIdx,dilationVec.at(tickIdx));
            averageHist->Fill(tickIdx,averageVec.at(tickIdx));
            differenceHist->Fill(tickIdx,differenceVec.at(tickIdx));
        }
    }

    size_t aveIdx(0);
    
    while(++aveIdx < averageVec.size())
    {
        // Have we gone over threshold?
        if (differenceVec.at(aveIdx) > diffThreshold)
        {
            size_t lastIdx = aveIdx;
            
            // We are over threshold, back up to find the start point
//            while(lastIdx > 0 && differenceVec.at(lastIdx) > diffMean + 3. * diffRMS) lastIdx--;
            if (lastIdx > 0) lastIdx--;
            
            size_t nextIdx = aveIdx;
            
            // Look for the end point which will be when we return to baseline
//            while(nextIdx < corValVec.size() && differenceVec.at(nextIdx) > diffMean + diffRMS) nextIdx++;
            while(nextIdx < corValVec.size() && differenceVec.at(nextIdx) > diffMean + 2.*diffRMS) nextIdx++;
            
            if (nextIdx >= corValVec.size())
            {
                nextIdx = corValVec.size() - 1;
            }
            
            if (nextIdx > lastIdx+1)
            {
                float lastVal = corValVec.at(lastIdx);
                float nextVal = corValVec.at(nextIdx);
                float slope   = (nextVal - lastVal) / float(nextIdx - lastIdx);
                
                if (nextIdx >= corValVec.size() - 1) slope = 0.;
                
                for(size_t idx = lastIdx; idx < nextIdx; idx++)
                    corValVec.at(idx) = slope * float(idx - lastIdx) + lastVal;
                
                aveIdx = nextIdx;
            }
            
        }
    }

    std::map<int,TH1D*>::const_iterator corrFixedMapItr = fCorrFixedHistVec.at(plane).find(wire);
    
    if (corrFixedMapItr != fCorrFixedHistVec.at(plane).end())
    {
        
        TH1D* corFixed = corrFixedMapItr->second;
        
        for(size_t tickIdx = 0; tickIdx < corValVec.size(); tickIdx++)
        {
            corFixed->Fill(tickIdx,corValVec.at(tickIdx));
        }
    }

    return;
}

void RawDigitCorrelatedCorrectionAlg::removeCorrelatedNoise(RawDigitAdcIdxPair& digitIdxPair,
                                                            raw::ChannelID_t    channel,
                                                            std::vector<float>& truncMeanWireVec,
                                                            std::vector<float>& truncRmsWireVec,
                                                            std::vector<short>& minMaxWireVec,
                                                            std::vector<short>& meanWireVec,
                                                            std::vector<float>& skewnessWireVec,
                                                            std::vector<float>& neighborRatioWireVec,
                                                            std::vector<float>& pedCorWireVec) const
{
    // This method represents an enhanced implementation of "Corey's Algorithm" for correcting the
    // correlated noise across a group of wires. The primary enhancement involves using a FFT to
    // "fit" for the underlying noise as a way to reduce the impact on the signal.
    WireToRawDigitVecMap& wireToRawDigitVecMap = digitIdxPair.first;
    WireToAdcIdxMap&      wireToAdcIdxMap      = digitIdxPair.second;
    
    // Recover plane/wire from channel
    std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);

    unsigned int plane = wids[0].Plane;
//    unsigned int wire  = wids[0].Wire;
    
    size_t maxTimeSamples(wireToRawDigitVecMap.begin()->second.size());
    size_t baseWireIdx(wireToRawDigitVecMap.begin()->first - wireToRawDigitVecMap.begin()->first % fNumWiresToGroup[plane]);
    
    // Should we turn on our diagnostics blocks?
    bool doFFTCorrection(fFillFFTHistograms && baseWireIdx / fNumWiresToGroup[plane] == fFFTHistsWireGroup[plane]);
    
    // Keep track of the "waveform" for the correction
    std::vector<float> origCorValVec;

    // Don't try to do correction if too few wires unless they have gaps
    if (wireToAdcIdxMap.size() > 2) // || largestGapSize > 2)
    {
        std::vector<float> corValVec;
        
        corValVec.resize(maxTimeSamples, 0.);

        // Build the vector of corrections for each time bin
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Try to do the calculations for histogramming outside of the loop over channels
            bool   fillHists(false);
            size_t histIdx(0);

            // Diagnostics
            if (doFFTCorrection)
            {
                fillHists = sampleIdx >= fFFTHistsStartTick[plane] && sampleIdx < fFFTHistsStartTick[plane] + fFFTNumHists[plane];
                histIdx   = sampleIdx - fFFTHistsStartTick[plane];
            }
            
            // Define a vector for accumulating values...
            // Loop over the wires at this time bin and get their pedestal corrected ADC values
            // We'll use a simple stl vector for this
            std::vector<float> adcValuesVec;
            
            for(const auto& wireAdcItr : wireToAdcIdxMap)
            {
                // Check that we should be doing something in this range
                // Note that if the wire is not to be considered then the "start" bin will be after the last bin
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second) continue;
                
                int wireIdx(wireAdcItr.first - baseWireIdx);

                // Accumulate
                adcValuesVec.push_back(float(wireToRawDigitVecMap.at(wireAdcItr.first)[sampleIdx]) - truncMeanWireVec[wireIdx]);
                
                // Make hists if requested
                if (fillHists)
                {
                    fWaveCorHists[plane][histIdx]->Fill(wireIdx, adcValuesVec.back());
                }
            }
            
            float medianValue = getMedian(adcValuesVec, float(-10000.));
            //float medianValue = getMostProbable(adcValuesVec, float(-10000.));
            
            corValVec[sampleIdx] = medianValue;
        }
        
        // Try to eliminate any real outliers
        if (fApplyCorSmoothing) smoothCorrectionVec(corValVec, channel);
        
        // Diagnostics block
        if (doFFTCorrection)
        {
            origCorValVec = corValVec;
            
            for(size_t tick = 0; tick < maxTimeSamples; tick++)
            {
                double corVal = origCorValVec[tick];
                
                for(const auto& wireAdcItr : wireToRawDigitVecMap)
                {
                    short  inputWaveformVal = wireAdcItr.second[tick];
                    size_t wireIdx          = wireAdcItr.first % fNumWiresToGroup[plane];
                    double inputVal         = inputWaveformVal - truncMeanWireVec[plane];
                    double outputVal        = std::round(inputVal - corVal);
                    
                    fInputWaveHists[plane][wireIdx]->Fill(tick, inputVal, 1.);
                    fStdCorWaveHists[plane][wireIdx]->Fill(tick, outputVal, 1.);
                    fWaveformProfHists[plane][wireIdx]->Fill(std::round(inputVal) + 0.5, 1.);
                }
            }
        }
        
        // Get the FFT correction
        if (fApplyFFTCorrection) fFFTAlg.getFFTCorrection(corValVec,fFFTMinPowerThreshold[plane]);
        
        // Now go through and apply the correction
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Now run through and apply correction
            for (const auto& wireAdcItr : wireToAdcIdxMap)
            {
                float corVal(corValVec[sampleIdx]);
                int   wireIdx(wireAdcItr.first);
                
                // If the "start" bin is after the "stop" bin then we are meant to skip this wire in the averaging process
                // Or if the sample index is in a chirping section then no correction is applied.
                // Both cases are handled by looking at the sampleIdx
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second)
                    corVal = 0.;
                
                short& rawDataTimeVal = wireToRawDigitVecMap.at(wireIdx).at(sampleIdx);
                
                // Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
                float newAdcValueFloat = float(rawDataTimeVal) - corVal - pedCorWireVec[wireIdx - baseWireIdx];
                
                rawDataTimeVal = std::round(newAdcValueFloat);
                
                if (doFFTCorrection)
                {
                    float stdCorWaveVal = std::round(newAdcValueFloat + pedCorWireVec[wireIdx - baseWireIdx] - truncMeanWireVec[wireIdx - baseWireIdx]);
                    
                    fOutputWaveHists[plane][wireIdx - baseWireIdx]->Fill(sampleIdx, stdCorWaveVal, 1.);
                    
                    fFFTCorHist[plane]->Fill(sampleIdx,corVal);
                }
            }
        }
        
        // Try working (again) with erosion/dilation vectors
        for(const auto& wireAdcItr : wireToAdcIdxMap)
        {
            RawDigitVector     erosionVec;
            RawDigitVector     dilationVec;
            std::vector<float> averageVec;
            RawDigitVector     differenceVec;
            
            RawDigitVector& rawDigitVec = wireToRawDigitVecMap.at(wireAdcItr.first);
            
            getErosionDilationAverageDifference(rawDigitVec,
                                                erosionVec,
                                                dilationVec,
                                                averageVec,
                                                differenceVec);
            
            // Get the mean and rms of the difference vector
            float diffMean;
            float diffRMS;
            
            getTruncatedMeanRMS(differenceVec, diffMean, diffRMS);
            
            // set a low bar...
            float diffThreshold = diffMean + 2. * diffRMS;
            
            size_t aveIdx(0);
            
            while(++aveIdx < averageVec.size())
            {
                // Have we gone over threshold?
                if (differenceVec.at(aveIdx) > diffThreshold)
                {
                    size_t lastIdx = aveIdx;
                    
                    // We are over threshold, back up to find the start point
                    while(lastIdx > 0 && differenceVec.at(lastIdx) > diffMean) lastIdx--;
                    
                    size_t nextIdx = aveIdx;
                    
                    // Look for the end point which will be when we return to baseline
                    while(nextIdx < averageVec.size() && differenceVec.at(nextIdx) > diffMean) nextIdx++;
                    
                    if (nextIdx >= averageVec.size())
                    {
                        nextIdx = averageVec.size() - 1;
                    }
                    
                    if (nextIdx > lastIdx)
                    {
                        float lastVal = averageVec.at(lastIdx);
                        float nextVal = averageVec.at(nextIdx);
                        float slope   = (nextVal - lastVal) / float(nextIdx - lastIdx);
                        
                        if (nextIdx >= averageVec.size() - 1) slope = 0.;
                        
                        for(size_t idx = lastIdx; idx < nextIdx; idx++)
                            averageVec.at(idx) = slope * float(idx - lastIdx) + lastVal;
                        
                        aveIdx = nextIdx;
                    }
                    
                }
            }
            
            // Now adjust the input vector...
            float truncMean = truncMeanWireVec[wireAdcItr.first - baseWireIdx];
            
            std::transform(rawDigitVec.begin()+wireAdcItr.second.first,
                           rawDigitVec.begin()+wireAdcItr.second.second,
                           averageVec.begin()+wireAdcItr.second.first,
                           rawDigitVec.begin()+wireAdcItr.second.first,
                           [truncMean](const auto& left, const auto& right){return std::round(left - right + truncMean);});
        }
    }
    
    // Final diagnostics block
    if (doFFTCorrection && !origCorValVec.empty())
    {
        double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        int    fftDataSize = origCorValVec.size();
        
        TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
        
        std::vector<double> fftInputArray(fftDataSize,0.);
        
        for(size_t tick = 0; tick < size_t(fftDataSize); tick++)
        {
            fftInputArray[tick] = origCorValVec[tick];
            
            fCorValHist[plane]->Fill(tick,origCorValVec[tick]);
        }
        
        fftr2c->SetPoints(fftInputArray.data());
        fftr2c->Transform();
        
        // Recover the power spectrum...
        std::vector<double> realVals(fftDataSize,0.);
        std::vector<double> imaginaryVals(fftDataSize,0.);
        double realPart(0.);
        double imaginaryPart(0.);
        
        fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
            
        for(size_t idx = 0; idx < size_t(fftDataSize)/2; idx++)
        {
            realPart      = realVals[idx];
            imaginaryPart = imaginaryVals[idx];
                
            double bin   = (idx * sampleFreq) / readOutSize;
            double power = std::sqrt(realPart*realPart + imaginaryPart*imaginaryPart);
            
            fFFTCorValHist[plane]->Fill(bin, power, 1.);
        }
    }
    
    // Run an FFT here to check our "corrected" wires
    if (fRunFFTCorrected)
    {
        double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        
        for(const auto& wireAdcItr : wireToAdcIdxMap)
        {
            if (wireAdcItr.second.first > wireAdcItr.second.second) continue;

            RawDigitVector& rawDataTimeVec = wireToRawDigitVecMap.at(wireAdcItr.first);
            
            size_t wireIdx = wireAdcItr.first  - baseWireIdx;
            
            int fftDataSize = rawDataTimeVec.size();
        
            TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
        
            std::vector<double> fftInputArray(fftDataSize,0.);
        
            for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = rawDataTimeVec[idx] - truncMeanWireVec[wireIdx];
        
            fftr2c->SetPoints(fftInputArray.data());
            fftr2c->Transform();
        
            // Recover the power spectrum...
            double realPart(0.);
            double imaginaryPart(0.);
        
            for(size_t idx = 0; idx < size_t(fftDataSize)/2; idx++)
            {
                fftr2c->GetPointComplex(idx+1, realPart, imaginaryPart);
            
                double bin   = (idx * sampleFreq) / readOutSize;
                double power = realPart*realPart + imaginaryPart*imaginaryPart;
            
                if (power > 0.) power = std::sqrt(power);
            
                fFFTHistCor[plane]->Fill(bin, power, 1.);
            }
        }
    }
    
    return;
}
    
template<class T> T RawDigitCorrelatedCorrectionAlg::getMedian(std::vector<T>& valuesVec, T defaultValue) const
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
    
float RawDigitCorrelatedCorrectionAlg::getMostProbable(const std::vector<float>& valuesVec, float defaultValue) const
{
    float mostProbableValue(defaultValue);
    
    if (!valuesVec.empty())
    {
        std::map<int,int> frequencyMap;
        int               mpCount(0);
        int               mpVal(0);
        
        for(const auto& val : valuesVec)
        {
            int intVal = std::round(2.*val);
            
            frequencyMap[intVal]++;
            
            if (frequencyMap.at(intVal) > mpCount)
            {
                mpCount = frequencyMap.at(intVal);
                mpVal   = intVal;
            }
        }
        
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
        
        mostProbableValue = 0.5 * float(meanSum) / float(meanCnt);
    }
    
    return std::max(mostProbableValue,defaultValue);
}
    
void RawDigitCorrelatedCorrectionAlg::getErosionDilationAverageDifference(const RawDigitVector& inputWaveform,
                                                                          RawDigitVector&       erosionVec,
                                                                          RawDigitVector&       dilationVec,
                                                                          std::vector<float>&   averageVec,
                                                                          RawDigitVector&       differenceVec) const
{
    // Set the window size
    int halfWindowSize(30);
    
    // Define the "smallest" function
    //    auto smaller = [](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);};
    
    // Initialize min and max elements
    std::pair<RawDigitVector::const_iterator,RawDigitVector::const_iterator> minMaxItr = std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);
    
    RawDigitVector::const_iterator minElementItr = minMaxItr.first;
    RawDigitVector::const_iterator maxElementItr = minMaxItr.second;
    
    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    averageVec.resize(inputWaveform.size());
    differenceVec.resize(inputWaveform.size());
    
    // Now loop through remaining elements and complete the vectors
    RawDigitVector::iterator     minItr = erosionVec.begin();
    RawDigitVector::iterator     maxItr = dilationVec.begin();
    std::vector<float>::iterator aveItr = averageVec.begin();
    RawDigitVector::iterator     difItr = differenceVec.begin();
    
    for (RawDigitVector::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
        *aveItr++ = 0.5 * (*maxElementItr + *minElementItr);
        *difItr++ = *maxElementItr - *minElementItr;
    }
    
    return;
}
    
void RawDigitCorrelatedCorrectionAlg::getTruncatedMeanRMS(const RawDigitVector& waveform, float& mean, float& rms) const
{
    // We need to get a reliable estimate of the mean and can't assume the input waveform will be ~zero mean...
    // The input waveform will in in (short) int form, so one can employ a map method to find the
    // most probable value, then develop mean about that
    std::map<int,int> frequencyMap;
    int               mpCount(0);
    int               mpVal(0);
    
    for(const auto& val : waveform)
    {
        frequencyMap[val]++;
        
        if (frequencyMap.at(val) > mpCount)
        {
            mpCount = frequencyMap.at(val);
            mpVal   = val;
        }
    }
    
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
    
    mean = float(meanSum) / float(meanCnt);
    
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform(waveform.size());
    
    // convert input waveform to real and zero suppress all at once...
    std::transform(waveform.begin(), waveform.end(), locWaveform.begin(),std::bind2nd(std::minus<float>(),mean));
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + locWaveform.size()/2, locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size()/2)));
    
    float threshold = std::min(4.,3. * localRMS);
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(0.5 * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // recalculate the rms
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    rms = std::sqrt(std::max(float(0.),localRMS / float(minNumBins)));
    
    return;
}
    
void RawDigitCorrelatedCorrectionAlg::getErosionDilationAverageDifference(const std::vector<float>& inputWaveform,
                                                                          std::vector<float>&       erosionVec,
                                                                          std::vector<float>&       dilationVec,
                                                                          std::vector<float>&       averageVec,
                                                                          std::vector<float>&       differenceVec) const
{
    // Set the window size
    int halfWindowSize(10);
    
    // Define the "smallest" function
    //    auto smaller = [](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);};
    
    // Initialize min and max elements
    std::pair<std::vector<float>::const_iterator,std::vector<float>::const_iterator> minMaxItr = std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);
    
    std::vector<float>::const_iterator minElementItr = minMaxItr.first;
    std::vector<float>::const_iterator maxElementItr = minMaxItr.second;
    
    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    averageVec.resize(inputWaveform.size());
    differenceVec.resize(inputWaveform.size());
    
    // Now loop through remaining elements and complete the vectors
    std::vector<float>::iterator minItr = erosionVec.begin();
    std::vector<float>::iterator maxItr = dilationVec.begin();
    std::vector<float>::iterator aveItr = averageVec.begin();
    std::vector<float>::iterator difItr = differenceVec.begin();
    
    for (std::vector<float>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
        *aveItr++ = 0.5 * (*maxElementItr + *minElementItr);
        *difItr++ = *maxElementItr - *minElementItr;
    }
    
    return;
}
    
void RawDigitCorrelatedCorrectionAlg::getTruncatedMeanRMS(const std::vector<float>& waveform, float& mean, float& rms) const
{
    // We need to get a reliable estimate of the mean and can't assume the input waveform will be ~zero mean...
    // The input waveform will in in (short) int form, so one can employ a map method to find the
    // most probable value, then develop mean about that
    std::map<int,int> frequencyMap;
    int               mpCount(0);
    int               mpVal(0);
    
    for(const auto& val : waveform)
    {
        int intVal = std::round(2 * val);
        
        frequencyMap[intVal]++;
        
        if (frequencyMap.at(intVal) > mpCount)
        {
            mpCount = frequencyMap.at(intVal);
            mpVal   = intVal;
        }
    }
    
    // take a weighted average of two neighbor bins
    int meanCnt = 0;
    int meanSum = 0;
    
    for(int idx = -3; idx < 4; idx++)
    {
        std::map<int,int>::iterator neighborItr = frequencyMap.find(mpVal+idx);
        
        if (neighborItr != frequencyMap.end() && 5 * neighborItr->second > mpCount)
        {
            meanSum += neighborItr->first * neighborItr->second;
            meanCnt += neighborItr->second;
        }
    }
    
    mean = 0.5 * float(meanSum) / float(meanCnt);
    
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform(waveform.size());
    
    // convert input waveform to real and zero suppress all at once...
    std::transform(waveform.begin(), waveform.end(), locWaveform.begin(),std::bind2nd(std::minus<float>(),mean));
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + locWaveform.size()/2, locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size()/2)));
    
    float threshold = 6. * localRMS;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(0.6 * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // recalculate the rms
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    rms = std::sqrt(std::max(float(std::fabs(0.1*mean)),localRMS / float(minNumBins)));
    
    return;
}

}
