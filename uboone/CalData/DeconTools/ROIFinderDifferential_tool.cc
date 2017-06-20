////////////////////////////////////////////////////////////////////////
/// \file   ROIFinderDifferential.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/CalData/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{

class ROIFinderDifferential : public IROIFinder
{
public:
    explicit ROIFinderDifferential(const fhicl::ParameterSet& pset);
    
    ~ROIFinderDifferential();
    
    void configure(const fhicl::ParameterSet& pset)                          override;
    void outputHistograms(art::TFileDirectory&)                        const override;
    
    void FindROIs(const Waveform&, size_t, double, CandidateROIVec&)   const override;
    
private:
    // The actual ROI finding algorithm
    void findROICandidates(Waveform::const_iterator       startItr,
                           Waveform::const_iterator       stopItr,
                           const geo::PlaneID::PlaneID_t& plane,
                           size_t                         roiStartTick,
                           float                          roiThreshold,
                           CandidateROIVec&               roiCandidateVec) const;
    
    // Member variables from the fhicl file
    std::vector<float>          fNumSigma;                   ///< "# sigma" rms noise for ROI threshold
    std::vector<float>          fZeroPoint;                  ///< Tolerance for return to zero
    std::vector<size_t>         fMax2MinDistance;            ///< Maxmimum allow peak to peak distance
    std::vector<unsigned short> fPreROIPad;                  ///< ROI padding
    std::vector<unsigned short> fPostROIPad;                 ///< ROI padding
    
    // Services
    const geo::GeometryCore*    fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFinderDifferential::ROIFinderDifferential(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIFinderDifferential::~ROIFinderDifferential()
{
}
    
void ROIFinderDifferential::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    std::vector<unsigned short> uin;
    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fNumSigma        = pset.get< std::vector<float> >         ("NumSigma"       );
    fZeroPoint       = pset.get< std::vector<float> >         ("ZeroPoint"      );
    fMax2MinDistance = pset.get< std::vector<size_t> >        ("Max2MinDistance");
    uin              = pset.get< std::vector<unsigned short> >("uPlaneROIPad"   );
    vin              = pset.get< std::vector<unsigned short> >("vPlaneROIPad"   );
    zin              = pset.get< std::vector<unsigned short> >("zPlaneROIPad"   );
    
    if(uin.size() != 2 || vin.size() != 2 || zin.size() != 2) {
        throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }
    
    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = uin[0];
    fPostROIPad[0] = uin[1];
    fPreROIPad[1]  = vin[0];
    fPostROIPad[1] = vin[1];
    fPreROIPad[2]  = zin[0];
    fPostROIPad[2] = zin[1];
    
    return;
}

void ROIFinderDifferential::FindROIs(const Waveform& waveform, size_t channel, double rmsNoise, CandidateROIVec& roiVec) const
{
    // The idea here is to consider the input waveform - if an induction plane then it is already in differential form,
    // if a collection plane then we form the differential - then smooth and look for peaks. The technique will be to
    // consider the peak signature as a maximum followed some bins later by a mininum and those whose difference between
    // max and min is more than the threshold are kept.
    
    // First up, determine what kind of wire we have
    std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
    const geo::PlaneID&      planeID = wids[0].planeID();
    geo::SigType_t           sigType = fGeometry->SignalType(planeID);
    
    // Local copy of the input waveform
    Waveform localWaveform(waveform.size());
    
    // If we have a collection plane then take the derivative
    if (sigType == geo::kCollection)
    {
        localWaveform[0]                 = 0;
        localWaveform[waveform.size()-1] = 0;
        
        for(size_t idx = 1; idx < waveform.size()-1; idx++)
            localWaveform[idx] = 0.5 * (waveform.at(idx+1) - waveform.at(idx-1));
    }
    // Otherwise a straight copy
    else std::copy(waveform.begin(),waveform.end(),localWaveform.begin());
    
    // Now smooth the differential waveform
    // This version for testing, we'll need a better version for production
    // We don't want the smoothing procedure to factor into the output
    // So start by making a local copy of the input vector
    Waveform tempVec = localWaveform;
    
    // Now run the "triangle" smoothing operation
    for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
        localWaveform[idx] = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;
    
    // At this point ready to search for candidate ROI's...
    // Idea will be to follow technique similar to hit finding which recursively searches for the next largest peak in the waveform.
    float roiThreshold = fNumSigma[planeID.Plane] * rmsNoise;
    
    findROICandidates(localWaveform.begin(),localWaveform.end(),planeID.Plane,0,roiThreshold,roiVec);
    
    // pad the ROIs
    for(auto& roi : roiVec)
    {
        // low ROI end
        roi.first  = std::max(int(roi.first - fPreROIPad[planeID.Plane]),0);
        // high ROI end
        roi.second = std::min(roi.second + fPostROIPad[planeID.Plane], waveform.size() - 1);
    }
    
    // merge overlapping (or touching) ROI's
    if(roiVec.size() > 1)
    {
        // temporary vector for merged ROIs
        CandidateROIVec tempRoiVec;
        
        // Loop through candidate roi's
        size_t startRoi = roiVec.front().first;
        size_t stopRoi  = startRoi;
        
        for(auto& roi : roiVec)
        {
            if (roi.first <= stopRoi) stopRoi = roi.second;
            else
            {
                tempRoiVec.push_back(CandidateROI(startRoi,stopRoi));
                
                startRoi = roi.first;
                stopRoi  = roi.second;
            }
        }
        
        // Make sure to get the last one
        tempRoiVec.push_back(CandidateROI(startRoi,stopRoi));
        
        roiVec = tempRoiVec;
    }
    
    return;
}
    
void ROIFinderDifferential::findROICandidates(Waveform::const_iterator       startItr,
                                              Waveform::const_iterator       stopItr,
                                              const geo::PlaneID::PlaneID_t& plane,
                                              size_t                         roiStartTick,
                                              float                          roiThreshold,
                                              CandidateROIVec&               roiCandidateVec) const
{
    // Require a minimum length
    size_t roiLength = std::distance(startItr,stopItr);
    
    if (roiLength > 2)
    {
        // The idea here is to find the largest positive value in the input waveform and use this as the basis of
        // our search for a candidate hit
        std::pair<Waveform::const_iterator,Waveform::const_iterator> minMaxItr = std::minmax_element(startItr,  stopItr);
        
        // It makes no sense to continue if the minimum is the first bin?
        // This is the bare mininum otherwise we could overrun the end of the waveform
        if (std::distance(minMaxItr.second, stopItr) > 2 && std::distance(startItr, minMaxItr.first) > 2)
        {
            Waveform::const_iterator maxItr = minMaxItr.second;
            Waveform::const_iterator minItr = minMaxItr.first;
            
            // Reset either the max or min iterator depending on which is bigger
            if (*maxItr > std::fabs(*minItr))
            {
                // The maximum is larger so search forward from here for the true minimum
                minItr = maxItr;
                
                float prevValue = *minItr++;
                float lastValue = *minItr++;
                
                while(minItr != stopItr)
                {
                    // Decreasing for two bins
                    if (prevValue < 0. && lastValue > prevValue && *minItr > lastValue)
                    {
                        // reset to what was the actual minimum value
                        minItr -= 2;
                        break;
                    }
                    
                    prevValue = lastValue;
                    lastValue = *minItr++;
                }
            }
            else
            {
                // Otherwise, we are starting at the minimum and searching backwards to find the max
                maxItr = minItr;
                
                float prevValue = *maxItr--;
                float lastValue = *maxItr--;
                
                while(maxItr != startItr)
                {
                    // Decreasing for two bins
                    if (prevValue > 0. && lastValue < prevValue && *maxItr < lastValue)
                    {
                        // reset to what was the actual minimum value
                        maxItr += 2;
                        break;
                    }
                    
                    prevValue = lastValue;
                    lastValue = *maxItr--;
                }
            }
        
            // Check that the range from maximum to minimum is over threshold
            float  maxMinRange    = *maxItr - *minItr;
            size_t maxMinDistance = std::distance(maxItr,minItr);
        
            if (maxMinRange > roiThreshold && maxMinDistance < fMax2MinDistance.at(plane))
            {
                // to complete the edges of the ROI, search both sides for the point which is essentially back to zero
                while(maxItr != startItr)
                {
                    if (std::fabs(*maxItr) < fZeroPoint.at(plane)) break;
                    maxItr--;
                }
            
                while(minItr != stopItr)
                {
                    if (std::fabs(*minItr) < fZeroPoint.at(plane)) break;
                    minItr++;
                }
            
                // Before saving this ROI, look for candidates preceeding this one
                // Note that preceeding snippet will reference to the current roiStartTick
                findROICandidates(startItr, maxItr, plane, roiStartTick, roiThreshold, roiCandidateVec);
            
                // Save this ROI
                size_t newStartTick = std::distance(startItr,maxItr) + roiStartTick;
                size_t newStopTick  = std::distance(startItr,minItr) + roiStartTick;
            
                roiCandidateVec.push_back(CandidateROI(newStartTick, newStopTick));
            
                // Now look for candidate ROI's downstream of this one
                findROICandidates(minItr, stopItr, plane, newStopTick, roiThreshold, roiCandidateVec);
            }
        }
    }
    
    return;
}

    
void ROIFinderDifferential::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "ROIFinderPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fROIFinderVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "ROIFinderPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "ROIFinder;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fROIFinderVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(ROIFinderDifferential)
}
