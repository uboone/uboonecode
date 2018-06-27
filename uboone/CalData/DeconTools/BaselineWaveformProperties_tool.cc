////////////////////////////////////////////////////////////////////////
/// \file   Baseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/CalData/DeconTools/IBaseline.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "uboone/CalData/DeconTools/WaveformPropertiesAlg.h"

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{

class BaselineWaveformProperties : public IBaseline
{
public:
    explicit BaselineWaveformProperties(const fhicl::ParameterSet& pset);
    
    ~BaselineWaveformProperties();
    
    void configure(const fhicl::ParameterSet& pset)                                override;
    void outputHistograms(art::TFileDirectory&)                              const override;
    
    float GetBaseline(std::vector<float> const&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
    mutable util::WaveformPropertiesAlg<float> fROIPropertiesAlg;
};
    
//----------------------------------------------------------------------
// Constructor.
BaselineWaveformProperties::BaselineWaveformProperties(const fhicl::ParameterSet& pset) :
    fROIPropertiesAlg(pset.get<fhicl::ParameterSet>("ROIPropertiesAlg"))
{
    configure(pset);
}
    
BaselineWaveformProperties::~BaselineWaveformProperties()
{
}
    
void BaselineWaveformProperties::configure(const fhicl::ParameterSet& pset)
{
    return;
}

    
float BaselineWaveformProperties::GetBaseline(std::vector<float> const& holder,
                                    raw::ChannelID_t    channel,
                                    size_t              roiStart,
                                    size_t              roiLen) const
{
    float base = fROIPropertiesAlg.GetWaveformPedestal(holder);
    
    return base;
}
    
void BaselineWaveformProperties::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "BaselinePlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fBaselineVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "BaselinePlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "Baseline;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fBaselineVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BaselineWaveformProperties)
}
