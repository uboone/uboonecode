////////////////////////////////////////////////////////////////////////
/// \file   Baseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IBaseline.h"
#include "uboone/CalData/DeconTools/WaveformPropertiesAlg.h"

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
    
}
