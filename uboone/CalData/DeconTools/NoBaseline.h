////////////////////////////////////////////////////////////////////////
/// \file   NoBaseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IBaseline.h"

namespace uboone_tool
{

class NoBaseline : public IBaseline
{
public:
    explicit NoBaseline(const fhicl::ParameterSet& pset);
    
    ~NoBaseline();
    
    void configure(const fhicl::ParameterSet& pset)                                override;
    void outputHistograms(art::TFileDirectory&)                              const override;
    
    float GetBaseline(std::vector<float> const&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
};
    
}
