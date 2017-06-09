////////////////////////////////////////////////////////////////////////
/// \file   BaselineStandard.h
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IBaseline.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"

namespace uboone_tool
{

class BaselineStandard : public IBaseline
{
public:
    explicit BaselineStandard(const fhicl::ParameterSet& pset);
    
    ~BaselineStandard();
    
    void configure(const fhicl::ParameterSet& pset)                                override;
    void outputHistograms(art::TFileDirectory&)                              const override;
    
    float GetBaseline(std::vector<float> const&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
    // fhicl parameters
    int    fNumBinsToAverage;
    
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> fSignalShaping;
};

}
