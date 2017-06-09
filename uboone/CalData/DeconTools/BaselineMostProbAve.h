////////////////////////////////////////////////////////////////////////
/// \file   Baseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IBaseline.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"

namespace uboone_tool
{

class BaselineMostProbAve : public IBaseline
{
public:
    explicit BaselineMostProbAve(const fhicl::ParameterSet& pset);
    
    ~BaselineMostProbAve();
    
    void configure(const fhicl::ParameterSet& pset)                                override;
    void outputHistograms(art::TFileDirectory&)                              const override;
    
    float GetBaseline(std::vector<float> const&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> fSignalShaping;
};
    
}
