/**
  Eats up CRTData and spits out BernZMQFragments

**/

#ifndef CRTDaqSim_HH_
#define CRTDaqSim_HH_

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "fhiclcpp/ParameterSet.h"

#include <string>
#include <vector>
#include <array>

namespace art{
  class Event;
  class SubRun;
  class Run;
}

// TODO: either use include or forward def
namespace bernfebdaq
{
  class BernZMQFragmentMetadata;
  class BernZMQEvent;
}

#define N_CRT_FEBS 76

namespace crt
{
  class CRTDaqSim :  public art:: EDProducer
  {

    /// Configuration parameters
    //  Name of the producer whose products DaqSim will consume
    std::string fProducerName;
    /// To Make the daq and sim identical
    std::string fInstanceName;
    // The time window in which the boards will be polled in s.
    uint32_t fPollingTime;
    // The static time correction diff in s?
    uint32_t fTimeCorrectionDiff;
    // The static time offset in s?
    uint32_t fTimeOffset;
    // The inter-channel coincidence window time in ns
    uint32_t fCoincidenceTime;

    /// State machine objects
    std::array<bernfebdaq::BernZMQFragmentMetadata, N_CRT_FEBS> fMetadata;
    std::array<std::vector<bernfebdaq::BernZMQEvent>, N_CRT_FEBS> fEvents;

    uint32_t fSeqID;

  public:

    /// Default ctor
    CRTDaqSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDaqSim();

    /// art::EDProducer::produce implementation
    void produce (art::Event&);

    void reconfigure(fhicl::ParameterSet const&);
    void beginRun(art::Run &);
    void beginSubRun(art::SubRun &);
    void endRun(art::Run &);
    void endSubRun(art::SubRun &);

  };
}


#endif  //CRTDaqSim_HH_
