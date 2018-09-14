#include "uboone/CRT/CRTProducts/CRTSimData.hh"
#include "uboone/CRT/CRTDaqSim.hh"


#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "artdaq-core/Data/Fragments.hh"
#include "bernfebdaq-core/Overlays/FragmentType.hh"


#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Run.h"

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{

  CRTDaqSim::CRTDaqSim(const fhicl::ParameterSet& pSet):
    fProducerName(pSet.get<std::string>("InputTag","crtdetsim")) ,
    fInstanceName(pSet.get<std::string>("InstanceLabel","BernZMQ")),
    fPollingTime(pSet.get<unsigned>("PollingTime",1)),
    fTimeCorrectionDiff(pSet.get<unsigned>("TimeCorrectionDiff",0)),
    fTimeOffset(pSet.get<unsigned>("TimeOffset",0)),
    fCoincidenceTime(pSet.get<unsigned>("CoincidenceTime",0))
  {
    //TODO: add in either some time jitter (POISSON), or adc noise (Gaussian)
    //art::ServiceHandle<rndm::NuRandomService> Seeds;
    //Seeds->createEngine(*this, "", "CRTDaq", pSet, "Seed");
    produces< std::vector<artdaq::Fragment> >(
        fInstanceName
      );
  }

  CRTDaqSim::~CRTDaqSim()
  {

  }

  void CRTDaqSim::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<artdaq::Fragment> > frags(
        new std::vector<artdaq::Fragment>);

    art::Handle<std::vector<crt::CRTSimData> > crtDataHandle;
    evt.getByLabel(fProducerName, crtDataHandle);

    for (auto crtDat : *crtDataHandle)
    {

      //TODO: Check the channel numbering again.
      unsigned feb_n = crtDat.Channel()/32;
      unsigned adc_n = crtDat.Channel()%32;

      //channel0ID = 32 * moduleID + 2 * stripID + 0
      mf::LogInfo("CRTDaqSim")<<"Converting CRTSimData for FEB: "<<feb_n<<" ADC_N: "<<adc_n<<std::endl<<"For input channel number: "<<crtDat.Channel();

      if(feb_n>=N_CRT_FEBS || adc_n>=32)
      {
        mf::LogWarning("CRTDaqSim")<<"Cannot Convert CRTSimData";
        continue;
      }

      bool is_part_of_previous_event = false;

      for(auto event_it = fEvents[feb_n].begin() ; 
          event_it != fEvents[feb_n].end(); ++event_it)
      {
        if ( crtDat.T0() >=event_it->Time_TS0()-fCoincidenceTime && crtDat.T0() <= event_it->Time_TS0()-fCoincidenceTime  ){
          is_part_of_previous_event = true;
          //add the ADC to this one
          mf::LogInfo("CRTDaqSim")<<"Adding data to existing Fragment";
          event_it->adc[adc_n] = crtDat.ADC();
          break;
        }
      }
      if (! is_part_of_previous_event){
        bernfebdaq::BernZMQEvent zmqEvent;
        mf::LogInfo("CRTDaqSim")<<"Creating New Fragment";
        // TODO: Put this into the configuration
        zmqEvent.mac5 = 5;
        zmqEvent.flags = 3;
        zmqEvent.lostcpu=0;
        zmqEvent.lostfpga=0;
        zmqEvent.ts0 = crtDat.T0();
        zmqEvent.ts1 = crtDat.T1();
        for(int tmp_adc =0; tmp_adc<32; tmp_adc++)
        {
          zmqEvent.adc[tmp_adc]=0;
        }
        zmqEvent.adc[adc_n]=crtDat.ADC();
        fEvents[feb_n].push_back(zmqEvent);
        fMetadata[feb_n].increment(0,0,0,1,0);
      }
    }

    // Now, Poll, if it's time to poll
    if(evt.time().value()>=fMetadata[0].time_end_seconds()){

      uint32_t ts_s = fMetadata[0].time_end_seconds();
      uint32_t ts_ns = 0;

      uint32_t te_s = ts_s + fPollingTime;
      uint32_t te_ns = 0;
      mf::LogInfo("CRTDaqSim")<<"Polling Data at time: "<<ts_s;

      for(int feb_n=0; feb_n<N_CRT_FEBS; feb_n++){

        auto frag_buf = artdaq::Fragment::FragmentBytes(fMetadata[feb_n].n_events()*sizeof(bernfebdaq::BernZMQEvent),  
                  fMetadata[feb_n].sequence_number(),feb_n,
                  bernfebdaq::detail::FragmentType::BernZMQ, fMetadata[feb_n]);

        frags->emplace_back( *frag_buf );
        std::copy(fEvents[feb_n].begin(),fEvents[feb_n].end(),(bernfebdaq::BernZMQEvent*)(frags->back().dataBegin()));

        //TODO: This might have to be moved to after the evt.put function

        fMetadata[feb_n] = bernfebdaq::BernZMQFragmentMetadata(ts_s, ts_ns, 
          te_s, te_ns,
          fTimeCorrectionDiff, fTimeOffset,
          evt.id().run(), 0,
          feb_n, 0,
          32, 32);
        fEvents[feb_n] = std::vector<bernfebdaq::BernZMQEvent>();
      }
    }
    evt.put(std::move(frags), fInstanceName);
  }

  void CRTDaqSim::reconfigure(fhicl::ParameterSet const& pSet)
  {
    fProducerName = pSet.get<std::string>("InputTag","crt_detsim");
    fPollingTime = pSet.get<unsigned>("PollingTime",1);
    fTimeCorrectionDiff = pSet.get<unsigned>("TimeCorrectionDiff",0);
    fTimeOffset = pSet.get<unsigned>("TimeOffset",0);
    fCoincidenceTime = pSet.get<unsigned>("CoincidenceTime",0);
  }

  void CRTDaqSim::beginRun(art::Run& run)
  {

  }

  void CRTDaqSim::beginSubRun(art::SubRun& subrun)
  {
    //TODO: Verify that the time conversion here is correct
    uint32_t ts_s = subrun.beginTime().value();
    uint32_t ts_ns = 0;

    uint32_t te_s = ts_s + fPollingTime;
    uint32_t te_ns = 0;

    uint32_t rid=0; //TODO All readers are 0 at the moment

    uint32_t nch=32;


    for(int feb_n=0;feb_n<N_CRT_FEBS;feb_n++)
    {
      fMetadata[feb_n] = bernfebdaq::BernZMQFragmentMetadata(ts_s, ts_ns, 
        te_s, te_ns,
        fTimeCorrectionDiff, fTimeOffset,
        subrun.id().run(), 0,
        feb_n, rid,
        nch, 32);
      fEvents[feb_n] = std::vector<bernfebdaq::BernZMQEvent>();
    }

  }

  void CRTDaqSim::endRun(art::Run& run)
  {
  }

  void CRTDaqSim::endSubRun(art::SubRun& subrun)
  {
  }


}

DEFINE_ART_MODULE(crt::CRTDaqSim)
