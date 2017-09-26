#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/CRTSimData.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"

namespace crt{

  CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): 
  fFileNames(pset.get< std::vector<std::string> >("InputFilenames")),
  fCRTEvent(fFileNames)
  {
    this->reconfigure(pset);
    produces< std::vector<CRTSimData> >();
  }

  CRTMerger::~CRTMerger()
  {

  }


  void beginRun(art::Run &)
  {
    
  }
 
  void endRun(art::Run &)
  {
    
  }

  void CRTMerger::produce(art::Event& event)
  {
    std::unique_ptr<std::vector<crt::CRTSimData> > crtHits(
        new std::vector<crt::CRTSimData>);

    //For this event
    auto event_t =  event.time().value();
    mf::LogInfo("CRTMerger")<<"Finding CRT Events for time : "<<event_t<<std::endl;;

    if(fCRTEvent.atEnd()) return;

    auto const& crthits = *(fCRTEvent.getValidHandle< std::vector<artdaq::Fragment> >(fTag));

    while(!fCRTEvent.atEnd()){
      for(auto crtFrag : crthits){
        bernfebdaq::BernZMQFragment bernfrag(crtFrag);
        auto feb_id = bernfrag.metadata()->feb_id();
        auto pc_id = bernfrag.metadata()->reader_id();

        for(unsigned i = 0; i<bernfrag.metadata()->n_events() ; i++){
          auto crt_time = bernfrag.eventdata(i)->Time_TS0();

          if (crt_time>= event_t-fTimeWindow/2 && crt_time< event_t+fTimeWindow/2 ){

            for(auto adc_channel =0; adc_channel<32; adc_channel++){

              unsigned channel = pc_id*1000+feb_id*100+adc_channel;
              crt::CRTSimData myNewHit(channel, bernfrag.eventdata(i)->Time_TS0(), bernfrag.eventdata(i)->Time_TS1(), bernfrag.eventdata(i)->ADC(adc_channel));
              crtHits->push_back(myNewHit);
            }
          }
          else if (crt_time>=event_t+fTimeWindow/2){
            event.put(std::move(crtHits));
            return;
          }
        }
      }
    fCRTEvent.next();
    }
    event.put(std::move(crtHits));
  }

  void CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCRTEvent.toBegin();
    fTag = {pset.get<std::string>("InputTagName","daq")};
    fTimeWindow = pset.get<unsigned>("TimeWindow",10);
  }


/*
  void CRTMerger::beginRun(art::Run const&)
  {

  }
  */
}

DEFINE_ART_MODULE(crt::CRTMerger)
