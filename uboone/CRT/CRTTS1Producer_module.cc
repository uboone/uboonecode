////////////////////////////////////////////////////////////////////////
// Class:       CRTTS1Producer
// Plugin Type: producer (art v2_06_03)
// File:        CRTTS1Producer_module.cc
//
// Generated at Tue Jul  4 11:13:52 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

namespace crt {
  class CRTTS1Producer;
}


class crt::CRTTS1Producer : public art::EDProducer {
public:
  explicit CRTTS1Producer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTS1Producer(CRTTS1Producer const &) = delete;
  CRTTS1Producer(CRTTS1Producer &&) = delete;
  CRTTS1Producer & operator = (CRTTS1Producer const &) = delete;
  CRTTS1Producer & operator = (CRTTS1Producer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  art::InputTag  fCRTDataLabel;


  TTree*       fTree;
  uint32_t     ev_time_s;
  uint32_t     ev_time_ns;
  uint32_t     ts1_time_s;
  uint32_t     ts1_time_ns;
  long int     ev_time;
  long int     ts1_time;
  long int     diff_time;
  uint32_t     n_crt_frags;
  

};


crt::CRTTS1Producer::CRTTS1Producer(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  produces< std::vector<bernfebdaq::BernZMQEvent> >();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tr_ts1","TS1Producer Tree");
  fTree->Branch("evt_s",&ev_time_s,"evt_s/i");
  fTree->Branch("evt_ns",&ev_time_ns,"evt_ns/i");
  fTree->Branch("evt_time",&ev_time,"evt_time/L");
  fTree->Branch("ts1_s",&ts1_time_s,"ts1_s/i");
  fTree->Branch("ts1_ns",&ts1_time_ns,"ts1_ns/i");
  fTree->Branch("ts1_time",&ts1_time,"ts1_time/L");
  fTree->Branch("diff_time",&diff_time,"diff_time/L");
  fTree->Branch("n_crt_frags",&n_crt_frags,"n_crt_frags/i");


}

void crt::CRTTS1Producer::produce(art::Event & e)
{

  ev_time_s = e.time().timeHigh();
  ev_time_ns = e.time().timeLow();
  ev_time = ((long int)ev_time_s)*1000000000 + (long int)ev_time_ns;

  ts1_time_s = 0;
  ts1_time_ns = 0;
  ts1_time = 0;
  n_crt_frags=0;
  diff_time = -999999999999;
  
  std::unique_ptr<std::vector<bernfebdaq::BernZMQEvent> > ts1Col(new std::vector<bernfebdaq::BernZMQEvent>);

  art::Handle< std::vector<artdaq::Fragment> > fragHandle;
  e.getByLabel(fCRTDataLabel,fragHandle);
  std::vector<artdaq::Fragment> const& fragVector(*fragHandle);

  bool first=true;
  n_crt_frags = fragVector.size();
  for(auto const& frag : fragVector){
    bernfebdaq::BernZMQFragment bfrag(frag);

    for(size_t i_e=0; i_e<bfrag.metadata()->n_events(); ++i_e)
      if(bfrag.eventdata(i_e)->IsReference_TS1()){
	ts1Col->push_back(bfrag.Event(i_e));
	if(first){
	  ts1_time_s  = bfrag.metadata()->time_start_seconds();
	  ts1_time_ns = ts1Col->back().Time_TS0();
	  ts1_time = ((long int)ts1_time_s)*1000000000 + (long int)ts1_time_ns;
	  diff_time = ts1_time - ev_time;
	  first=false;
	}
      }
  }

  fTree->Fill();

  e.put(std::move(ts1Col));
}

void crt::CRTTS1Producer::reconfigure(fhicl::ParameterSet const & p)
{
  fCRTDataLabel = p.get<art::InputTag>("CRTDataLabel");
}

DEFINE_ART_MODULE(crt::CRTTS1Producer)
