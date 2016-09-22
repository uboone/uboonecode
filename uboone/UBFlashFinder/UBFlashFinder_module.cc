////////////////////////////////////////////////////////////////////////
// Class:       UBFlashFinder
// Module Type: producer
// File:        UBFlashFinder_module.cc
//
// Generated at Tue Sep 13 22:30:26 2016 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/RecoBase/OpHit.h"
#include "lardata/RecoBase/OpFlash.h"

#include <memory>
#include <string>
#include "FlashFinderManager.h"
#include "FlashFinderFMWKInterface.h"
#include "PECalib.h"

class UBFlashFinder;

class UBFlashFinder : public art::EDProducer {
public:
  explicit UBFlashFinder(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBFlashFinder(UBFlashFinder const &) = delete;
  UBFlashFinder(UBFlashFinder &&) = delete;
  UBFlashFinder & operator = (UBFlashFinder const &) = delete;
  UBFlashFinder & operator = (UBFlashFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  // Declare member data here.
  ::pmtana::FlashFinderManager _mgr;
  ::pmtana::PECalib _pecalib;
  std::string _flash_producer;
  std::string _hit_producer;

  bool _beam_flash;
};


UBFlashFinder::UBFlashFinder(pmtana::Config_t const & p)
// :
// Initialize member data here.
{
  _hit_producer   = p.get<std::string>("OpHitProducer");
  _flash_producer = p.get<std::string>("OpFlashProducer");
  
  auto const flash_algo  = p.get<std::string>("FlashFinderAlgo");
  auto const flash_pset = p.get<pmtana::Config_t>("AlgoConfig");
  auto algo_ptr = ::pmtana::FlashAlgoFactory::get().create(flash_algo,flash_algo);
  algo_ptr->Configure(flash_pset);
  _mgr.SetFlashAlgo(algo_ptr);

  _pecalib.Configure(p.get<pmtana::Config_t>("PECalib"));

  _beam_flash = p.get<bool>("BeamFlash");

  produces< std::vector<recob::OpFlash>   >();

}

void UBFlashFinder::produce(art::Event & e)
{

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<recob::OpFlash> > opflashes(new std::vector<recob::OpFlash>);

  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit> > ophit_h;
  e.getByLabel(_hit_producer,ophit_h);

  // make sure hits look good
  if(!ophit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate OpHit!"<<std::endl;
    throw std::exception();
  }

  ::pmtana::LiteOpHitArray_t ophits;
  double trigger_time=1.1e20;
  for(auto const& oph : *ophit_h) {
    ::pmtana::LiteOpHit_t loph;
    if(trigger_time > 1.e20) trigger_time = oph.PeakTimeAbs() - oph.PeakTime();
    loph.peak_time = oph.PeakTime();

    size_t opdet = ::pmtana::OpDetFromOpChannel(oph.OpChannel());
    loph.pe = ( _beam_flash ? _pecalib.BeamPE(opdet,oph.Area(),oph.Amplitude()) : _pecalib.CosmicPE(opdet,oph.Area(),oph.Amplitude()) );

    loph.channel = oph.OpChannel();
    ophits.emplace_back(std::move(loph));
  }
  
  auto const flash_v = _mgr.RecoFlash(ophits);

  for(const auto& lflash :  flash_v) {

    recob::OpFlash flash(lflash.time, lflash.time_err,
			 trigger_time + lflash.time,
			 (trigger_time + lflash.time) / 1600.,
			 lflash.channel_pe);
    opflashes->emplace_back(std::move(flash));
  }

  e.put(std::move(opflashes));
}

DEFINE_ART_MODULE(UBFlashFinder)
