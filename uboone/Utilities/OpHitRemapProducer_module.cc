////////////////////////////////////////////////////////////////////////
// Class:       OpHitRemapProducer
// Plugin Type: producer (art v2_05_01)
// File:        OpHitRemapProducer_module.cc
//
// Generated at Wed Mar 28 17:17:14 2018 by Giuseppe Cerati using cetskelgen
// from cetlib version v1_21_00.
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

#include "lardataobj/RecoBase/OpHit.h"

#include <memory>

class OpHitRemapProducer;


class OpHitRemapProducer : public art::EDProducer {
public:
  explicit OpHitRemapProducer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpHitRemapProducer(OpHitRemapProducer const &) = delete;
  OpHitRemapProducer(OpHitRemapProducer &&) = delete;
  OpHitRemapProducer & operator = (OpHitRemapProducer const &) = delete;
  OpHitRemapProducer & operator = (OpHitRemapProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  art::InputTag inputTag_;
  std::vector<int> opchanmap_;
  unsigned int offset_;
};


OpHitRemapProducer::OpHitRemapProducer(fhicl::ParameterSet const & p)
{
  inputTag_ = p.get<art::InputTag>("OpHitsInputTag");
  opchanmap_ = p.get<std::vector<int> >("OpChannelSwapMap");
  offset_ = p.get<unsigned int>("MapOffset");
  if (offset_>0) {
    for (size_t i = 0; i<opchanmap_.size(); ++i) opchanmap_[i] = opchanmap_[i]+offset_;
    opchanmap_.insert(opchanmap_.begin(),offset_,0);
  }
  //
  produces<std::vector<recob::OpHit> >();
}

void OpHitRemapProducer::produce(art::Event & e)
{
  auto output = std::make_unique<std::vector<recob::OpHit> >();
  const auto& input = e.getValidHandle<std::vector<recob::OpHit> >(inputTag_);
  for (const auto& oph : *input){
    //std::cout << "original ophit OpChannel=" << oph.OpChannel() << " PeakTimeAbs=" << oph.PeakTimeAbs() << " PeakTime=" << oph.PeakTime() << " Frame=" << oph.Frame()
    //      << " Width=" << oph.Width() << " Area=" << oph.Area() << " Amplitude=" << oph.Amplitude() << " PE=" << oph.PE() << " FastToTotal=" << oph.FastToTotal() 
    //	      << std::endl;
    output->emplace_back( recob::OpHit(opchanmap_[oph.OpChannel()], 
				       oph.PeakTime(), 
				       oph.PeakTimeAbs(), 
				       oph.Frame(), 
				       oph.Width(), 
				       oph.Area(), 
				       oph.Amplitude(), 
				       oph.PE(), 
				       oph.FastToTotal() ) );
    //std::cout << "remapped ophit OpChannel=" << output->back().OpChannel() << " PeakTimeAbs=" << output->back().PeakTimeAbs() << " PeakTime=" << output->back().PeakTime() << " Frame=" << output->back().Frame()
    //	      << " Width=" << output->back().Width() << " Area=" << output->back().Area() << " Amplitude=" << output->back().Amplitude() << " PE=" << output->back().PE() << " FastToTotal=" << output->back().FastToTotal() 
    //	      << std::endl;
  }
  e.put(std::move(output));
}

DEFINE_ART_MODULE(OpHitRemapProducer)
