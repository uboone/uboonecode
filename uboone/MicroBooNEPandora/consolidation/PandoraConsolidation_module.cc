////////////////////////////////////////////////////////////////////////
// Class:       PandoraConsolidation
// Plugin Type: producer (art v2_07_03)
// File:        PandoraConsolidation_module.cc
//
// Generated at Thu Sep 21 13:56:07 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
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

class PandoraConsolidation;


class PandoraConsolidation : public art::EDProducer {
public:
  explicit PandoraConsolidation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraConsolidation(PandoraConsolidation const &) = delete;
  PandoraConsolidation(PandoraConsolidation &&) = delete;
  PandoraConsolidation & operator = (PandoraConsolidation const &) = delete;
  PandoraConsolidation & operator = (PandoraConsolidation &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};

PandoraConsolidation::PandoraConsolidation(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

void PandoraConsolidation::produce(art::Event & e)
{
  // Implementation of required member function here.
}

void PandoraConsolidation::beginJob()
{
  // Implementation of optional member function here.
}

void PandoraConsolidation::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PandoraConsolidation)
