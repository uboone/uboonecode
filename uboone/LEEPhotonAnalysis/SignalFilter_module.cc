////////////////////////////////////////////////////////////////////////
// Class:       SignalFilter
// Plugin Type: filter (art v2_05_00)
// File:        SignalFilter_module.cc
//
// Generated at Wed Oct 18 03:11:59 2017 by Robert Murrells using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "nusimdata/SimulationBase/MCTruth.h"

#include "FilterSignal.h"

class SignalFilter;


class SignalFilter : public art::EDFilter {

  FilterSignal ffs;

  unsigned int counter;
  unsigned int counterold;

public:
  explicit SignalFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SignalFilter(SignalFilter const &) = delete;
  SignalFilter(SignalFilter &&) = delete;
  SignalFilter & operator = (SignalFilter const &) = delete;
  SignalFilter & operator = (SignalFilter &&) = delete;

  void endJob();  

  // Required functions.
  bool filter(art::Event & e) override;
  void cout_particles(art::Event & e);

};


SignalFilter::SignalFilter(fhicl::ParameterSet const & p) :
  counter(0),
  counterold(0){}


void SignalFilter::cout_particles(art::Event & e) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  for(simb::MCTruth const mct : *ev_mct) {

    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout << mcp.TrackId() << " " << mcp.PdgCode() << " " << mcp.Mother() << " " << mcp.StatusCode() << "\n";
    }  
    
  }

}


bool SignalFilter::filter(art::Event & e) {

  size_t temp;
  
  //if(ffs.Run(e, temp)) cout_particles(e);

  if(ffs.Run(e, temp)) ++counter;
  if(ffs.RunOld(e, temp)) ++counterold;

  return ffs.Run(e, temp);
  
}



void SignalFilter::endJob() {

  std::cout << counter << " " << counterold << "\n";

}



DEFINE_ART_MODULE(SignalFilter)
