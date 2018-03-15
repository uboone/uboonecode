

#ifndef FILTERSIGNAL_CXX
#define FILTERSIGNAL_CXX


#include "FilterSignal.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "nusimdata/SimulationBase/MCTruth.h"



bool FilterSignal::RunOld(art::Event const & e, size_t & delta_mct_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  size_t signal_counter = 0;

  for(size_t mct_index = 0; mct_index < ev_mct->size(); ++mct_index) {

    simb::MCTruth const & mct = ev_mct->at(mct_index);

    if(mct.GetNeutrino().CCNC() == 0) continue;

    size_t external_photon_parent_index = SIZE_MAX;
    bool continue_bool = false;

    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      if(mcp.TrackId() != i) {
	std::cout << "ERROR: " << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nTrackId does not match index\n";
	exit(1);
      }
      if(mcp.StatusCode() != 1) continue;
      switch(abs(mcp.PdgCode())) {
      case 22: {
	if(external_photon_parent_index == SIZE_MAX) {
	  external_photon_parent_index = mcp.Mother();
	}
	else {
	  continue_bool = true;
	}
	break;
      }
      case 12:
      case 14:
      case 2112:
      case 2212:
	break;
      default:
	continue_bool = true;
      }
      if(continue_bool) break;
    }

    if(external_photon_parent_index == SIZE_MAX || continue_bool) 
      continue;   

    if(mct.GetParticle(external_photon_parent_index).PdgCode() == 22)
      external_photon_parent_index = mct.GetParticle(external_photon_parent_index).Mother();

    if(!(abs(mct.GetParticle(external_photon_parent_index).PdgCode()) == 2214 || abs(mct.GetParticle(external_photon_parent_index).PdgCode()) == 2114))
      continue;

    delta_mct_index = mct_index;
    ++signal_counter;

  }

  return signal_counter > 0;

}



bool FilterSignal::RunSingle(art::Event const & e, size_t const & mct_index, size_t & exiting_photon_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  simb::MCTruth const & mct = ev_mct->at(mct_index);
  if(mct.GetNeutrino().CCNC() == 0) return false;
  size_t external_photon_parent_index = SIZE_MAX;

  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(mcp.TrackId() != i) {
      std::cout << "ERROR: " << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nTrackId does not match index\n";
      exit(1);
    }
    if(mcp.StatusCode() != 1) continue;
    switch(abs(mcp.PdgCode())) {
    case 22: {
      if(external_photon_parent_index == SIZE_MAX) {
	exiting_photon_index = i;
	external_photon_parent_index = mcp.Mother();
      }
      else return false;
      break;
    }
    case 12:
    case 14:
    case 2112:
    case 2212:
      break;
    default:
      return false;
    }
  }

  if(external_photon_parent_index == SIZE_MAX) return false;

  if(mct.GetParticle(external_photon_parent_index).PdgCode() == 22)
    external_photon_parent_index = mct.GetParticle(external_photon_parent_index).Mother();

  if(!(abs(mct.GetParticle(external_photon_parent_index).PdgCode()) == 2214 || abs(mct.GetParticle(external_photon_parent_index).PdgCode()) == 2114))
    return false;

  return true;

}



bool FilterSignal::Run(art::Event const & e, size_t & delta_mct_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  size_t signal_counter = 0;
  size_t temp = 0;

  for(size_t mct_index = 0; mct_index < ev_mct->size(); ++mct_index) {
    if(!RunSingle(e, mct_index, temp)) continue;
    delta_mct_index = mct_index;
    ++signal_counter;
  }

  return signal_counter > 0;

}

#endif
