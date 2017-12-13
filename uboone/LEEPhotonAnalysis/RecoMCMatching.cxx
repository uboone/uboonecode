
#include "RecoMCMatching.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"




RecoMCMatch::RecoMCMatch(int imc_type, size_t imc_index, double iratio) :
  mc_type(imc_type),
  mc_index(imc_index),
  ratio(iratio){}



RecoMCMatching::RecoMCMatching() :
  fconsider_mcparticles(false),
  fhit_tree(nullptr),
  fsimch_tree(nullptr) {

}



void RecoMCMatching::FillHitTree(art::Event const & e) {

  if(!fhit_tree) {
    art::ServiceHandle<art::TFileService> tfs;
    fhit_tree = tfs->make<TTree>("hit_tree", "");
    fhit_tree->Branch("hit_from_reco_track", &fhit_from_reco_track, "hit_from_reco_track/I");
    fhit_tree->Branch("hit_time", &fhit_time, "hit_time/D");
    fhit_tree->Branch("hit_tdc", &fhit_tdc, "hit_tdc/D");
  }

  art::ValidHandle<std::vector<recob::Track>> const & ev_track =
    e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);
  art::FindManyP<recob::Hit> TrackToHit(ev_track, e, ftrack_producer);  

  art::ValidHandle<std::vector<recob::Shower>> const & ev_shower =
    e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);
  art::FindManyP<recob::Hit> ShowerToHit(ev_shower, e, fshower_producer);  
  
  detinfo::DetectorClocks const * ts = lar::providerFrom<detinfo::DetectorClocksService>();

  for(size_t i = 0; i < ev_track->size(); ++i) {
    for(art::Ptr<recob::Hit> const & hit : TrackToHit.at(i)) {
      fhit_from_reco_track = 1;
      fhit_time = hit->PeakTime();
      fhit_tdc = ts->TPCTick2TDC(fhit_time);
      fhit_tree->Fill();
    }
  }

  for(size_t i = 0; i < ev_shower->size(); ++i) {
    for(art::Ptr<recob::Hit> const & hit : ShowerToHit.at(i)) {
      fhit_from_reco_track = 0;
      fhit_time = hit->PeakTime();
      fhit_tdc = ts->TPCTick2TDC(fhit_time);
      fhit_tree->Fill();
    }
  }

}



void RecoMCMatching::FillSimchTree(art::Event const & e) {

  if(!fsimch_tree) {
    art::ServiceHandle<art::TFileService> tfs;
    fsimch_tree = tfs->make<TTree>("simch_tree", "");
    fsimch_tree->Branch("simch_tdc", &fsimch_tdc, "simch_tdc/I");
  }

  art::ValidHandle<std::vector<sim::SimChannel>> const & ev_simch =
    e.getValidHandle<std::vector<sim::SimChannel>>(fsimch_producer);
  
  for(sim::SimChannel const & simch : *ev_simch) {
    for(std::pair<unsigned short, std::vector<sim::IDE> > const & time_ide : simch.TDCIDEMap()) {
      fsimch_tdc = time_ide.first;
      fsimch_tree->Fill();
    }	  
  }
}



void RecoMCMatching::FillRatio(std::vector<RecoMCMatch> & ratio_v,
			       size_t const ev_mctrack_size,
			       size_t const ev_mcshower_size,
			       size_t const ratio_index,
			       size_t const best_mc_index,
			       double const ratio,
			       bool const last) {

  RecoMCMatch rmcm(0, best_mc_index, ratio);

  if(last) {
    rmcm.mc_index = SIZE_MAX;
  }
  else if(best_mc_index < ev_mctrack_size) {
    rmcm.mc_type = 1;
  }
  else if(best_mc_index >= ev_mctrack_size && best_mc_index < ev_mctrack_size + ev_mcshower_size) {
    rmcm.mc_type = 2;
    rmcm.mc_index = best_mc_index - ev_mctrack_size;
  }
  else if(best_mc_index != SIZE_MAX) {
    rmcm.mc_type = 3;
  }

  ratio_v.push_back(rmcm);

}



void RecoMCMatching::MatchSimchInfo(art::Event const & e,
				    std::vector< std::vector< art::Ptr< recob::Hit > > > const & reco_to_hit_v) {

  fmcq_vv.clear();
  ftrack_matches.clear();
  fshower_matches.clear();

  if(reco_to_hit_v.size() == 0) return;

  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack =
    e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower =
    e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle =
    e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  if(ev_mctrack->empty() && ev_mcshower->empty() && !(fconsider_mcparticles && !ev_mcparticle->empty())) return;

  art::ValidHandle<std::vector<recob::Track>> const & ev_track =
    e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);  
  art::ValidHandle<std::vector<sim::SimChannel>> const & ev_simchannel =
    e.getValidHandle<std::vector<sim::SimChannel>>(fsimch_producer);
  
  std::vector<std::vector<unsigned int>> g4_trackid_v;
  std::vector<unsigned int> entered_trackid_v;

  for(size_t mc_index = 0; mc_index < ev_mctrack->size(); ++mc_index) {
    sim::MCTrack const & mctrack = ev_mctrack->at(mc_index);
    std::vector<unsigned int> id_v;
    id_v.push_back(mctrack.TrackID());
    g4_trackid_v.push_back(id_v);
    entered_trackid_v.push_back(mctrack.TrackID());
  }
  
  std::vector<unsigned int> entered_trackid_v2;
  if(fconsider_mcparticles) entered_trackid_v2 = entered_trackid_v;
  for(size_t mc_index = 0; mc_index < ev_mcshower->size(); ++mc_index){
    sim::MCShower const & mcshower = ev_mcshower->at(mc_index);
    std::vector<unsigned int> id_v;
    id_v.reserve(mcshower.DaughterTrackID().size());
    for(auto const& id : mcshower.DaughterTrackID()) {
      if(id == mcshower.TrackID()) continue;
      if(std::find(entered_trackid_v.begin(), entered_trackid_v.end(), id) != entered_trackid_v.end())
	continue;
      id_v.push_back(id);
      if(fconsider_mcparticles) entered_trackid_v2.push_back(id);
    }
    id_v.push_back(mcshower.TrackID());
    if(fconsider_mcparticles) entered_trackid_v2.push_back(mcshower.TrackID());
    g4_trackid_v.push_back(id_v);
  }  

  std::vector<size_t> mcparticle_index_v;
  if(fconsider_mcparticles) {
    for(size_t mc_index = 0; mc_index < ev_mcparticle->size(); ++mc_index) {
      simb::MCParticle const & mcp = ev_mcparticle->at(mc_index);
      auto const id = mcp.TrackId();
      if(std::find(entered_trackid_v2.begin(), entered_trackid_v2.end(), id) != entered_trackid_v2.end())
	continue; 
      std::vector<unsigned int> id_v;
      id_v.push_back(id);
      mcparticle_index_v.push_back(mc_index);
      g4_trackid_v.push_back(id_v);
    }
  }

  simchi.FillSimchInfo(*ev_simchannel,
		       g4_trackid_v,
		       reco_to_hit_v);
  fmcq_vv = simchi.GetMCQVV();
  
  for(size_t reco_index = 0; reco_index < reco_to_hit_v.size(); ++reco_index) {
       
    std::vector<double> const  mcq_v = fmcq_vv.at(reco_index);
    size_t const mcsize = mcq_v.size();

    double total_mcq = 0;
    double max_mcq = 0;
    bool last = false;    
    size_t best_mc_index = SIZE_MAX;

    for(size_t mcq_index = 0; mcq_index < mcsize; ++mcq_index) {
      double const mcq = mcq_v.at(mcq_index);
      total_mcq += mcq;
      if(mcq > max_mcq) {
	max_mcq = mcq;
	best_mc_index = mcq_index;
	if(mcq_index == mcsize - 1) last = true;
      }
    }  

    if(fconsider_mcparticles && !last && best_mc_index != SIZE_MAX && best_mc_index >= ev_mctrack->size() + ev_mcshower->size()) {
      best_mc_index = mcparticle_index_v.at(best_mc_index - ev_mctrack->size() - ev_mcshower->size());
    }

    double ratio;
    if(total_mcq == 0) ratio = 0;
    else ratio = max_mcq / total_mcq;
    if(reco_index < ev_track->size()) FillRatio(ftrack_matches, ev_mctrack->size(), ev_mcshower->size(), reco_index, best_mc_index, ratio, last);
    else FillRatio(fshower_matches, ev_mctrack->size(), ev_mcshower->size(), reco_index - ev_track->size(), best_mc_index, ratio, last);

  }

}



void RecoMCMatching::MatchAll(art::Event const & e) {  

  art::ValidHandle<std::vector<recob::Track>> const & ev_track =
    e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);
  art::FindManyP<recob::Hit> TrackToHit(ev_track, e, ftrack_producer);
  art::ValidHandle<std::vector<recob::Shower>> const & ev_shower =
    e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);
  art::FindManyP<recob::Hit> ShowerToHit(ev_shower, e, fshower_producer);

  std::vector<std::vector<art::Ptr<recob::Hit>>> reco_to_hit_v;
  reco_to_hit_v.reserve(TrackToHit.size() + ShowerToHit.size());

  for(size_t i = 0; i < TrackToHit.size(); ++i) {
    reco_to_hit_v.push_back(TrackToHit.at(i));
  }

  for(size_t i = 0; i < ShowerToHit.size(); ++i) {
    reco_to_hit_v.push_back(ShowerToHit.at(i));
  }

  if(fmatch_type == "" || fmatch_type == "SimchInfo") MatchSimchInfo(e, reco_to_hit_v);

}



void RecoMCMatching::CoutMatches(art::Event const & e) {

  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle =
    e.getValidHandle<std::vector<simb::MCParticle>>("largeant");  
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack =
    e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower =
    e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  art::ValidHandle<std::vector<recob::Track>> const & ev_track =
    e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);
  art::ValidHandle<std::vector<recob::Shower>> const & ev_shower =
    e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);

  std::cout << "==================================================================================\n\n"
	    << "Event: " << e.id().event() << "\n"
	    << "Reco Tracks: " << ev_track->size() << " Reco Showers: " << ev_shower->size() << "\n"
	    << "MCTracks: " << ev_mctrack->size() << " MCShowers: " << ev_mcshower->size();
  if(fconsider_mcparticles) std::cout << " MCParticles: " << ev_mcparticle->size();
  std::cout << "\n\n";

  for(size_t i = 0; i < fmcq_vv.size(); ++i) {
    if(i < ev_track->size()) std::cout << "Track: " << i;
    else std::cout << "Shower: " << i - ev_track->size();
    std::cout << " - MCTracks: ";
    bool temp_mcshower = true;
    bool temp_mcparticle = true;
    double total_mcq = 0;
    for(size_t j = 0; j < fmcq_vv.at(i).size(); ++j) {
      if(j == fmcq_vv.at(i).size()-1) {
	std::cout << "Other: "; 
      }
      std::cout << fmcq_vv.at(i).at(j) << " ";
      total_mcq += fmcq_vv.at(i).at(j);
      if(j >= ev_mctrack->size()-1 && temp_mcshower) {
	std::cout << "MCShowers: ";
	temp_mcshower = false;
      }
      if(fconsider_mcparticles && j >= ev_mctrack->size() + ev_mcshower->size() - 1 && temp_mcparticle) {
	std::cout << "MCParticles: ";
	temp_mcparticle = false;
      }
    }
    std::cout << "Total: " << total_mcq << "\n";
  }

  std::cout << "\nTrack Ratios:\n";
  for(size_t i = 0; i < ftrack_matches.size(); ++i) {
    RecoMCMatch const & rmcm  = ftrack_matches.at(i);
    std::cout << "Track: " << i << " - ";
    if(rmcm.mc_type == 0) std::cout << "No match: ";
    else if(rmcm.mc_type == 1) std::cout << "MCTrack index: ";
    else if(rmcm.mc_type == 2) std::cout << "MCShower index: ";
    else if(rmcm.mc_type == 3) std::cout << "MCParticle index: ";
    std::cout << rmcm.mc_index << " ratio: " << rmcm.ratio << "\n";
  }

  std::cout << "\nShower Ratios:\n";
  for(size_t i = 0; i < fshower_matches.size(); ++i) {
    RecoMCMatch const & rmcm  = fshower_matches.at(i);
    std::cout << "Shower: " << i << " - ";
    if(rmcm.mc_type == 0) std::cout << "No match: ";
    else if(rmcm.mc_type == 1) std::cout << "MCTrack index: ";
    else if(rmcm.mc_type == 2) std::cout << "MCShower index: ";
    else if(rmcm.mc_type == 3) std::cout << "MCParticle index: ";
    std::cout << rmcm.mc_index << " ratio: " << rmcm.ratio << "\n";
  }

  std::cout << "\n";

}
