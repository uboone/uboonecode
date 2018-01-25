

#ifndef VERTEXQUALITY_CXX
#define VERTEXQUALITY_CXX

#include "VertexQuality.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"



VertexQuality::VertexQuality(std::string const & track_producer,
			     std::string const & shower_producer,
			     RecoMCMatching const & rmcm) :
  ftrack_producer(track_producer),
  fshower_producer(shower_producer),
  frmcm(&rmcm),
  ftpc_volume(0,
	      -lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
	      0,
	      2*lar::providerFrom<geo::Geometry>()->DetHalfWidth(),
	      lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
	      lar::providerFrom<geo::Geometry>()->DetLength()) {

  art::ServiceHandle< art::TFileService > tfs;
  fvertex_tree = tfs->make<TTree>("vertex_quality_tree", "");  

  fvertex_tree->Branch("tpc_volume_contained", &ftpc_volume_contained, "tpc_volume_contained/I");

  fvertex_tree->Branch("dist", &fdist, "dist/D");
  fvertex_tree->Branch("true_track_total", &ftrue_track_total, "true_track_total/I");
  fvertex_tree->Branch("true_shower_total", &ftrue_shower_total, "true_shower_total/I");
  fvertex_tree->Branch("reco_track_total", &freco_track_total, "reco_track_total/I");
  fvertex_tree->Branch("reco_shower_total", &freco_shower_total, "reco_shower_total/I");
  fvertex_tree->Branch("correct_track_total", &fcorrect_track_total, "correct_track_total/I");
  fvertex_tree->Branch("correct_shower_total", &fcorrect_shower_total, "correct_shower_total/I");

  fvertex_tree->Branch("track_matching_ratio_v", &ftrack_matching_ratio_v);
  fvertex_tree->Branch("track_true_pdg_v", &ftrack_true_pdg_v);
  fvertex_tree->Branch("track_true_origin_v", &ftrack_true_origin_v);

  fvertex_tree->Branch("shower_matching_ratio_v", &fshower_matching_ratio_v);
  fvertex_tree->Branch("shower_true_pdg_v", &fshower_true_pdg_v);
  fvertex_tree->Branch("shower_true_origin_v", &fshower_true_origin_v);

}



art::Ptr<simb::MCTruth> VertexQuality::TrackIDToMCTruth(art::Event const & e, int const geant_track_id) {

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, "largeant", truthToParticles, particlesToTruth);

  for(auto iter : particlesToTruth) {
    if(iter.first->TrackId() == geant_track_id) {
      return iter.second;
    }
  }

  std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nERROR: no mctruth found\n";
  exit(1);

}



void VertexQuality::GetTrueObjects(art::Event const & e,
				   std::vector<size_t> & mctrack_v,
				   std::vector<size_t> & mcshower_v,
				   std::vector<size_t> & mcparticle_v) {
  
  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");  
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack = e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  
  double const threshold = 1;

  geoalgo::Point_t const true_nu_vtx = ev_mctruth->front().GetNeutrino().Nu().Position(0);
  std::vector<int> used_trackid_v;

  for(size_t i = 0; i < ev_mctrack->size(); ++i) {
    sim::MCTrack const & mctr = ev_mctrack->at(i);
    if(true_nu_vtx.Dist(mctr.Start().Position()) > threshold) continue;
    mctrack_v.push_back(i);
    used_trackid_v.push_back(mctr.TrackID());
  }

  for(size_t i = 0; i < ev_mcshower->size(); ++i) {
    sim::MCShower const & mcs = ev_mcshower->at(i);
    if(std::find(used_trackid_v.begin(), used_trackid_v.end(), mcs.TrackID()) != used_trackid_v.end() || 
       true_nu_vtx.Dist(mcs.Start().Position()) > threshold) continue;
    mcshower_v.push_back(i);
    used_trackid_v.push_back(mcs.TrackID());
    for(unsigned int const trackid : mcs.DaughterTrackID()) {
      if(trackid == mcs.TrackID() ||
	 std::find(used_trackid_v.begin(), used_trackid_v.end(), mcs.TrackID()) != used_trackid_v.end()) continue;
      used_trackid_v.push_back(trackid);
    }
  }

  for(size_t i = 0; i < ev_mcparticle->size(); ++i) {
    simb::MCParticle const & mcp = ev_mcparticle->at(i);
    if(true_nu_vtx.Dist(mcp.Position(0)) > threshold || 
       std::find(used_trackid_v.begin(), used_trackid_v.end(), mcp.TrackId()) != used_trackid_v.end()) 
      continue;
    mcparticle_v.push_back(i);
    used_trackid_v.push_back(mcp.TrackId());
  }  

}



void VertexQuality::GetTrueObjects(art::Event const & e,
				   std::vector<int> & mcparticle_v) {
  
  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth = e.getValidHandle<std::vector<simb::MCTruth>>("generator");  
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");  

  double const threshold = 1;

  geoalgo::Point_t const true_nu_vtx = ev_mctruth->front().GetNeutrino().Nu().Position(0);

  for(size_t i = 0; i < ev_mcshower->size(); ++i) {
    sim::MCShower const & mcs = ev_mcshower->at(i);
    if(true_nu_vtx.Dist(mcs.Start().Position()) > threshold) continue;
    mcparticle_v.push_back(mcs.TrackID());
    for(unsigned int const trackid : mcs.DaughterTrackID()) {
      if(trackid == mcs.TrackID()) continue;
      mcparticle_v.push_back(trackid);
    }
  }

  for(size_t i = 0; i < ev_mcparticle->size(); ++i) {
    simb::MCParticle const & mcp = ev_mcparticle->at(i);
    if(true_nu_vtx.Dist(mcp.Position(0)) > threshold || 
       std::find(mcparticle_v.begin(), mcparticle_v.end(), mcp.TrackId()) != mcparticle_v.end()) 
      continue;
    mcparticle_v.push_back(mcp.TrackId());
  }  

}



double VertexQuality::GetTrueTotal(std::vector<int> const & mcparticle_v) {

  std::vector<RecoMCMatch> const & shower_matches = frmcm->GetShowerMatches();
  std::vector<RecoMCMatch> const & track_matches = frmcm->GetTrackMatches();

  double total = 0;

  for(int const trkid : mcparticle_v) {  
    for(RecoMCMatch const & rmcm : shower_matches) {
      std::unordered_map<int, double> const & trkide_map = rmcm.trkide_map;
      auto const tm_it = trkide_map.find(trkid);
      if(tm_it == trkide_map.end()) continue;
      total += tm_it->second;
    }
    for(RecoMCMatch const & rmcm : track_matches) {
      std::unordered_map<int, double> const & trkide_map = rmcm.trkide_map;
      auto const tm_it = trkide_map.find(trkid);
      if(tm_it == trkide_map.end()) continue;
      total += tm_it->second;
    }
  }

  return total;

}



void VertexQuality::GetTrueRecoObjects(art::Event const & e,
				       std::vector<size_t> & track_v,
				       std::vector<size_t> & shower_v) {

  std::vector<size_t> mctrack_v;
  std::vector<size_t> mcshower_v;
  std::vector<size_t> mcparticle_v;
  GetTrueObjects(e, mctrack_v, mcshower_v, mcparticle_v);

  std::vector<RecoMCMatch> const & shower_matches = frmcm->GetShowerMatches();
  std::vector<RecoMCMatch> const & track_matches = frmcm->GetTrackMatches(); 

  for(size_t i = 0; i < shower_matches.size(); ++i) {
    RecoMCMatch const & rmcm = shower_matches.at(i);
    if(rmcm.mc_type == frmcm->fmc_type_shower &&
       std::find(mcshower_v.begin(), mcshower_v.end(), rmcm.mc_index) != mcshower_v.end()) {
      shower_v.push_back(i);
    }
    else if(rmcm.mc_type == frmcm->fmc_type_track &&
	    std::find(mctrack_v.begin(), mctrack_v.end(), rmcm.mc_index) != mctrack_v.end()) {
      shower_v.push_back(i);
    }
    else if(rmcm.mc_type == frmcm->fmc_type_particle &&
	    std::find(mcparticle_v.begin(), mcparticle_v.end(), rmcm.mc_index) != mcparticle_v.end()) {
      shower_v.push_back(i);
    }
  }

  for(size_t i = 0; i < track_matches.size(); ++i) {
    RecoMCMatch const & rmcm = track_matches.at(i);
    if(rmcm.mc_type == frmcm->fmc_type_shower &&
       std::find(mcshower_v.begin(), mcshower_v.end(), rmcm.mc_index) != mcshower_v.end()) {
      track_v.push_back(i);
    }
    else if(rmcm.mc_type == frmcm->fmc_type_track &&
	    std::find(mctrack_v.begin(), mctrack_v.end(), rmcm.mc_index) != mctrack_v.end()) {
      track_v.push_back(i);
    }
    else if(rmcm.mc_type == frmcm->fmc_type_particle &&
	    std::find(mcparticle_v.begin(), mcparticle_v.end(), rmcm.mc_index) != mcparticle_v.end()) {
      track_v.push_back(i);
    }
  }

}



void VertexQuality::Reset() {

  freco_track_total = 0;
  fcorrect_track_total = 0;
  freco_shower_total = 0;
  fcorrect_shower_total = 0;
  
  ftrack_matching_ratio_v.clear();
  ftrack_true_pdg_v.clear();
  ftrack_true_origin_v.clear();

  fshower_matching_ratio_v.clear();
  fshower_true_pdg_v.clear();
  fshower_true_origin_v.clear();

}
 


void VertexQuality::RunDist(art::Event const & e,
			    ParticleAssociations const & pas) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack = e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  geoalgo::Point_t const true_nu_vtx = ev_mctruth->front().GetNeutrino().Nu().Position(0);

  ftpc_volume_contained = 0;
  if(ftpc_volume.Contain(true_nu_vtx)) ftpc_volume_contained = 1;

  std::vector<RecoMCMatch> const & shower_matches = frmcm->GetShowerMatches();
  std::vector<RecoMCMatch> const & track_matches = frmcm->GetTrackMatches(); 

  std::vector<size_t> track_v;
  std::vector<size_t> shower_v;
  GetTrueRecoObjects(e, track_v, shower_v);
  ftrue_track_total = track_v.size(); 
  ftrue_shower_total = shower_v.size(); 

  DetectorObjects const & deto_v = pas.GetDetectorObjects();
  std::vector<ParticleAssociation> const & pa_v = pas.GetAssociations();

  for(size_t const pa_index : pas.GetSelectedAssociations()) {

    Reset();

    ParticleAssociation const & pa = pa_v.at(pa_index);
    fdist = true_nu_vtx.Dist(pa.GetRecoVertex());

    for(size_t const object_index : pa.GetObjectIndices()) {
      DetectorObject const & deto = deto_v.GetDetectorObject(object_index);
      size_t const original_index = deto.foriginal_index;
      if(deto.freco_type == deto_v.ftrack_reco_type) {
	RecoMCMatch const & rmcm = track_matches.at(original_index);
	ftrack_matching_ratio_v.push_back(rmcm.ratio);
	if(rmcm.mc_type == frmcm->fmc_type_track) {
	  sim::MCTrack const & mctr = ev_mctrack->at(rmcm.mc_index);
	  ftrack_true_pdg_v.push_back(mctr.PdgCode());
	  ftrack_true_origin_v.push_back(mctr.Origin());
	}
	else if(rmcm.mc_type == frmcm->fmc_type_shower) {
	  sim::MCShower const & mcs = ev_mcshower->at(rmcm.mc_index);
	  ftrack_true_pdg_v.push_back(mcs.PdgCode());
	  ftrack_true_origin_v.push_back(mcs.Origin());
	}
	else if(rmcm.mc_type == frmcm->fmc_type_particle) {
	  simb::MCParticle const & mcp = ev_mcparticle->at(rmcm.mc_index);
	  ftrack_true_pdg_v.push_back(mcp.PdgCode());
	  ftrack_true_origin_v.push_back(TrackIDToMCTruth(e, mcp.TrackId())->Origin());
	}
	++freco_track_total;
	if(std::find(track_v.begin(), track_v.end(), original_index) != track_v.end()) {
	  ++fcorrect_track_total;
	}
      }
      if(deto.freco_type == deto_v.fshower_reco_type) {
	RecoMCMatch const & rmcm = shower_matches.at(original_index);
	fshower_matching_ratio_v.push_back(rmcm.ratio);
	if(rmcm.mc_type == frmcm->fmc_type_track) {
	  sim::MCTrack const & mctr = ev_mctrack->at(rmcm.mc_index);
	  fshower_true_pdg_v.push_back(mctr.PdgCode());
	  fshower_true_origin_v.push_back(mctr.Origin());
	}
	else if(rmcm.mc_type == frmcm->fmc_type_shower) {
	  sim::MCShower const & mcs = ev_mcshower->at(rmcm.mc_index);
	  fshower_true_pdg_v.push_back(mcs.PdgCode());
	  fshower_true_origin_v.push_back(mcs.Origin());
	}
	else if(rmcm.mc_type == frmcm->fmc_type_particle) {
	  simb::MCParticle const & mcp = ev_mcparticle->at(rmcm.mc_index);
	  fshower_true_pdg_v.push_back(mcp.PdgCode());
	  fshower_true_origin_v.push_back(TrackIDToMCTruth(e, mcp.TrackId())->Origin());
	}
	++freco_shower_total;
	if(std::find(shower_v.begin(), shower_v.end(), original_index) != shower_v.end()) {
	  ++fcorrect_shower_total;
	}
      }

    }

    fvertex_tree->Fill();

  }

}



void VertexQuality::Run(art::Event const & e,
			ParticleAssociations const & pas) {

  std::vector<int> mcparticle_v;
  GetTrueObjects(e, mcparticle_v);
  double const true_total = GetTrueTotal(mcparticle_v);

  std::vector<RecoMCMatch> const & shower_matches = frmcm->GetShowerMatches();
  std::vector<RecoMCMatch> const & track_matches = frmcm->GetTrackMatches();

  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcparticle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack = e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

  if(true_total == 0 && track_matches.size() != 0 && shower_matches.size() != 0) {

    art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth = e.getValidHandle<std::vector<simb::MCTruth>>("generator");

    std::cout << "true_total: " << true_total << " ev_mctruth->size() " << ev_mctruth->size() << " shower_matches.size(): " << shower_matches.size() << " track_matches.size(): " << track_matches.size() << "\ntrue_nu_vtx: " <<   geoalgo::Point_t(ev_mctruth->front().GetNeutrino().Nu().Position(0)) << "\n";

    std::cout << "\n";
    for(RecoMCMatch const & rmcm : shower_matches) {
      std::cout << "reco shower: " << rmcm.ratio << " origin: ";
      if(rmcm.mc_type == 1) {
	std::cout << ev_mcshower->at(rmcm.mc_index).Origin() << " start: " << geoalgo::Point_t(ev_mcshower->at(rmcm.mc_index).Start().Position());
      }
      else if(rmcm.mc_type == 2) {
	std::cout << ev_mctrack->at(rmcm.mc_index).Origin() << " start: " << geoalgo::Point_t(ev_mctrack->at(rmcm.mc_index).Start().Position());
      }
      else if(rmcm.mc_type == 3) {
	std::cout << TrackIDToMCTruth(e, ev_mcparticle->at(rmcm.mc_index).TrackId())->Origin();
      }
      else {
	std::cout << "wah";
      }
      std::cout << "\n";
    }

    std::cout << "\n";
    for(RecoMCMatch const & rmcm : track_matches) {
      std::cout << "reco track: " << rmcm.ratio << " origin: ";
      if(rmcm.mc_type == 1) {
	std::cout << ev_mcshower->at(rmcm.mc_index).Origin();
      }
      else if(rmcm.mc_type == 2) {
	std::cout << ev_mctrack->at(rmcm.mc_index).Origin();
      }
      else if(rmcm.mc_type == 3) {
	std::cout << TrackIDToMCTruth(e, ev_mcparticle->at(rmcm.mc_index).TrackId())->Origin();
      }
      else {
	std::cout << "wah";
      }
      std::cout << "\n";
    }
  
    std::cout << "\n";
    for(sim::MCTrack const & mctr : *ev_mctrack) {
      if(mctr.Origin() == 1) {
	std::cout << "MCTrack Origin: " << mctr.Origin() << " " << mctr.Start().E() << " " << mctr.size() << " " << geoalgo::Point_t(mctr.Start().Position()) << "\n";
      }
    }
    
    for(sim::MCShower const & mcs : *ev_mcshower) {
      if(mcs.Origin() == 1) {
	std::cout << "MCShower Origin: " << mcs.Origin() << " " << mcs.Start().E() << " " << mcs.DetProfile().E() << " " << geoalgo::Point_t(mcs.Start().Position()) << "\n";
      }
    }
    
    std::cout << "\n";
    
  }

  /*

  DetectorObjects const & deto_v = pas.GetDetectorObjects();
  std::vector<ParticleAssociation> const & pa_v = pas.GetAssociations();

  for(size_t const pa_index : pas.GetSelectedAssociations()) {

    double true_total = 0;
    double reco_total = 0;

    for(size_t const object_index : pa_v.at(pa_index).GetObjectIndices()) {
      DetectorObject const & deto = deto_v.GetDetectorObject(object_index);
      if(deto.freco_type == deto_v.fshower_reco_type) {
	RecoMCMatch const & shower_match = shower_matches.at(deto.foriginal_index);
	std::unordered_map<int, double> const & trkide_map = shower_match.trkide_map;
	for(int const trkid : mcparticle_v) 
      }
      else if(deto.freco_type == deto_v.ftrack_reco_type) {
	RecoMCMatch const & track_match = track_matches.at(deto.foriginal_index);
	std::cout << deto.freco_type << " " << track_match.mc_type << "\n";
      }
      else {
	std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nERROR: Unrecognized detector object type\n";
	exit(1);
      }
    }

  }

  */

}



#endif
