////////////////////////////////////////////////////////////////////////
// Class:       LEEPhotonAnalysis
// Plugin Type: analyzer (art v2_05_00)
// File:        LEEPhotonAnalysis_module.cc
//
// Generated at Wed Jun  7 10:49:14 2017 by Robert Murrells using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"
#include "TFile.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "../LLBasicTool/GeoAlgo/GeoVector.h"
#include "../LLBasicTool/GeoAlgo/GeoSphere.h"
#include "../LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "../LLBasicTool/GeoAlgo/GeoCone.h"
#include "../LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "../LLBasicTool/GeoAlgo/GeoAABox.h"


#include "DetectorObjects.h"
#include "ParticleAssociations.h"

#include "VertexBuilder.h"
#include "VertexQuality.h"
#include "FillTreeVariables.h"
#include "RecoMCMatching.h"

#include "RecoMCMatching.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

class LEEPhotonAnalysis : public art::EDAnalyzer {

  double fstart_prox;
  double fshower_prox;
  double fmax_bp_dist;
  double fcpoa_vert_prox; 
  double fcpoa_trackend_prox;
  bool fverbose;

  std::string fmcordata;
  bool fmcrecomatching;
  bool frun_vertex_quality;
  std::string fpot_producer;
  std::string fpfp_producer;
  std::string ftrack_producer;
  std::string ftrack_calo_producer;
  std::string fshower_producer;
  std::string fhit_producer;
  std::string frmcmassociation_producer;
  std::string fopflash_producer;
  std::string fswtrigger_product;

  RecoMCMatching frmcm;

  TTree * fPOTtree;
  int fnumber_of_events;
  double fpot;

  unsigned int fspec_event;

  VertexQuality fvq;
  FillTreeVariables fftv;

public:

  explicit LEEPhotonAnalysis(fhicl::ParameterSet const & p);

  LEEPhotonAnalysis(LEEPhotonAnalysis const &) = delete;
  LEEPhotonAnalysis(LEEPhotonAnalysis &&) = delete;
  LEEPhotonAnalysis & operator = (LEEPhotonAnalysis const &) = delete;
  LEEPhotonAnalysis & operator = (LEEPhotonAnalysis &&) = delete;

  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob();
  void beginSubRun(art::SubRun const & sr);
  void fillwpandora(art::Event const & e,
		    ParticleAssociations & pas);
  void AddTracks(DetectorObjects & detos,
		 art::ValidHandle<std::vector<recob::Track>> const & ev_t,
		 bool const track_original_indices = false);
  void AddShowers(DetectorObjects & detos,
		  art::ValidHandle<std::vector<recob::Shower>> const & ev_s,
		  bool const track_original_indices = false);
  void analyze(art::Event const & e) override;
  void endJob();  

};


LEEPhotonAnalysis::LEEPhotonAnalysis(fhicl::ParameterSet const & p) : 
  EDAnalyzer(p),
  fstart_prox(-1),
  fshower_prox(-1),
  fcpoa_vert_prox(-1),
  fcpoa_trackend_prox(-1),
  fverbose(false),
  fmcrecomatching(false),
  frun_vertex_quality(false),
  fPOTtree(nullptr),
  fnumber_of_events(0),
  fpot(0),
  fspec_event(UINT_MAX),
  fvq(ftrack_producer,
      fshower_producer,
      frmcm) {

  art::ServiceHandle<art::TFileService> tfs;
  fPOTtree = tfs->make<TTree>("get_pot", "");
  fPOTtree->Branch("number_of_events", &fnumber_of_events, "number_of_events/I");

  this->reconfigure(p);

  }


void LEEPhotonAnalysis::reconfigure(fhicl::ParameterSet const & p) {

  fmcordata = p.get<std::string>("mcordata");
  if(fmcordata != "mc" && fmcordata != "data") {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "ERROR: mcordata must be set to \"mc\" or \"data\"\n";
    exit(1);
  }
  p.get_if_present<bool>("mcrecomatching", fmcrecomatching);
  p.get_if_present<bool>("run_vertex_quality", frun_vertex_quality);
  p.get_if_present<std::string>("pot_producer", fpot_producer);
  if(fpot_producer != ""){
    fPOTtree->Branch("pot", &fpot, "pot/D");
  }
  p.get_if_present<std::string>("pfp_producer", fpfp_producer);
  ftrack_producer = p.get<std::string>("track_producer");
  ftrack_calo_producer = p.get<std::string>("track_calo_producer");
  fshower_producer = p.get<std::string>("shower_producer");
  p.get_if_present<std::string>("hit_producer", fhit_producer);
  p.get_if_present<std::string>("rmcmassociation_producer", frmcmassociation_producer);
  if(fmcrecomatching) {
    frmcm.Configure(fhit_producer,
		    ftrack_producer,
		    fshower_producer,
		    frmcmassociation_producer);
  }
  fopflash_producer = p.get<std::string>("opflash_producer");
  fswtrigger_product = p.get<std::string>("trigger_product");

  fftv.SetProducers(fmcordata,
		    ftrack_producer,
		    ftrack_calo_producer,
		    fshower_producer,
		    fhit_producer,
		    fopflash_producer,
		    fswtrigger_product,
		    frmcmassociation_producer,
		    &frmcm);

  fstart_prox = p.get<double>("start_prox");
  fshower_prox = p.get<double>("shower_prox");
  fmax_bp_dist = p.get<double>("max_bp_dist");
  fcpoa_vert_prox = p.get<double>("cpoa_vert_prox");
  fcpoa_trackend_prox = p.get<double>("cpoa_trackend_prox");

  p.get_if_present<bool>("verbose", fverbose);
  frmcm.SetVerbose(fverbose);
  fftv.SetVerbose(fverbose);

  p.get_if_present<unsigned int>("spec_event", fspec_event);

}


void LEEPhotonAnalysis::beginJob() {

  fftv.SetUpTreeBranches();

}


void LEEPhotonAnalysis::beginSubRun(art::SubRun const & sr) {

  if(fpot_producer != ""){
    if(fpot_producer == "generator"){
      fpot += sr.getValidHandle<sumdata::POTSummary>(fpot_producer)->totgoodpot;
    }else{
      art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
      if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
	fpot += potSummaryHandlebnbETOR875->totpot; 
      }
    }
  }


  /*art::Handle< sumdata::POTSummary > potListHandle;
    if(sr.getByLabel("generator",potListHandle))
    fpot +=potListHandle->totpot;

    art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR860;
    if (sr.getByLabel("beamdata","bnbETOR860",potSummaryHandlebnbETOR860)){

    std::cout<<"B0: "<<potSummaryHandlebnbETOR860->totpot<<std::endl;
    }

    art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
    if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
    std::cout<<"B1: "<<potSummaryHandlebnbETOR875->totpot<<std::endl;
    }

    art::Handle<sumdata::POTSummary> potSummaryHandlenumiETORTGT;
    if (sr.getByLabel("beamdata","numiETORTGT",potSummaryHandlenumiETORTGT)){
    std::cout<<"B2: "<< potSummaryHandlenumiETORTGT->totpot<<std::endl;
    }
  */


}


void LEEPhotonAnalysis::fillwpandora(art::Event const & e, ParticleAssociations & pas) {

  art::ValidHandle<std::vector<recob::PFParticle>> const & ev_pfp =e.getValidHandle<std::vector<recob::PFParticle>>(fpfp_producer);

  art::FindManyP<recob::Vertex> PFPToVertex(ev_pfp, e, fpfp_producer);
  art::FindManyP<recob::Shower> PFPToShower(ev_pfp, e, fpfp_producer);
  art::FindManyP<recob::Track> PFPToTrack(ev_pfp, e, fpfp_producer);

  DetectorObjects const & detos = pas.GetDetectorObjects();

  for(size_t i = 0; i < ev_pfp->size(); ++i) {

    recob::PFParticle const & pfp = ev_pfp->at(i);
    if(abs(pfp.PdgCode()) != 12 && abs(pfp.PdgCode()) != 14) continue;
    if(fverbose) std::cout << "PFP: " << i << " Parent pdg: " << pfp.PdgCode() << "\n";

    std::vector<art::Ptr<recob::Vertex>> asso_vertices = PFPToVertex.at(i);
    if(asso_vertices.size() > 1) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one vertex associated with pfp\n";
    }
    else if(asso_vertices.size() == 0) {
      if(fverbose) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: no vertex associated with pfp\n";
      continue;
    }
    std::vector<double> xyz(3);
    asso_vertices.front()->XYZ(&(xyz[0]));
    geoalgo::Point_t const vertex(xyz[0], xyz[1], xyz[2]);
    if(fverbose) std::cout << "Vertex:"  << vertex << "\n";

    std::vector<size_t> associated_indices;

    for(size_t const s : pfp.Daughters()) {

      recob::PFParticle const & daughter = ev_pfp->at(s);

      if(fverbose) std::cout << "\tPFP: " << s << " Daughter pdg: " << daughter.PdgCode() << "\n";

      if(abs(daughter.PdgCode()) == 13) {
	std::vector<art::Ptr<recob::Track>> asso_tracks = PFPToTrack.at(s);
	if(asso_tracks.size() > 1) {
	  std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one track associated with pfp\n";
	}
	else if(asso_tracks.size() == 0) {
	  if(fverbose) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: no track associated with pfp\n";
	  continue;
	}
	else {
	  associated_indices.push_back(detos.GetTrackIndexFromOriginalIndex(asso_tracks.front().key()));
	}
      }

      else if(abs(daughter.PdgCode()) == 11) {
	std::vector<art::Ptr<recob::Shower>> asso_showers = PFPToShower.at(s);
	if(asso_showers.size() > 1) {
	  std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one shower associated with pfp\n";
	}
	else if(asso_showers.size() == 0) {
	  if(fverbose) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: no shower associated with pfp\n";
	  continue;
	}
	else {
	  associated_indices.push_back(detos.GetShowerIndexFromOriginalIndex(asso_showers.front().key()));
	}
      }

      else {
	std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: pfp pdg: " << daughter.PdgCode() << "\n";
      }

    }

    pas.AddAssociation(associated_indices, std::vector<geoalgo::Point_t>(), vertex, 0);

  }

  for(size_t i = 0; i < pas.GetAssociations().size(); ++i) {
    ParticleAssociation const & pa = pas.GetAssociations().at(i);
    if(pa.GetConnections().size() > 0) {
      std::cout << "Warning, pandora filled ParticleAssociation: " << i << " has connections\n";
      exit(1);
    }
  }

  pas.GetShowerAssociations();

}


void LEEPhotonAnalysis::AddTracks(DetectorObjects & detos,
				  art::ValidHandle<std::vector<recob::Track>> const & ev_t,
				  bool const track_original_indices) {

  for(size_t i = 0; i < ev_t->size(); ++i) {

    recob::Track const & t = ev_t->at(i);

    geoalgo::Trajectory traj;
    traj.reserve(t.NumberTrajectoryPoints());
    for(size_t i = 0; i < t.NumberTrajectoryPoints(); ++i)
      traj.push_back(t.LocationAtPoint(i)); 

    detos.AddTrack(i, traj, track_original_indices);

  }

}


void LEEPhotonAnalysis::AddShowers(DetectorObjects & detos,
				   art::ValidHandle<std::vector<recob::Shower>> const & ev_s,
				   bool const track_original_indices) {

  for(size_t i = 0; i < ev_s->size(); ++i) {
    recob::Shower const & s = ev_s->at(i);
    detos.AddShower(i, geoalgo::Cone_t(s.ShowerStart(), s.Direction(), s.Length(), 0), track_original_indices);
  }

}  


void LEEPhotonAnalysis::analyze(art::Event const & e) {

  if(fspec_event != UINT_MAX && e.id().event() != fspec_event) {
    return;
  }

  ++fnumber_of_events;

  //right get a list of all tracks and showers from the named producers (here mostly pandoraNu)
  art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);
  art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);

  if(fverbose)
    std::cout << "Event: " << e.id().event() << std::endl
	      << "=======================================================\n";

  ///////////////// Runs the vertex builder

  VertexBuilder vb;
  vb.SetVerbose(fverbose);

  vb.SetMaximumTrackEndProximity(fstart_prox);
  vb.SetMaximumShowerIP(fshower_prox);
  vb.SetMaximumBackwardsProjectionDist(fmax_bp_dist);
  vb.CPOAToVert(fcpoa_vert_prox);
  vb.SetMaximumTrackEndProx(fcpoa_trackend_prox);

  ParticleAssociations pas;
  pas.SetVerbose(fverbose);

  if(fpfp_producer == "") {
    if(fverbose) std::cout << "Run vertex builder\n";
    AddShowers(pas.GetDetectorObjects(), ev_s);
    AddTracks(pas.GetDetectorObjects(), ev_t);
    vb.Run(pas);
  }
  else {
    if(fverbose) std::cout << "Run pandora\n";
    AddShowers(pas.GetDetectorObjects(), ev_s, true);
    AddTracks(pas.GetDetectorObjects(), ev_t, true);
    fillwpandora(e, pas);
  }

  /////////////////

  if(fmcrecomatching) {
    frmcm.MatchWAssociations(e);
    if(frun_vertex_quality) fvq.RunSig(e, pas);
  }
  
  /////////////////

  if(fverbose){ std::cout << "Fill tree variables\n";}
  fftv.Fill(e, pas);

  if(fspec_event != UINT_MAX && e.id().event() == fspec_event) {
    std::cout << "Found event: " << fspec_event << "\n";
    exit(0);
  }

}


void LEEPhotonAnalysis::endJob() {

  fPOTtree->Fill();

}


DEFINE_ART_MODULE(LEEPhotonAnalysis)
