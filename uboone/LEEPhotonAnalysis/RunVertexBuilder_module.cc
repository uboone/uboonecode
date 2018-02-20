////////////////////////////////////////////////////////////////////////
// Class:       RunVertexBuilder
// Plugin Type: analyzer (art v2_05_01)
// File:        RunVertexBuilder_module.cc
//
// Generated at Tue Feb  6 13:53:54 2018 by Robert Murrells using cetskelgen
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

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "GetPermutations.h"

#include "VertexBuilder.h"
#include "RecoMCMatching.h"
#include "VertexQuality.h"

class RunVertexBuilder;



class RunVertexBuilder : public art::EDAnalyzer {

  bool fis_nc_delta_rad;
  std::string ftrack_producer;
  std::string fshower_producer;
  bool ftrack_only;
  bool fverbose;

  double fstart_prox;
  double fstart_prox_min;
  double fstart_prox_max;
  double fstart_prox_inc;
  double fshower_prox;
  double fshower_prox_min;
  double fshower_prox_max;
  double fshower_prox_inc;
  double fmax_bp_dist;
  double fmax_bp_dist_min;
  double fmax_bp_dist_max;
  double fmax_bp_dist_inc;
  double fcpoa_vert_prox;
  double fcpoa_vert_prox_min; 
  double fcpoa_vert_prox_max;
  double fcpoa_vert_prox_inc;
  double fcpoa_trackend_prox;
  double fcpoa_trackend_prox_min;
  double fcpoa_trackend_prox_max;
  double fcpoa_trackend_prox_inc;

  int fmin_index;
  int fmax_index;

  double fstart_prox_inc_size;
  double fshower_prox_inc_size;
  double fmax_bp_dist_inc_size;
  double fcpoa_vert_prox_inc_size;
  double fcpoa_trackend_prox_inc_size;

  std::vector<double> fparameter_set;
  std::vector<double> fparameter_min;
  std::vector<double> fparameter_max;
  std::vector<double> fparameter_inc;
  std::vector<double> fparameter_inc_size;

  std::vector< std::vector<double> > fpermutation_v; 

  VertexBuilder fvb;

  std::string fhit_producer;
  std::string frmcmassociation_producer;

  RecoMCMatching frmcm;  

  VertexQuality fvq;

  TTree * fevent_tree;

  int tevent;

  bool ffill_parameter_tree;

  TTree * fparameter_tree;

  double tstart_prox;
  double tshower_prox;
  double tmax_bp_dist;
  double tcpoa_vert_prox;
  double tcpoa_trackend_prox;  

public:

  explicit RunVertexBuilder(fhicl::ParameterSet const & p);
  RunVertexBuilder(RunVertexBuilder const &) = delete;
  RunVertexBuilder(RunVertexBuilder &&) = delete;
  RunVertexBuilder & operator = (RunVertexBuilder const &) = delete;
  RunVertexBuilder & operator = (RunVertexBuilder &&) = delete;

  void AddTracks(DetectorObjects & detos,
		 art::ValidHandle<std::vector<recob::Track>> const & ev_t,
		 bool const track_original_indices = false);
  void AddShowers(DetectorObjects & detos,
		 art::ValidHandle<std::vector<recob::Shower>> const & ev_t,
		 bool const track_original_indices = false);
  void analyze(art::Event const & e) override;

};



RunVertexBuilder::RunVertexBuilder(fhicl::ParameterSet const & p) :
  EDAnalyzer(p),
  ftrack_only(false),
  fverbose(false),
  fstart_prox(-1),
  fstart_prox_min(-1),
  fstart_prox_max(-1),
  fstart_prox_inc(-1),
  fshower_prox(-1),
  fshower_prox_min(-1),
  fshower_prox_max(-1),
  fshower_prox_inc(-1),
  fmax_bp_dist(-1),
  fmax_bp_dist_min(-1),
  fmax_bp_dist_max(-1),
  fmax_bp_dist_inc(-1),
  fcpoa_vert_prox(-1),
  fcpoa_vert_prox_min(-1), 
  fcpoa_vert_prox_max(-1),
  fcpoa_vert_prox_inc(-1),
  fcpoa_trackend_prox(-1),
  fcpoa_trackend_prox_min(-1),
  fcpoa_trackend_prox_max(-1),
  fcpoa_trackend_prox_inc(-1),
  fmin_index(-1),
  fmax_index(-1),
  fvq(ftrack_producer,
      fshower_producer,
      frmcm) {

  fis_nc_delta_rad = p.get<bool>("NCDeltaRadiative");
  ftrack_producer = p.get<std::string>("track_producer");
  fshower_producer = p.get<std::string>("shower_producer");
  p.get_if_present<bool>("track_only", ftrack_only);
  p.get_if_present<bool>("verbose", fverbose);

  p.get_if_present<int>("min_index", fmin_index);
  p.get_if_present<int>("max_index", fmax_index);

  p.get_if_present<double>("start_prox", fstart_prox);
  p.get_if_present<double>("start_prox_min", fstart_prox_min);
  p.get_if_present<double>("start_prox_max", fstart_prox_max);
  p.get_if_present<double>("start_prox_inc", fstart_prox_inc);
  if(fstart_prox == -1 && ((fstart_prox_min == -1 || fstart_prox_max == -1 || fstart_prox_inc == -1) || fmin_index < 0 || fmax_index == -1)) {
    std::cout << "start_prox not set correctly: "
	      << fstart_prox << " "
	      << fstart_prox_min << " "
	      << fstart_prox_max << " "
	      << fstart_prox_inc << "\n";
    exit(1);
  }
  p.get_if_present<double>("shower_prox", fshower_prox);
  p.get_if_present<double>("shower_prox_min", fshower_prox_min);
  p.get_if_present<double>("shower_prox_max", fshower_prox_max);
  p.get_if_present<double>("shower_prox_inc", fshower_prox_inc);
  if(fshower_prox == -1 && ((fshower_prox_min == -1 || fshower_prox_max == -1 || fshower_prox_inc == -1) || fmin_index < 0 || fmax_index == -1)) {
    std::cout << "shower_prox not set correctly: "
	      << fshower_prox << " "
	      << fshower_prox_min << " "
	      << fshower_prox_max << " "
	      << fshower_prox_inc << "\n";
    exit(1);
  }
  p.get_if_present<double>("max_bp_dist", fmax_bp_dist);
  p.get_if_present<double>("max_bp_dist_min", fmax_bp_dist_min);
  p.get_if_present<double>("max_bp_dist_max", fmax_bp_dist_max);
  p.get_if_present<double>("max_bp_dist_inc", fmax_bp_dist_inc);
  if(fmax_bp_dist == -1 && ((fmax_bp_dist_min == -1 || fmax_bp_dist_max == -1 || fmax_bp_dist_inc == -1) || fmin_index < 0 || fmax_index == -1)) {
    std::cout << "max_bp_dist not set correctly: "
	      << fmax_bp_dist << " "
	      << fmax_bp_dist_min << " "
	      << fmax_bp_dist_max << " "
	      << fmax_bp_dist_inc << "\n";
    exit(1);
  }
  p.get_if_present<double>("cpoa_vert_prox", fcpoa_vert_prox);
  p.get_if_present<double>("cpoa_vert_prox_min", fcpoa_vert_prox_min);
  p.get_if_present<double>("cpoa_vert_prox_max", fcpoa_vert_prox_max);
  p.get_if_present<double>("cpoa_vert_prox_inc", fcpoa_vert_prox_inc);
  if(fcpoa_vert_prox == -1 && ((fcpoa_vert_prox_min == -1 || fcpoa_vert_prox_max == -1 || fcpoa_vert_prox_inc == -1) || fmin_index < 0 || fmax_index == -1)) {
    std::cout << "cpoa_vert_prox not set correctly: "
	      << fcpoa_vert_prox << " "
	      << fcpoa_vert_prox_min << " "
	      << fcpoa_vert_prox_max << " "
	      << fcpoa_vert_prox_inc << "\n";
    exit(1);
  }
  p.get_if_present<double>("cpoa_trackend_prox", fcpoa_trackend_prox);
  p.get_if_present<double>("cpoa_trackend_prox_min", fcpoa_trackend_prox_min);
  p.get_if_present<double>("cpoa_trackend_prox_max", fcpoa_trackend_prox_max);
  p.get_if_present<double>("cpoa_trackend_prox_inc", fcpoa_trackend_prox_inc);
  if(fcpoa_trackend_prox == -1 && ((fcpoa_trackend_prox_min == -1 || fcpoa_trackend_prox_max == -1 || fcpoa_trackend_prox_inc == -1) || fmin_index < 0 || fmax_index == -1)) {
    std::cout << "cpoa_trackend_prox not set correctly: "
	      << fcpoa_trackend_prox << " "
	      << fcpoa_trackend_prox_min << " "
	      << fcpoa_trackend_prox_max << " "
	      << fcpoa_trackend_prox_inc << "\n";
    exit(1);
  }

  fstart_prox_inc_size = (fstart_prox_max - fstart_prox_min) / fstart_prox_inc;
  fshower_prox_inc_size = (fshower_prox_max - fshower_prox_min) / fshower_prox_inc;
  fmax_bp_dist_inc_size = (fmax_bp_dist_max - fmax_bp_dist_min) / fmax_bp_dist_inc;
  fcpoa_vert_prox_inc_size = (fcpoa_vert_prox_max - fcpoa_vert_prox_min) / fcpoa_vert_prox_inc;
  fcpoa_trackend_prox_inc_size = (fcpoa_trackend_prox_max - fcpoa_trackend_prox_min) / fcpoa_trackend_prox_inc;

  fparameter_set = {fstart_prox,
		    fshower_prox,
		    fmax_bp_dist,
		    fcpoa_vert_prox,
		    fcpoa_trackend_prox};

  fparameter_min = {fstart_prox_min,
		    fshower_prox_min,
		    fmax_bp_dist_min,
		    fcpoa_vert_prox_min,
		    fcpoa_trackend_prox_min};

  fparameter_max = {fstart_prox_max,
		    fshower_prox_max,
		    fmax_bp_dist_max,
		    fcpoa_vert_prox_max,
		    fcpoa_trackend_prox_max};

  fparameter_inc = {fstart_prox_inc,
		    fshower_prox_inc,
		    fmax_bp_dist_inc,
		    fcpoa_vert_prox_inc,
		    fcpoa_trackend_prox_inc};
  
  fparameter_inc_size = {fstart_prox_inc_size,
			 fshower_prox_inc_size,
			 fmax_bp_dist_inc_size,
			 fcpoa_vert_prox_inc_size,
			 fcpoa_trackend_prox_inc_size};

  if(fmin_index != -1 && fmax_index != -1) {
    GetPermutations(fparameter_set,
		    fparameter_min,
		    fparameter_max,
		    fparameter_inc_size,
		    fpermutation_v);
  }

  if(fmin_index > -1 && fmax_index != -1 && fpermutation_v.size() - 1 < size_t(fmax_index)) {
    std::cout << "max_index: " << fmax_index << " is greater than permutation_v allows: " << fpermutation_v.size() - 1 << "\n";
    exit(1);
  }

  fhit_producer = p.get<std::string>("hit_producer");
  frmcmassociation_producer = p.get<std::string>("rmcmassociation_producer");

  fvb.SetVerbose(fverbose);

  frmcm.Configure(fhit_producer,
		  ftrack_producer,
		  fshower_producer,
		  frmcmassociation_producer);

  art::ServiceHandle< art::TFileService > tfs;
  fevent_tree = tfs->make<TTree>("event_tree", "");  
  fevent_tree->Branch("event", &tevent, "event/I");

  ffill_parameter_tree = true;
  fparameter_tree = tfs->make<TTree>("parameter_tree", "");  
  fparameter_tree->Branch("start_prox", &tstart_prox, "start_prox/D");
  fparameter_tree->Branch("shower_prox", &tshower_prox, "shower_prox/D");
  fparameter_tree->Branch("max_bp_dist", &tmax_bp_dist, "max_bp_dist/D");  
  fparameter_tree->Branch("cpoa_vert_prox", &tcpoa_vert_prox, "cpoa_vert_prox/D");
  fparameter_tree->Branch("cpoa_trackend_prox", &tcpoa_trackend_prox, "cpoa_trackend_prox/D");

}



void RunVertexBuilder::AddTracks(DetectorObjects & detos,
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



void RunVertexBuilder::AddShowers(DetectorObjects & detos,
				  art::ValidHandle<std::vector<recob::Shower>> const & ev_s,
				  bool const track_original_indices) {

  for(size_t i = 0; i < ev_s->size(); ++i) {
    recob::Shower const & s = ev_s->at(i);
    detos.AddShower(i, geoalgo::Cone_t(s.ShowerStart(), s.Direction(), s.Length(), 0), track_original_indices);
  }

}  



void RunVertexBuilder::analyze(art::Event const & e) {

  art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producer);
  art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);

  tevent = e.id().event();
  fevent_tree->Fill();

  if(fmin_index < 0 || fmax_index == -1) {

    ParticleAssociations pas;
    AddShowers(pas.GetDetectorObjects(), ev_s);
    AddTracks(pas.GetDetectorObjects(), ev_t);
    if(fverbose) std::cout << "Run VB\n";
    
    fvb.SetMaximumTrackEndProximity(fstart_prox);
    fvb.SetMaximumShowerIP(fshower_prox);
    fvb.SetMaximumBackwardsProjectionDist(fmax_bp_dist);
    fvb.CPOAToVert(fcpoa_vert_prox);
    fvb.SetMaximumTrackEndProx(fcpoa_trackend_prox);  
    
    if(ftrack_only) {
      fvb.SetDetectorObjects(&pas.GetDetectorObjects());
      fvb.AssociateTracks(pas);
    }
    else fvb.Run(pas);
    frmcm.MatchWAssociations(e);
    if(fverbose) std::cout << "Run VQ\n";
    fvq.SetParameters(fstart_prox, fshower_prox, fmax_bp_dist, fcpoa_vert_prox, fcpoa_trackend_prox);
    if(fis_nc_delta_rad) fvq.RunSig(e, pas, ftrack_only);
    else fvq.RunDist(e, pas, ftrack_only);

    tstart_prox = fstart_prox;
    tshower_prox = fshower_prox;
    tmax_bp_dist = fmax_bp_dist;
    tcpoa_vert_prox = fcpoa_vert_prox;
    tcpoa_trackend_prox = fcpoa_trackend_prox;
    fparameter_tree->Fill();
    
  }

  else {

    for(size_t i = fmin_index; i <= size_t(fmax_index); ++i) {

      ParticleAssociations pas;
      AddShowers(pas.GetDetectorObjects(), ev_s);
      AddTracks(pas.GetDetectorObjects(), ev_t);
      if(fverbose) std::cout << "Run VB\n";

      std::vector<double> const & parameters = fpermutation_v.at(i);
      double const start_prox = parameters.at(0);
      double const shower_prox = parameters.at(1);
      double const max_bp_dist = parameters.at(2);
      double const cpoa_vert_prox = parameters.at(3);
      double const cpoa_trackend_prox = parameters.at(4);

      fvb.SetMaximumTrackEndProximity(start_prox);
      fvb.SetMaximumShowerIP(shower_prox);
      fvb.SetMaximumBackwardsProjectionDist(max_bp_dist);
      fvb.CPOAToVert(cpoa_vert_prox);
      fvb.SetMaximumTrackEndProx(cpoa_trackend_prox);  
  
      if(ftrack_only) {
	fvb.SetDetectorObjects(&pas.GetDetectorObjects());
	fvb.AssociateTracks(pas);
      }
      else fvb.Run(pas);
      frmcm.MatchWAssociations(e);
      if(fverbose) std::cout << "Run VQ\n";
      fvq.SetParameters(start_prox, shower_prox, max_bp_dist, cpoa_vert_prox, cpoa_trackend_prox);
      if(fis_nc_delta_rad) fvq.RunSig(e, pas, ftrack_only);
      else fvq.RunDist(e, pas, ftrack_only);

      if(ffill_parameter_tree) {
	tstart_prox = start_prox;
	tshower_prox = shower_prox;
	tmax_bp_dist = max_bp_dist;
	tcpoa_vert_prox = cpoa_vert_prox;
	tcpoa_trackend_prox = cpoa_trackend_prox;
	fparameter_tree->Fill();   
      }

      if(fverbose) {
	std::cout << "cpoa_trackend_prox: " << cpoa_trackend_prox << "\n"
		  << "cpoa_vert_prox: " << cpoa_vert_prox << "\n"
		  << "max_bp_dist: " << max_bp_dist << "\n"
		  << "shower_prox: " << shower_prox << "\n"
		  << "start_prox: " << start_prox << "\n";
      }

    }

    ffill_parameter_tree = false;

  }

}



DEFINE_ART_MODULE(RunVertexBuilder)
