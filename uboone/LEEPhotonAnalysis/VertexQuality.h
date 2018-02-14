

#ifndef VERTEXQUALITY_H
#define VERTEXQUALITY_H

#include "ParticleAssociations.h"
#include "RecoMCMatching.h"
#include "FilterSignal.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"


class VertexQuality {

  std::string ftrack_producer;
  std::string fshower_producer;
  
  RecoMCMatching const * frmcm;
  FilterSignal ffs;

  geoalgo::AABox ftpc_volume;

  TTree * fvertex_tree;
  TTree * fvertex_tree_event;
  TTree * fvertex_tree_event_signal;

  double fstart_prox;
  double fshower_prox;
  double fmax_bp_dist;
  double fcpoa_vert_prox; 
  double fcpoa_trackend_prox;

  int freco_vertex_present;
  int fis_nc_delta_rad;
  int fnc_delta_rad_split_shower;
  int ftpc_volume_contained;

  double fdist;
  double fdistx;
  double fdisty;
  double fdistz;

  int ftrue_track_total;
  int ftrue_shower_total;
  int freco_track_total;
  int freco_shower_total;
  int fcorrect_track_total;
  int fcorrect_shower_total;
  
  std::vector<double> ftrack_matching_ratio_v;
  std::vector<int> ftrack_true_pdg_v;
  std::vector<int> ftrack_true_origin_v;

  std::vector<double> fshower_matching_ratio_v;
  std::vector<int> fshower_true_pdg_v;
  std::vector<int> fshower_true_origin_v;

 public:

  VertexQuality(std::string const & track_producer,
		std::string const & shower_producer,
		RecoMCMatching const & rmcm);

  void SetParameters(double const start_prox, 
		     double const shower_prox,
		     double const max_bp_dist,
		     double const cpoa_vert_prox, 
		     double const cpoa_trackend_prox) {

    fstart_prox = start_prox;
    fshower_prox = shower_prox;
    fmax_bp_dist = max_bp_dist;
    fcpoa_vert_prox = cpoa_vert_prox; 
    fcpoa_trackend_prox = cpoa_trackend_prox;    

  }

  void SetupVertexQualityTree();
  void SetupVertexQualityTreeClosest();
  void SetupVertexQualityTreeSignal();

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const & e, int const geant_track_id);

  void GetTrueObjects(art::Event const & e,
		      size_t const mct_index,
		      std::vector<size_t> & mctrack_v,
		      std::vector<size_t> & mcshower_v,
		      std::vector<size_t> & mcparticle_v);
  void GetTrueRecoObjects(art::Event const & e,
			  size_t const mct_index,
			  std::vector<size_t> & track_v,
			  std::vector<size_t> & shower_v);
  void GetTrueObjects(art::Event const & e,
		      std::vector<int> & mcparticle_v);
  double GetTrueTotal(std::vector<int> const & mcparticle_v);

  void Reset();
  void FillTree(art::Event const & e,
		TTree * tree, 
		ParticleAssociations const & pas,
		size_t const pa_index,
		geoalgo::Point_t const & true_nu_vtx,
		std::vector<size_t> const & track_v,
		std::vector<size_t> const & shower_v);
  void RunDist(art::Event const & e,
	       ParticleAssociations const & pas,
	       bool const track_only = false);
  size_t GetNCDeltadRadPhoton(art::Event const & e, size_t const nc_delta_rad_mct_index, size_t const exiting_photon_index, int & mc_type);
  void RunSig(art::Event const & e,
	      ParticleAssociations const & pas,
	      bool const track_only = false);

};


#endif
