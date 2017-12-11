


#ifndef RECOMCMATCHING_H
#define RECOMCMATCHING_H



#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"

#include "SimchInfo.h"



struct RecoMCMatch {

  int mc_type;
  size_t mc_index;
  double ratio;
  
  RecoMCMatch(int imc_type, size_t imc_index, double iratio);
  
};



class RecoMCMatching {

  std::string fmatch_type;
  std::string fsimch_producer;
  std::string fhit_producer;
  std::string ftrack_producer;
  std::string fshower_producer;

  bool fconsider_mcparticles;

  TTree * fhit_tree;
  int fhit_from_reco_track;
  double fhit_time;
  double fhit_tdc;

  TTree * fsimch_tree;
  int fsimch_tdc;
  
  SimchInfo simchi;

  std::vector<std::vector<double>> fmcq_vv;
  std::vector<RecoMCMatch> ftrack_matches;
  std::vector<RecoMCMatch> fshower_matches;



 public:



  RecoMCMatching();

  

  void Configure(std::string const & match_type,
		 std::string const & simch_producer,
		 std::string const & hit_producer,
		 std::string const & track_producer,
		 std::string const & shower_producer) {

    fmatch_type = match_type;
    fsimch_producer = simch_producer;
    fhit_producer = hit_producer;
    ftrack_producer = track_producer;
    fshower_producer = shower_producer;

  }

  void ConsiderMCParticles(bool const consider_mcparticles = true) {fconsider_mcparticles = consider_mcparticles;}

  void FillHitTree(art::Event const & e);
  void FillSimchTree(art::Event const & e);

  void Match(art::Event const & e,
	     std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void MatchMCBTAlg(art::Event const & e,
		    std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void MatchSimchInfo(art::Event const & e,
		      std::vector<std::vector<art::Ptr<recob::Hit>>> const & reco_to_hit_v);
  void FillRatio(std::vector<RecoMCMatch> & ratio_v,
		 size_t const ev_mctrack_size,
		 size_t const ev_mcshower_size,
		 size_t const ratio_index,
		 size_t const best_mc_index,
		 double const ratio,
		 bool const last);
  void MatchAll(art::Event const & e);
  void CoutMatches(art::Event const & e);

  std::vector<RecoMCMatch> const & GetTrackMatches() const {
    return ftrack_matches;
  }
  std::vector<RecoMCMatch> const & GetShowerMatches() const {
    return fshower_matches;
  }      
  std::vector< std::vector<double> > const & GetMatchVV() const {
    return fmcq_vv;
  }



};



#endif
