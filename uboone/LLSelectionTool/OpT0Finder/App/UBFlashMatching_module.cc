////////////////////////////////////////////////////////////////////////
// Class:       UBFlashMatching
// Module Type: producer
// File:        UBFlashMatching_module.cc
//
// Generated at Mon Feb 15 11:26:13 2016 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

//
// Basic stdlib includes
//
#include <string>
#include <iostream>

//
// LArSoft fmwk includes
//
#include "RecoBase/OpFlash.h"
#include "RecoBase/Track.h"

//
// OpT0Finder fmwk includes
//
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/QWeightPoint.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/TimeCompatMatch.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/MaxNPEWindow.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/NPtFilter.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

class UBFlashMatching;

class UBFlashMatching : public art::EDProducer {
public:
  explicit UBFlashMatching(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBFlashMatching(UBFlashMatching const &) = delete;
  UBFlashMatching(UBFlashMatching &&) = delete;
  UBFlashMatching & operator = (UBFlashMatching const &) = delete;
  UBFlashMatching & operator = (UBFlashMatching &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  /// FlashMatcManager instance
  ::flashana::FlashMatchManager _mgr;
  /// LightPath algorithm to convert recob::Track into flashana::QCluster_t
  ::flashana::LightPath _light_path_alg;
  std::string _track_producer_name; ///< Input recob::Track producer name
  std::string _flash_producer_name; ///< Input recob::OpFlash producer name

};


UBFlashMatching::UBFlashMatching(fhicl::ParameterSet const & p)
  : _mgr()
// :
// Initialize member data here.
{
  _track_producer_name = p.get<std::string>("TrackProducer");
  _flash_producer_name = p.get<std::string>("FlashProducer");
  
  //
  // Attach algorithms to be used
  //
  _mgr.SetAlgo( new ::flashana::NPtFilter()           );
  _mgr.SetAlgo( new ::flashana::MaxNPEWindow()        );
  _mgr.SetAlgo( new ::flashana::TimeCompatMatch()     );
  _mgr.SetAlgo( new ::flashana::QWeightPoint()        );
  _mgr.SetAlgo( new ::flashana::PhotonLibHypothesis() );
  // Also attach LightPath instance to be configured via Manager
  _mgr.AddCustomAlgo( &_light_path_alg                );

  //
  // Now configure FlashMatchManager (which configures algorithms as well)
  //
  _mgr.Configure(p);

  // Call appropriate produces<>() functions here.
}

void UBFlashMatching::produce(art::Event & e)
{
  // Define # PMTs here as const (we should retrieve from geo::Geometry for good practice)
  const size_t num_pmts = 32;

  //
  // Steps to be done:
  // -1) Get necessary data products (recob::OpFlash, recob::Track, etc)
  //  0) Provide input to FlashMatchManager: FlashArray_t and QClusterArray_t
  //  1) Run FlashMatchManager & retrieve matches
  //  2) Store data products (anab::FlashMatch and associations)
  //

  //
  // Step -1): Get necessary data products from fmwk
  //
  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(_track_producer_name,trackHandle);
  if(!trackHandle.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Track from " 
	      << _track_producer_name << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector<recob::OpFlash> > flashHandle;
  e.getByLabel(_flash_producer_name,flashHandle);
  if(!flashHandle.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::OpFlash from " 
	      << _flash_producer_name << std::endl;
    throw std::exception();
  }

  //
  //  0) Provide input to FlashMatchManager: FlashArray_t and QClusterArray_t
  //

  //  0-a) FlashArray_t
  for(size_t opflash_index=0; opflash_index < flashHandle->size(); ++opflash_index) {

    // Retrieve individual recob::OpFlash and construct flashana::Flash_t
    auto const& opf = (*flashHandle)[opflash_index];

    ::flashana::Flash_t flash;
    flash.pe_v.resize(num_pmts);
    for(size_t pmt_index=0; pmt_index<num_pmts; ++pmt_index)
      
      flash.pe_v[pmt_index] = opf.PE(pmt_index);
    
    flash.idx   = opflash_index;
    flash.time  = opf.Time();
    flash.x     = 128.;
    flash.x_err = 128.;
    flash.y     = opf.YCenter();
    flash.y_err = opf.YWidth();
    flash.z     = opf.ZCenter();
    flash.z_err = opf.ZWidth();
    
    // Register to a manager
    _mgr.Emplace(std::move(flash));
  }

  //  0-b) QClusterArray_t
  for(size_t track_index=0; track_index < trackHandle->size(); ++track_index) {

    // Retrieve individual recob::OpFlash and construct flashana::Flash_t
    auto const& track = (*trackHandle)[track_index];

    // Construct ::geoalgo::Trajectory (i.e. vector of points) to use for LightPath
    ::geoalgo::Trajectory trj;
    // Set # points same as input track object, and initialize each point as 3D point
    trj.resize(track.NumberTrajectoryPoints(),::geoalgo::Point_t(3,0.));

    // Now loop over points and set actual xyz values
    for(size_t point_index = 0; point_index < trj.size(); ++point_index) {
      
      // Get reference to be modified
      auto&       copy_pt = trj[point_index];
      // Get const reference to get values
      auto const& orig_pt = track.LocationAtPoint(point_index);

      copy_pt[0] = orig_pt[0];
      copy_pt[1] = orig_pt[1];
      copy_pt[2] = orig_pt[2];
    }
    
    auto qcluster = _light_path_alg.FlashHypothesis(trj);

    qcluster.idx = track_index;
    
    // Register to a manager
    _mgr.Emplace(std::move(qcluster));
  } 

  //
  //  1) Run FlashMatchManager & retrieve matches  
  //
  auto match_result_v = _mgr.Match();

  //
  //  2) Store data products (anab::FlashMatch and associations)
  //
  for(auto const& match : match_result_v)
    
    std::cout << "Match result ... "
	      << "Flash ID: " << match.flash_id
	      << " with "
	      << "TPC ID: " << match.tpc_id
	      << " ... Score: " << match.score
	      << std::endl;

}

DEFINE_ART_MODULE(UBFlashMatching)
