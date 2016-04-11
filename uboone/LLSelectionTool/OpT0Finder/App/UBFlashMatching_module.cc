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
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <memory>
//
// ROOT fmwk includes
//
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"


//
// Basic stdlib includes
//
#include <string>
#include <iostream>

//
// LArSoft fmwk includes
//
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/AnalysisBase/FlashMatch.h"
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

//
// Drift velocity includes
//
//#include "lardata/Utilities/LArPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "larcore//CoreUtils/ServiceUtil.h"

// Truth matching includes
//
#include "larsim/MCCheater/BackTracker.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTrajectory.h"

class UBFlashMatching;

//using namespace anab;

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
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob();

private:

  /// FlashMatcManager instance
  ::flashana::FlashMatchManager _mgr;
  /// LightPath algorithm to convert recob::Track into flashana::QCluster_t
  ::flashana::LightPath _light_path_alg;
  std::string _track_producer_name; ///< Input recob::Track producer name
  std::string _flash_producer_name; ///< Input recob::OpFlash producer name
  std::string fPFPModuleLabel;      ///< Input recob::PFParticle producer name
 
  //DetectorProperties
  //const util::LArProperties& detp;

 
  ///Histogram Initializations
  TH1D* fTrackIDCodeHist;
  TH1D* fTrackPhiHist;

  Double_t fTrackCharge;
  Double_t fTrackLight;
  Double_t fTrackPosY;
  Double_t fTrackPosZ;
  Double_t fTrackTrueTime;
  Double_t fTrackRecoTime;
  Double_t fFlashTime;
  Double_t fFlashTimeWidth;
  Double_t fFlashAbsTime;
  Double_t fFlashOnBeamTime;
  Int_t    fTrackMatched;
  Double_t fMatchScore;

  TTree *fMatchTree;

  Int_t fNPandoraTrees;
  Int_t fNTracks;
  Int_t fNFlashes;
  Int_t fNFilteredFlashes;
  Int_t fMultiMatch;
  Double_t fMMFlash_Light;
  Double_t fMMFlash_Time;
  Double_t fMMFlash_TimeWidth;
  Double_t fMMFlash_AbsTime;
  Double_t fMMFlash_OnBeamTime;
  //Double_t fMM_DeltaT;
  //Double_t fMM_DeltaR;

  TTree *fEventTree;

  Double_t fFlashLight;
  Double_t fFlashyTime;
  Int_t    fFlashNMatches;

  TTree *fFlashTree;

  //TH2D* fTrackIDvsFlashMatchScoreHist;
  int fTrackID;
  int fTrackPhi;
  std::string fSpillName;

};


UBFlashMatching::UBFlashMatching(fhicl::ParameterSet const & p)
  : EDProducer()
  , _mgr()
  //, detp(art::ServiceHandle<util::LArProperties>()) 
   
// :
// Initialize member data here.
{
  
  this->reconfigure(p); 
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
  produces< std::vector<anab::FlashMatch> >(fSpillName);
  produces< art::Assns <recob::Track, anab::FlashMatch> >(fSpillName);
  //produces< art::Assns <recob::OpFlash, anab::FlashMatch > >(fSpillName);
  // Call appropriate produces<>() functions here.
}

void UBFlashMatching::reconfigure(fhicl::ParameterSet const& p)
{
  _track_producer_name = p.get<std::string>("TrackProducer");
  _flash_producer_name = p.get<std::string>("FlashProducer");
  fPFPModuleLabel      = p.get<std::string>("PFPModulelabel","pandoraCosmic");
   fSpillName.clear();
   size_t pos = _track_producer_name.find(":");
   if( pos!=std::string::npos ) {
      fSpillName = _track_producer_name.substr( pos+1 );
      _track_producer_name = _track_producer_name.substr( 0, pos );
    }

  return;
}
void UBFlashMatching::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
    
    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fTrackIDCodeHist        = tfs->make<TH1D>("trackIDcodes",";Track ID Code;", 10, -10, 10);
    fTrackPhiHist           = tfs->make<TH1D>("trackPhi",";Track Phi;"        , 10, -5, 5);

    fMatchTree = tfs->make<TTree>("match_tree", "match_tree");
    fEventTree = tfs->make<TTree>("event_tree", "event_tree");
    fFlashTree = tfs->make<TTree>("flash_tree", "flash_tree");

    fMatchTree->Branch("TrackCharge"    , &fTrackCharge     , "TrackCharge/D"     );
    fMatchTree->Branch("TrackLight"     , &fTrackLight      , "TrackLight/D"      );
    fMatchTree->Branch("TrackPosY"      , &fTrackPosY       , "TrackPosY/D"       );
    fMatchTree->Branch("TrackPosZ"      , &fTrackPosZ       , "TrackPosZ/D"       );
    fMatchTree->Branch("TrackTrueTime"  , &fTrackTrueTime   , "TrackTrueTime/D"   );
    fMatchTree->Branch("TrackRecoTime"  , &fTrackRecoTime   , "TrackRecoTime/D"   );
    fMatchTree->Branch("FlashTime"      , &fFlashTime       , "FlashTime/D"       );
    fMatchTree->Branch("FlashTimeWidth" , &fFlashTimeWidth  , "FlashTimeWidth/D"  );
    fMatchTree->Branch("FlashAbsTime"   , &fFlashAbsTime    , "FlashAbsTime/D"    );
    fMatchTree->Branch("FlashOnBeamTime", &fFlashOnBeamTime , "FlashOnBeamTime/D" );

    fMatchTree->Branch("TrackMatched"  , &fTrackMatched  , "TrackMatched/I"   );
    fMatchTree->Branch("MatchScore"    , &fMatchScore    , "MatchScore/D"     );

    fEventTree->Branch("NPandoraTrees"       , &fNPandoraTrees      , "NPandoraTrees/I"     );
    fEventTree->Branch("NTracks"             , &fNTracks            , "NTracks/I"           );
    fEventTree->Branch("NFlashes"            , &fNFlashes           , "NFlashes/I"          );
    fEventTree->Branch("NFilteredFlashes"    , &fNFilteredFlashes   , "NFilteredFlashes/I"  );
    fEventTree->Branch("MultiMatch"          , &fMultiMatch         , "MultiMatch/I"        );
    fEventTree->Branch("MMFlash_Light"       , &fMMFlash_Light      , "MMFlash_Light/D"     );
    fEventTree->Branch("MMFlash_Time"        , &fMMFlash_Time       , "MMFlash_Time/D"      );
    fEventTree->Branch("MMFlash_TimeWidth"   , &fMMFlash_TimeWidth  , "MMFlash_TimeWidth/D" );
    fEventTree->Branch("MMFlash_AbsTime"     , &fMMFlash_AbsTime    , "MMFlash_AbsTime/D"   );
    fEventTree->Branch("MMFlash_OnBeamTime"  , &fMMFlash_OnBeamTime , "MMFlash_OnBeamTime/D");

    fFlashTree->Branch("FlashLight"   , &fFlashLight    , "FlashLight/D"   );
    fFlashTree->Branch("FlashTime"    , &fFlashyTime    , "FlashTime/D"    );

}

void UBFlashMatching::produce(art::Event & e)
{
  std::cout << "NEW EVENT" << std::endl;
  _mgr.Reset();
  // Define # PMTs here as const (we should retrieve from geo::Geometry for good practice)
  const size_t num_pmts = 32;

  //
  //pointer to put onto event
  //

  //art::Ptr<anab::FlashMatch>  flashmatch ( new anab::FlashMatch);
  std::unique_ptr<std::vector<anab::FlashMatch>> flashmatchtrack ( new std::vector<anab::FlashMatch>);
  std::unique_ptr<std::vector<anab::FlashMatch>> flashmatchopflash ( new std::vector<anab::FlashMatch>);
  //std::unique_ptr<std::vector<anab::FlashMatch>> flashmatch ( new std::vector<anab::FlashMatch>);
  std::unique_ptr<art::Assns<recob::Track, anab::FlashMatch>> flashTrackAssociations( new art::Assns<recob::Track, anab::FlashMatch>);
  std::unique_ptr<art::Assns<recob::OpFlash, anab::FlashMatch>> flashOpFlashAssociations( new art::Assns<recob::OpFlash, anab::FlashMatch >);

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
  e.getByLabel(_track_producer_name,fSpillName, trackHandle);
  if(!trackHandle.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Track from " 
	      << _track_producer_name << std::endl;
    throw std::exception();
  }

  art::FindManyP<recob::Hit> trackHitAssns(trackHandle, e, _track_producer_name);
  if(!trackHitAssns.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Hit from " 
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
  else
  {
    std::cout << "flashHandle has size " << flashHandle->size() << std::endl;
  }
 //-----------------------------PFParticle Extraction-----------------------------------------// 
   art::Handle< std::vector<recob::PFParticle> > pfpVecHandle;
   e.getByLabel(fPFPModuleLabel,  pfpVecHandle);
   if(!pfpVecHandle.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::PFParticle from " 
	      << fPFPModuleLabel << std::endl;
    throw std::exception();
  }

  art::FindManyP<recob::Track> trackPFPAssns(pfpVecHandle, e , _track_producer_name); 
  if(!trackPFPAssns.isValid()) {
    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Track from " 
	      << fPFPModuleLabel << std::endl;
    throw std::exception();
  }

  // Get a PFParticle-to-vertex map.
  lar_pandora::VertexVector allPfParticleVertices;
  lar_pandora::PFParticlesToVertices pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(e, fPFPModuleLabel, allPfParticleVertices, pfParticleToVertexMap);

  lar_pandora::PFParticleVector pfparticlelist;
  lar_pandora::PFParticlesToClusters pfParticleToClusterMap;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, fPFPModuleLabel, pfparticlelist, pfParticleToClusterMap);
  
  lar_pandora::TrackVector allPfParticleTracks;
  lar_pandora::PFParticlesToTracks pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, fPFPModuleLabel, allPfParticleTracks, pfParticleToTrackMap);
  
//  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(pfparticlelist, pfParticleToClusterMap,  &clustersToHits, PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, 1);
  std::cout<<"#allPFParticleTracks: "<<allPfParticleTracks.size()<<std::endl;  

  size_t NPFParticles = pfparticlelist.size();
  std::cout<<"#PFParticles: "<<NPFParticles<<std::endl;
  std::vector< int> selfID;
  std::vector <int> isPrimary;
  std::vector <int> numDaughters;
  std::vector <int> parentID;
  std::vector <int> daughterIDs;
  std::vector <int> pdgCode;
  std::vector <int> isNeutrino;
  std::vector <int> isTrack;
  std::vector <int> trackID;
  std::vector <int> vertexID;

for (size_t i = 0; i < NPFParticles; ++i){
     // std::cout<<"First for loop"<<std::endl;
      selfID.push_back(pfparticlelist[i]->Self());
     // std::cout<<"SelfID Set"<<std::endl;
      isPrimary.push_back((Short_t)pfparticlelist[i]->IsPrimary());
     // std::cout<<"Primary Set"<<std::endl;
      numDaughters.push_back(pfparticlelist[i]->NumDaughters());
     // std::cout<<"Daughters Set"<<std::endl;
      parentID.push_back(pfparticlelist[i]->Parent());
      //std::cout<<"Parent Set"<<std::endl;
      pdgCode.push_back(pfparticlelist[i]->PdgCode());
     // std::cout<<"PDG Set"<<std::endl;
         
      
      std::cout<<"particleListinfo: "<<" SelfID[i] "<<selfID.at(i)<<" isPrimary[i]: "<<isPrimary[i]<<" numDaughters[i]: "<<numDaughters[i]<<" parentID[i]: "<<parentID[i]<<" pdgCode: "<<pdgCode[i]<<std::endl;

      // Set the daughter IDs.
      std::vector<size_t> daughterIDs = pfparticlelist[i]->Daughters();
      std::cout<<"right after daughter loop: NumDaughters: "<<daughterIDs.size()<<std::endl; 
      
     /* for (size_t j = 0; j < daughterIDs.size(); ++j)
        daughterIDs.push_back(daughterIDs[j]);*/

      auto vertexMapIter = pfParticleToVertexMap.find(pfparticlelist[i]);
	std::cout<<"above vertexMapIter if statement:"<<std::endl;
      if (vertexMapIter != pfParticleToVertexMap.end()) {
          lar_pandora::VertexVector pfParticleVertices = vertexMapIter->second;
          
          if (pfParticleVertices.size() == 1)
            vertexID.push_back(pfParticleVertices.at(0)->ID());
 	    
          else
            std::cerr << "Warning: there was more than one vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
      }
      else
        std::cerr << "Warning: there was no vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;

 
  if (lar_pandora::LArPandoraHelper::IsTrack(pfparticlelist[i])){
     isTrack.push_back(1);

	std::cout<<"above trackMapIter if statement:"<<std::endl;
   auto trackMapIter = pfParticleToTrackMap.find(pfparticlelist[i]);
    if (trackMapIter != pfParticleToTrackMap.end()) {
       lar_pandora::TrackVector pfParticleTracks = trackMapIter->second;
            
       if (pfParticleTracks.size() == 1)
         trackID.push_back(pfParticleTracks.at(0)->ID());
       else
          std::cerr << "Warning: there was more than one track found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
     }
     else
       std::cerr << "Warning: there was no track found for track-like PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
   }
   else
       isTrack.push_back(0);
}
 //----------------------------End PFParticle Extraction---------------------------------------// 

  //
  //  0) Provide input to FlashMatchManager: FlashArray_t and QClusterArray_t
  //

  //  0-a) FlashArray_t
  std::vector<flashana::Flash_t> opflashVec;
  int count = 1;

  art::ServiceHandle<cheat::BackTracker> bt;

  int nfiltered = 0;
  for(size_t opflash_index=0; opflash_index < flashHandle->size(); ++opflash_index) {

    // Retrieve individual recob::OpFlash and construct flashana::Flash_t
    auto const& opf = (*flashHandle)[opflash_index];
  //  opflashPtrVec.push_back(opf);
    ::flashana::Flash_t flash;
    flash.pe_v.resize(num_pmts);
    for(size_t pmt_index=0; pmt_index<num_pmts; ++pmt_index)
    {
      flash.pe_v[pmt_index] = opf.PE(pmt_index);
    }
    
    flash.idx   = opflash_index;
    flash.time  = opf.Time();
    flash.x     = 128.;
    flash.x_err = 128.;
    flash.y     = opf.YCenter();
    flash.y_err = opf.YWidth();
    flash.z     = opf.ZCenter();
    flash.z_err = opf.ZWidth();

    ::flashana::Flash_t flash_clone = flash;
    opflashVec.push_back(flash_clone);
  
    fFlashLight = opf.TotalPE();
    fFlashyTime = opf.Time();
    fFlashTimeWidth  = opf.TimeWidth();
    fFlashAbsTime    = opf.AbsTime();
    fFlashOnBeamTime = opf.OnBeamTime();

    fFlashTree->Fill();
    if (fFlashLight > 5) nfiltered++;

    // Register to a manager
    _mgr.Emplace(std::move(flash));
    count++;
  }
  //  size_t opflashsize = opf.size();

  //  0-b) QClusterArray_t
   std::cout<<"finished opflash stuff"<<std::endl;
 
std::cout<<"#PFParticles: "<<pfpVecHandle->size()<<std::endl;
/*-------------------------------FLASHMATCH VIA primary PFParticle in Event------------------------- 
int pfp_index =0;
for (size_t i = 0; i < NPFParticles; ++i)
{
  if (pfparticlelist[i]->IsPrimary()==1) 
  {
    if(pfParticleToTrackMap.count(pfparticlelist[i])==0) { std::cout<<"NO TRACKS made for pfparticle#: "<<i<<std::endl; continue;}
    auto trackVec = pfParticleToTrackMap.find(pfparticlelist[i])->second;
    for( auto const& track: trackVec)
    { 
      
      fTrackID = track->ID(); 
      fTrackIDCodeHist->Fill(fTrackID);
      fTrackPhi = track->Phi();
      fTrackPhiHist->Fill(fTrackPhi);
    
      // Construct ::geoalgo::Trajectory (i.e. vector of points) to use for LightPath
      ::geoalgo::Trajectory trj;
      // Set # points same as input track object, and initialize each point as 3D point
      trj.resize(track->NumberTrajectoryPoints(),::geoalgo::Point_t(3,0.));
      // Now loop over points and set actual xyz values
      for(size_t point_index = 0; point_index < trj.size(); ++point_index) 
      {
        // Get reference to be modified
        auto&       copy_pt = trj[point_index];
        // Get const reference to get values
        auto const& orig_pt = track->LocationAtPoint(point_index);
        
        copy_pt[0] = orig_pt[0];
        copy_pt[1] = orig_pt[1];
        copy_pt[2] = orig_pt[2];
	std::cout<<"Trjpnts: "<< copy_pt[0]<<", "<<copy_pt[1]<<", "<<copy_pt[2]<<std::endl;
      }
      
      auto qcluster = _light_path_alg.FlashHypothesis(trj);
      qcluster.idx = pfp_index;
      ++pfp_index; 

      // Register to a manager
      _mgr.Emplace(std::move(qcluster));
     }
  }
}
*/

/*-------------------------------FLASHMATCH VIA Summing all tracks associated to all PFParticles in Event-------------------------*/ 
int pfp_index =0;
::flashana::QCluster_t summed_cluster;
for (size_t i = 0; i < NPFParticles; ++i)
{
  if(pfParticleToTrackMap.count(pfparticlelist[i])==0) { std::cout<<"NO TRACKS made for pfparticle#: "<<i<<std::endl; continue;}
  auto trackVec = pfParticleToTrackMap.find(pfparticlelist[i])->second;
  for( auto const& track: trackVec)
  { 
      
    fTrackID = track->ID(); 
    fTrackIDCodeHist->Fill(fTrackID);
    fTrackPhi = track->Phi();
    fTrackPhiHist->Fill(fTrackPhi);
    
    // Construct ::geoalgo::Trajectory (i.e. vector of points) to use for LightPath
    ::geoalgo::Trajectory trj;
    // Set # points same as input track object, and initialize each point as 3D point
    trj.resize(track->NumberTrajectoryPoints(),::geoalgo::Point_t(3,0.));

    // Now loop over points and set actual xyz values
    for(size_t point_index = 0; point_index < trj.size(); ++point_index) 
    {
      
      // Get reference to be modified
      auto&       copy_pt = trj[point_index];
      // Get const reference to get values
      auto const& orig_pt = track->LocationAtPoint(point_index);

      copy_pt[0] = orig_pt[0];
      copy_pt[1] = orig_pt[1];
      copy_pt[2] = orig_pt[2];
    }
    
    auto qcluster = _light_path_alg.FlashHypothesis(trj);
    summed_cluster += qcluster;
  }
    summed_cluster.idx = pfp_index;
    
    // Register to a manager
    _mgr.Emplace(std::move(summed_cluster));
}


//-------------------------------commenting out old Track Match Way for PFParticle Implementation------------------------- 
/* for(size_t track_index=0; track_index < trackHandle->size(); ++track_index) {


    // Retrieve individual recob::Track and construct flashana::Flash_t
    auto const& track = (*trackHandle)[track_index];
    // trackPtrVec.push_back(track);
    fTrackID = track.ID(); 
    fTrackIDCodeHist->Fill(fTrackID);
    fTrackPhi = track.Phi();
    fTrackPhiHist->Fill(fTrackPhi);
    
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
*/
//-------------------------------------------------------------------------------------------------------------------------
  

  //
  //  1) Run FlashMatchManager & retrieve matches  
  //
  auto match_result_v = _mgr.Match();

  //
  //  2) Store data products (anab::FlashMatch and associations)
  //
  
  bool inbeam= true;	//should get this info from the fcl parameters. 
 // const art::Ptr<recob::Track>  trackPtr;
 // size_t match_track_index = 0;
 // size_t match_flash_index = 0;

  detinfo::DetectorProperties const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detinfo::DetectorClocks const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();
  //auto const* detprop = lar::providerFrom< detinfo::DetectorProperties >();

 // auto const* detp = art::ServiceHandle<util::LArProperties>();
//  art::ServiceHandle<util::LArProperties> detp;
  
  
int NumPrimaries=0;
for(unsigned int i = 0; i <isPrimary.size(); i++)
    if(isPrimary[i]==1) ++NumPrimaries;
std::cout<<"NumberOfPrimaries: "<<NumPrimaries<<std::endl;

const double driftVelocity = detprop->DriftVelocity( detprop->Efield(), detprop->Temperature() );
const double triggertime   = detclock->TriggerTime();
const double triggeroffsettpc   = detclock->TriggerOffsetTPC();

const double beamgate = detclock->BeamGateTime();

std::cout<<"BeamGateTime:  "	<<beamgate<<std::endl;
std::cout<<"Trigger Time: "	<<triggertime<<std::endl;
std::cout<<"Trigger OffsetTPC: "<<triggeroffsettpc<<std::endl;
std::vector<::flashana::FlashMatch_t> match_v;

fNPandoraTrees = NPFParticles;
fNTracks = allPfParticleTracks.size();
std::cout << "Counting flashes" << std::endl;
fNFlashes = flashHandle->size();
fNFilteredFlashes = nfiltered;
std::cout << "Filling event tree" << std::endl;

const int n_flashes = opflashVec.size(); 
int n_matches[n_flashes];
for (int i = 0; i < n_flashes; i++) n_matches[i] = 0;

for (size_t i = 0; i < NPFParticles; ++i){
   if (pfparticlelist[i]->IsPrimary()==1) {
    if(pfParticleToTrackMap.count(pfparticlelist[i])==0) { std::cout<<"NO TRACKS made for pfparticle#: "<<i<<std::endl; continue;}
    auto trackVec = pfParticleToTrackMap.find(pfparticlelist[i])->second;
    std::cout<<"Size of trackVec: "<<trackVec.size()<<std::endl;
    for( auto const& track: trackVec)
    { 
    std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
    

    float charge=0;

std::cout<<"clearing charge"<<std::endl;    
/*
std::cout<<"track pointer size: "<< trackHitAssns.at(track.key()).size()<<std::endl;
    std::vector<art::Ptr<recob::Hit>> trackHitVec = trackHitAssns.at(track.key());
std::cout<<"creating recobHitVec: "<<trackHitVec.size()<<std::endl;    
    for(size_t hit_index=0; hit_index<trackHitVec.size(); ++hit_index)
    {
    	charge += trackHitVec.at(hit_index)->Integral();
    }
    std::cout<<"Hit Charge: "<<charge<<std::endl;
 */   
    TVector3 startpos = track->LocationAtPoint(0);
    TVector3 endpos = track->LocationAtPoint(track->NumberTrajectoryPoints() - 1);
    startpos += endpos;
    startpos *= 0.5;
    std::cout << "Average position: " << startpos.Y() << ", " << startpos.Z() << std::endl;
    
    ::flashana::FlashMatch_t match; 
    for(size_t match_index=0; match_index < match_result_v.size(); ++match_index) 
    {
      match = match_result_v.at(match_index);
      if((int)match.tpc_id==track->ID()) 
      {
        break; 
      }
    }

    std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
    if (match.tpc_id==::flashana::kINVALID_ID) 
    {
      std::cout<<"INVALID TPC_ID"<<std::endl;
      fTrackCharge = -100;
      fTrackLight = -100;
      fTrackPosY = -100;
      fTrackPosZ = -100;
      fTrackTrueTime = -100;
      fTrackRecoTime = -100;
      fFlashTime = -100;
      fTrackMatched = 0;
      fMatchScore = -100;

      fMatchTree->Fill();
    }
    else
    {
      fTrackMatched = 1;
      fMatchScore = match.score;
     std::cout<<"--------------------------------------------TESTING ELSE STATEMENT---------------------------"<<std::endl;
      double light = opflashVec.at(match.flash_id).TotalPE();
      fTrackLight = light;
      fTrackCharge = charge;
      fTrackPosY = startpos.Y();
      fTrackPosZ = startpos.Z();

      double driftpos = track->LocationAtPoint(0).X();
      double drifttime = (driftpos/10.)/driftVelocity; //converting from mm to cm
      fTrackRecoTime = drifttime;

      fTrackTrueTime = -100;
      std::cout << "BACKTRACKER AWAY" << std::endl;
      std::cout << "Track ID = " << track->ID() << ", track length = " << track->Length() << std::endl;
      const simb::MCParticle *true_particle = bt->TrackIDToParticle(track->ID()+1);
      std::cout << "BACKTRACKER YOU DID GOOD" << std::endl;
      if (true_particle != NULL)
      {
        double truetime = true_particle->T();
        std::cout << "BACKTRACKER THAT WAS A GOOD PARTICLE YOU FOUND, IT HAS TIME = " << truetime << std::endl;
        const double trigOffset = detclock->G4ToElecTime(truetime) - detclock->TriggerTime();
        std::cout<<"Calculated Trig Offset: "<<trigOffset<<" G4ToElecTim: "<<detclock->G4ToElecTime(truetime)<<" TriggerTime: "<<detclock->TriggerTime()<<" TrueTime: "<<truetime<<" BeamGate: "<<beamgate<<std::endl;
        fTrackTrueTime = truetime - (1000*trigOffset);
      }
      else
        std::cout << "BACKTRACKER THAT WAS NOT A GOOD PARTICLE THOUGH PLEASE TRY HARDER" << std::endl;


      double flashtime = opflashVec.at(match.flash_id).time;
      fFlashTime = flashtime;
      n_matches[match.flash_id]++;

      fMatchTree->Fill();
      std::cout<<"-------------------------------FINISHED FILLING TTREE---------------------------"<<std::endl;
      
    }
    anab::FlashMatch Flash((double)match.score, (int)match.flash_id, (int)match.tpc_id, (bool)inbeam);
    
    flashmatchtrack->push_back(Flash);
    util::CreateAssn(*this, e, *flashmatchtrack, track,*flashTrackAssociations,fSpillName);
  }
}
}

int max_matches = 0;
int max_match_index = -1;
for (int i = 0; i < n_flashes; i++)
{
  if (n_matches[i] > max_matches) max_matches = n_matches[i];
  max_match_index = i;
}
fMultiMatch = max_matches;
fMMFlash_Light      = 0;
fMMFlash_Time       = 0;
fMMFlash_TimeWidth  = 0;
fMMFlash_AbsTime    = 0;
fMMFlash_OnBeamTime = 0;
if (max_match_index != -1)
{
  auto const& opf = (*flashHandle)[max_match_index];

  fMMFlash_Light      = opf.TotalPE();
  fMMFlash_Time       = opf.Time();
  fMMFlash_TimeWidth  = opf.TimeWidth();
  fMMFlash_AbsTime    = opf.AbsTime();
  fMMFlash_OnBeamTime = opf.OnBeamTime();
}
fEventTree->Fill();



/* 
  for(size_t flash_index=0; flash_index < flashHandle->size(); ++flash_index) {

     art::Ptr<recob::OpFlash> flash(flashHandle,flash_index); 
      
      ::flashana::FlashMatch_t match; 
      std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
      std::cout<<"Size of flashmatch vector: "<<match_result_v.size()<<std::endl;
      std::cout<<"Size of opFlash vector: "     <<flashHandle->size()<<std::endl;
      for(size_t match_index=0; match_index < match_result_v.size(); ++match_index) {
     	 match = match_result_v.at(match_index);
	 std::cout<<"MatchID: "<<match.flash_id<< " trackID: "<<flash->ID()<<std::endl;
         if((int)match.tpc_id==flash->ID()) break; 
      }
      std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
      if (match.tpc_id==::flashana::kINVALID_ID) std::cout<<"INVALID TPC_ID"<<std::endl;  
      match = match_v.at(match_flash_index); 
      anab::FlashMatch Flash((double)match.score, (int)match.flash_id, (int)match.tpc_id, (bool)inbeam);
      flashmatchopflash->push_back(Flash);
      util::CreateAssn(*this, e, *flashmatchopflash, flash, *flashOpFlashAssociations,fSpillName);
   }
  
   for(auto const& match : match_result_v) {
	std::cout << "Match result ... "
	          << "Flash ID: " << match.flash_id
		  << " with "
		  << "TPC ID: " << match.tpc_id
		  << " ... Score: " << match.score
      		  << std::endl;
//      anab::FlashMatch Flash((double)match.score, (int)match.flash_id, (int)match.tpc_id, (bool)inbeam);
 //     flashmatch->push_back(Flash);
  }

  selfID.clear();
  isPrimary.clear();
  numDaughters.clear();
  parentID.clear();
  daughterIDs.clear();
  pdgCode.clear();
  isNeutrino.clear();
  isTrack.clear();
  trackID.clear();
  vertexID.clear();*/
   e.put(std::move(flashmatchtrack),fSpillName);   
   e.put(std::move(flashTrackAssociations),fSpillName);   
  // e.put(std::move(flashOpFlashAssociations),fSpillName);   
}

DEFINE_ART_MODULE(UBFlashMatching)
