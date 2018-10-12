////////////////////////////////////////////////
// LifetimeTrackByTrack module to be used for nearline monitoring
// Author: Sowjanya Gollapinni, first commit: July 11, 2016
// 
// The module selects tracks that cross anode and cathode
// and applies some angular cuts and produces dQ/ds vs drifttime
// scatter plot. The actual lifetime analysis is done in a 
// post processing script to extract QA/QC value.
// Associated fcl file is lifetime.fcl
///////////////////////////////////////////////

#ifndef LifetimeTrackByTrack_Module
#define LifetimeTrackByTrack_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"

const int nbin = 22; //split the total drift time into 22 bins
const double binsize = 100; //us
const int kMaxTrack  = 1000;  //maximum number of tracks
const int kNplanes   = 3;     //number of wire planes

using namespace std;

namespace {

// Local functions.

//========================================================================
// Length of reconstructed track, trajectory by trajectory.
double length(const recob::Track& track)
{
  double result = 0.;
  TVector3 disp = track.LocationAtPoint(0);
  int n = track.NumberTrajectoryPoints();

  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track.LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}

} // end namespace

//========================================================================

namespace microboone{
  
  class LifetimeTrackByTrack : public art::EDAnalyzer {
  public:
    
    explicit LifetimeTrackByTrack(fhicl::ParameterSet const& pset);
    virtual ~LifetimeTrackByTrack();
    
    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    
    void reconfigure(fhicl::ParameterSet const& pset);
    void reset();
    
  private:
    
    // the parameters we'll read from the .fcl
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
    
    TTree *fTree;
    
    // Event 
    Int_t    run;                  //run number
    Int_t    subrun;               //subrun number
    Int_t    event;                //event number
    Double_t evttime; 		   //event time in sec
    
    //Track information
    //Track plane data
    Short_t    _ntracks;               
    
    // more track info
    Float_t _trkstartx;     // starting x position.
    Float_t _trkstarty;     // starting y position.
    Float_t _trkstartz;     // starting z position.
    Float_t _trkendx;	   // ending x position.
    Float_t _trkendy;	   // ending y position.
    Float_t _trkendz;	   // ending z position.
    Float_t _trktheta;	   // theta.
    Float_t _trkphi;	   // phi.
    Float_t _trkstartdcosx;
    Float_t _trkstartdcosy;
    Float_t _trkstartdcosz;
    Float_t _trkenddcosx;
    Float_t _trkenddcosy;
    Float_t _trkenddcosz;
    Float_t _trkthetaxz;    // theta_xz.
    Float_t _trkthetayz;    // theta_xz.    
    Float_t _trkmom;	   // momentum.
    Float_t _trklen;	   // length.
    std::vector<float> _trk_x_v;
    std::vector<float> _trk_t_v;
    std::vector<float> _trk_y_v;
    std::vector<float> _trk_z_v;
    std::vector<float> _trk_dqds_v;
    
    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    art::ServiceHandle<geo::Geometry> geom;
    
  }; // class LifetimeTrackByTrack
  
  //========================================================================
  LifetimeTrackByTrack::LifetimeTrackByTrack(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }
  //========================================================================
  LifetimeTrackByTrack::~LifetimeTrackByTrack(){
    //destructor
  }
  //========================================================================
  void LifetimeTrackByTrack::reconfigure(fhicl::ParameterSet const& p){
    
    fTrackModuleLabel         = p.get<std::string>("TrackModuleLabel");
    fCalorimetryModuleLabel   = p.get<std::string>("CalorimetryModuleLabel");
    fSaveTrackInfo            = p.get< bool>("SaveTrackInfo",false);  
    
  }
  //========================================================================
  void LifetimeTrackByTrack::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    
    art::ServiceHandle<art::TFileService> tfs;
    
    fTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fTree->Branch("event", &event);
    fTree->Branch("run", &run);
    fTree->Branch("subrun", &subrun); 
    fTree->Branch("evttime", &evttime);
    fTree->Branch("_trklen",&_trklen,"trklen/F");
    fTree->Branch("_trkthetaxz", &_trkthetaxz, "trkthetaxz/F");
    fTree->Branch("_trkthetayz", &_trkthetayz, "trkthetayz/F");
    fTree->Branch("_trkstartx", &_trkstartx, "trkstartx/F");
    fTree->Branch("_trkstarty", &_trkstarty, "trkstarty/F");
    fTree->Branch("_trkstartz", &_trkstartz, "trkstartz/F");
    fTree->Branch("_trkendx", &_trkendx, "trkendx/F");
    fTree->Branch("_trkendy", &_trkendy, "trkendy/F");
    fTree->Branch("_trkendz", &_trkendz, "trkendz/F");
    fTree->Branch("_trktheta", &_trktheta, "trktheta/F");
    fTree->Branch("_trkphi", &_trkphi, "trkphi/F");
    fTree->Branch("_trkstartdcosx", &_trkstartdcosx,"trkstartdcosx/F");
    fTree->Branch("_trkstartdcosy", &_trkstartdcosy,"trkstartdcosy/F");
    fTree->Branch("_trkstartdcosz", &_trkstartdcosz,"trkstartdcosz/F");
    fTree->Branch("_trkenddcosx", &_trkenddcosx, "trkenddcosx/F");
    fTree->Branch("_trkenddcosy", &_trkenddcosy, "trkenddcosy/F");
    fTree->Branch("_trkenddcosz", &_trkenddcosz, "trkenddcosz/F");
    fTree->Branch("_trk_x_v","std::vector<float>",&_trk_x_v);
    fTree->Branch("_trk_t_v","std::vector<float>",&_trk_t_v);
    fTree->Branch("_trk_y_v","std::vector<float>",&_trk_y_v);
    fTree->Branch("_trk_z_v","std::vector<float>",&_trk_z_v);
    fTree->Branch("_trk_dqds_v","std::vector<float>",&_trk_dqds_v);
  }
  
  //========================================================================
  void LifetimeTrackByTrack::endJob(){     
    
  }
  
  //========================================================================
  void LifetimeTrackByTrack::beginRun(const art::Run& /*run*/){
    mf::LogInfo("LifetimeTrackByTrack")<<"begin run..."<<std::endl;
  }
  //========================================================================
  void LifetimeTrackByTrack::analyze( const art::Event& evt ){
    if (!evt.isRealData()) return;
    
    
    event  = evt.id().event(); 
    run    = evt.run();
    subrun = evt.subRun();
    
    art::Timestamp ts = evt.time();   
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();        
    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
       art::fill_ptr_vector(tracklist, trackListHandle);
       
    size_t NTracks = tracklist.size(); 
    _ntracks = (int) NTracks;
    
    for(unsigned int i=0; i<NTracks;++i){//loop over tracks
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;    
      
      reset();

      TVector3 pos, dir_start, dir_end, end;  
      
      double tlen = 0.,mom = 0.;
      int ntraj = 0;	
      ntraj = track.NumberTrajectoryPoints();

      if (ntraj > 0) {
        pos	 = track.Vertex();
        dir_start = track.VertexDirection();
        dir_end   = track.EndDirection();
        end	 = track.End();
        tlen	 = length(track);
	_trklen  = tlen;
	if(track.NumberTrajectoryPoints() > 0)
     	     mom = track.VertexMomentum();

	if (_trklen < 10.) continue; // skip tracks < 10 cm long

	double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
        double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());	
	
	//save track information of crossing tracks
	_trkthetaxz = theta_xz;
	_trkthetayz = theta_yz;
	_trkstartx        = pos.X();
	_trkstarty        = pos.Y();
	_trkstartz        = pos.Z();
	_trkendx	    = end.X();
	_trkendy	    = end.Y();
	_trkendz	    = end.Z();
	_trktheta	    = dir_start.Theta();
	_trkphi	    = dir_start.Phi();
	_trkstartdcosx    = dir_start.X();
	_trkstartdcosy    = dir_start.Y();
	_trkstartdcosz    = dir_start.Z();
	_trkenddcosx      = dir_end.X();
	_trkenddcosy      = dir_end.Y();
	_trkenddcosz      = dir_end.Z();
	_trkthetaxz       = theta_xz;
	_trkthetayz       = theta_yz;
	_trkmom	    = mom;
	_trklen	    = tlen;  

	art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
	if (fmcal.isValid()){
	  std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
	  if (calos.size() > 3) {
	    // if you get this message, there is probably a bug somewhere since
	    // the calorimetry planes should be 3.
	    mf::LogError("LifetimeTrackByTrack:limits")
	      << "the " << fTrackModuleLabel << " track #" << i
	      << " has " << calos.size() << " planes for calorimetry , only 3"
	      << " stored in tree";
	  }//if (calos.size() > 3)
	  for (size_t ical = 0; ical<calos.size(); ++ical){
	    if (!calos[ical]) continue;
	    if (!calos[ical]->PlaneID().isValid) continue;
	    int planenum = calos[ical]->PlaneID().Plane;
	    if (planenum!=2) continue; //planenum<0||planenum>2) continue;
	    const size_t NHits = calos[ical] -> dEdx().size();
	    if (planenum == 2){ //only interested in plane 2 for now
	      double minx = 1e10;
	      for(size_t iHit = 0; iHit < NHits; ++iHit) {	      
		const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
		if (TrkPos.X()<minx)
		  minx = TrkPos.X();
	      }// loop NHits
	      for(size_t iHit = 0; iHit < NHits; ++iHit) {				
		const auto& TrkPos1 = (calos[ical] -> XYZ())[iHit];
		float x = TrkPos1.X()-minx; //subtract the minx to get correct t0
		float y = TrkPos1.Y();
		float z = TrkPos1.Z();
		float t = x/(XDriftVelocity*1000); //change the velocity units to cm/ns to cm/us
		float dqds = calos[ical]->dQdx()[iHit];
		_trk_x_v.push_back(x);
		_trk_t_v.push_back(t);
		_trk_y_v.push_back(y);
		_trk_z_v.push_back(z);
		_trk_dqds_v.push_back(dqds);
	      }// loop NHits 
	    } // if planenum ==2	  
	  }// loop over ical	
	}// if fmcal.isValid()
      }// if ntraj>0

      fTree->Fill();

    }// loop over tracks 
    
}// end of analyze function




//========================================================================
void LifetimeTrackByTrack::reset(){
  
  _ntracks = 0;
  _trkstartx = -99999.;   
  _trkstarty = -99999.;   
  _trkstartz = -99999.;   
  _trkendx = -99999.;     
  _trkendy = -99999.;     
  _trkendz = -99999.;     
  _trktheta = -99999.;    
  _trkphi = -99999.; 
  _trkstartdcosx = -99999.;     
  _trkstartdcosy = -99999.;     
  _trkstartdcosz = -99999.;  
  _trkenddcosx = -99999.;
  _trkenddcosy = -99999.;
  _trkenddcosz = -99999.;
  _trkthetaxz = -99999.;   
  _trkthetaxz = -99999.;
  _trkthetayz = -99999.;        
  _trkthetayz = -99999.; 
  _trkmom = -99999.;      
  _trklen = -99999.;
  
}

//========================================================================
DEFINE_ART_MODULE(LifetimeTrackByTrack)

} 

#endif // LifetimeTrackByTrack_Module
