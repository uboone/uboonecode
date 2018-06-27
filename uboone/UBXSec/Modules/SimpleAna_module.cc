////////////////////////////////////////////////////////////////////////
// Class:       SimpleAna
// Plugin Type: analyzer (art v2_05_00)
// File:        SimpleAna_module.cc
//
// Generated at Tue Oct 31 16:36:23 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class SimpleAna
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module with simple example to retrieve results
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version analyzer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Tue Oct 31 16:36:23 2017
 *
 */
// Art include
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "lardata/DetectorInfo/DetectorClocks.h"
#include "uboone/EventWeight/MCEventWeight.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"


#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"


#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"




// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"
/*
#include "larsim/MCCheater/BackTrackerService.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
*/

#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/MCBase/MCHitCollection.h"

// Algorithms include
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/UBXSec/Algorithms/VertexCheck.h"
//#include "uboone/UBXSec/Algorithms/McPfpMatch.h"
#include "uboone/UBXSec/Algorithms/FindDeadRegions.h"
#include "uboone/UBXSec/Algorithms/MuonCandidateFinder.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"
#include "uboone/UBXSec/Algorithms/NuMuCCEventSelection.h"
#include "uboone/UBXSec/Algorithms/TrackQuality.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

//include header files for new backtracker test
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// Root include
#include "TString.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoVector.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <fstream>




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class SimpleAna;


class SimpleAna : public art::EDAnalyzer {
public:
  explicit SimpleAna(fhicl::ParameterSet const & p);
  virtual ~ SimpleAna();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleAna(SimpleAna const &) = delete;
  SimpleAna(SimpleAna &&) = delete;
  SimpleAna & operator = (SimpleAna const &) = delete;
  SimpleAna & operator = (SimpleAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  //add the new defined function here
  bool isContained(double x, double y, double z) const;
  int Topology(int nmuons, int nelectrons, int npions, int npi0, int nprotons);
  void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
 
  double GetEnergy(std::vector<art::Ptr<anab::Calorimetry>> calos);
  std::vector<double> GetdEdx(std::vector<art::Ptr<anab::Calorimetry>> calos);
  std::vector<double> GetRR(std::vector<art::Ptr<anab::Calorimetry>> calos);


private:
  std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
  //::ubana::McPfpMatch mcpfpMatcher;
  ::ubana::MuonCandidateFinder _muon_finder;
  ::ubana::FiducialVolume _fiducial_volume;




  std::string _tpcobject_producer;
  std::string _shower_producer;
  std::string _track_producer;
  std::string _particle_id_producer;
  std::string _mc_ghost_producer;
  std::string _pfp_producer;
  std::string _acpt_producer;
  std::string _calorimetry_producer;
  std::string _trigger_label;
  std::string _hit_producer;
  std::string fTrackMCSFitLabel;

  double fDistToEdgeX;
  double fDistToEdgeY;
  double fDistToEdgeZ;
  double fMinTrk2VtxDist;
  //Database to understand particle pdg
  const TDatabasePDG* _database_pdg = TDatabasePDG::Instance();

  // Detector info service
  ::detinfo::DetectorProperties const* fDetectorProperties;
  
  

  // Declare member data here.
  TTree* _tree;
  int _run, _subrun, _event;
  bool _status_ccpi0, _status_ccincl;
  double _bnb_correction;
  double _nu_energy, _nu_ccnc, _nu_mode, _nu_pdg;
  double _lep_costheta, _lep_phi, _lep_mom;

  int Nmuons_from_genie, Nprotons_from_genie, Nelectrons_from_genie;
  int Npionpms_from_genie, Npion0s_from_genie; 
  int Nothers_from_genie;
  int Nprotons_from_genie_200thresh, Nprotons_from_genie_300thresh, Nprotons_from_genie_400thresh;
  //variable for topology check
  float _fTruenuvrtxx, _fTruenuvrtxy, _fTruenuvrtxz;

  TLorentzVector *fHitNucP4;

  std::vector<double> *trueProtonsTrueMomentum;
  std::vector<double> *trueProtonsTrueTheta;
  std::vector<double> *trueProtonsTruePhi;
  std::vector<double> *trueProtonsEndMomentum;
  std::vector<std::string> *trueProtonsEndProcess;

 
  
  int truthtop;
  int truthtop_200thresh;
  int truthtop_300thresh;
  int truthtop_400thresh;

  int truthtop_genie;
  int truthtop_genie_200thresh;
  int truthtop_genie_300thresh;
  int truthtop_genie_400thresh;


  std::vector<int> *fhg4parpdg;
  std::vector<int> *fhg4parstatus;
  std::vector<float> *fhg4parpx;
  std::vector<float> *fhg4parpy;
  std::vector<float> *fhg4parpz;
  std::vector<float> *fhg4partheta;
  std::vector<float> *fhg4parphi;
  std::vector<float> *fhg4parp;




  std::vector<string>   *_fg4processname; //Physics process by which the particle was created
  std::vector<int>   *_fg4TrackId;
  std::vector<int>   *_fg4Mother;


  std::vector<int> *_fg4parpdg; 
  std::vector<int> *_fg4parstatus; //status code of secondary particles assigned by geant 
  std::vector<float> *_fg4EngS; //energy of particles at started in GeV; 
  std::vector<float> *_fg4EngE; //energy of particles at ended in GeV; 
  std::vector<float> *_fg4px; 
  std::vector<float> *_fg4py; 
  std::vector<float> *_fg4pz; 
  std::vector<float> *_fg4P; 
  std::vector<float> *_fg4Mass; 

  std::vector<float> *_fg4StartPointx; 
  std::vector<float> *_fg4StartPointy; 
  std::vector<float> *_fg4StartPointz; 
  std::vector<float> *_fg4EndPointx; 
  std::vector<float> *_fg4EndPointy; 
  std::vector<float> *_fg4EndPointz; 


  std::vector<float> *_fg4StartT;
  std::vector<float> *_fg4EndT;
  std::vector<float> *_fg4theta;
  std::vector<float> *_fg4phi;
  std::vector<float> *_fg4theta_xz;
  std::vector<float> *_fg4theta_yz;

  std::vector<simb::Origin_t>   *_fg4origin;
  std::vector<int> *_fg4MCTruthIndex;
  std::vector<int>   *_fg4NumberDaughters;
  //============================================================ 
  int _npfps, _npfp_tracks, _npfp_showers;
  bool _pfp_ntrk2flag;
  bool _pfp_pdqdxflag, _pfp_upinFVflag;



  std::vector<bool> *pfp_isShower;
  std::vector<bool> *pfp_isTrack;
  std::vector<bool> *pfp_isPrimary;
  std::vector<int> *pfp_ndaughters;

  std::vector<int> *pfp_pdg;
  std::vector<int> *pfp_origin;
  std::vector<int> *pfp_status;
  std::vector<int> *pfp_parId;
  std::vector<float> *pfp_theta;
  std::vector<float> *pfp_costheta;
  std::vector<float> *pfp_phi;
  std::vector<float> *pfp_mom;
  std::vector<float> *pfp_startx;
  std::vector<float> *pfp_starty;
  std::vector<float> *pfp_startz;
  std::vector<float> *pfp_endx;
  std::vector<float> *pfp_endy;
  std::vector<float> *pfp_endz;
  std::vector<float> *pfp_endE;
  std::vector<float> *pfp_KE;
  std::vector<float> *pfp_Mass;

 
  std::vector<bool> *track_pfp_upflag;
  std::vector<int> *track_pfp_Id;
  std::vector<float> *track_pfp_length;
  std::vector<float> *track_pfp_theta;
  std::vector<float> *track_pfp_costheta;
  std::vector<float> *track_pfp_phi;
  std::vector<float> *track_pfp_startx;
  std::vector<float> *track_pfp_starty;
  std::vector<float> *track_pfp_startz;
  std::vector<float> *track_pfp_endx;
  std::vector<float> *track_pfp_endy;
  std::vector<float> *track_pfp_endz;
  std::vector<float> *track_pfp_Mom; 
  std::vector<float> *track_pfp_trunmeandqdx;
  std::vector<float> *track_pfp_trunmeandqdx_U;
  std::vector<float> *track_pfp_trunmeandqdx_V;

  std::vector<float> *track_pfp_pida;
  std::vector<int> *track_pfp_nhits;
  std::vector<std::vector<double>> *track_pfp_dEdx;
  std::vector<std::vector<double>> *track_pfp_RR;

  float track_pfp_mom_mucand; 
  std::vector<float> *track_pfp_mom_pcand;

 //===============================================================
  //----add the variables for CC1uNP analysis here-----------




  int _ntrks;
  bool _ccinclflag=false;
  bool _ntrk2flag=false;
  bool _noextrkflag=false;
  bool _upinFVflag=false;
  bool _pdqdxflag=false;

  bool OOFVflag=false;
 
  
  //vertex info
  float _nuvtxx_reco, _nuvtxy_reco, _nuvtxz_reco;

  //tracks reco info
  std::vector<int> *upflag;
  std::vector<int> *trackId;
  std::vector<float> *tracklength;
  std::vector<float> *trackmom;
  std::vector<float> *trackcostheta;
  std::vector<float> *trackphi;

  std::vector<double> *trackstartx;
  std::vector<double> *trackstarty;
  std::vector<double> *trackstartz;

  std::vector<double> *trackendx;
  std::vector<double> *trackendy;
  std::vector<double> *trackendz;

  std::vector<float>  *tracktrunmeandqdx;
  std::vector<float>  *tracktrunmeandqdx_U;
  std::vector<float>  *tracktrunmeandqdx_V;
  std::vector<std::vector<double>>  *trackdEdx;
  std::vector<std::vector<double>> *trackRR;
  std::vector<int> *tracknhits;

  //tracks sim info
  std::vector<int> *trackcand_origin;
  std::vector<int> *trackcand_nuset;
  std::vector<int> *trackcand_parPDG;
  std::vector<int> *trackcand_parStatusCode;
  std::vector<float> *trackcand_parTheta;                  
  std::vector<float> *trackcand_parCosTheta;
  std::vector<float> *trackcand_parSinTheta;                  
  std::vector<float> *trackcand_parE;        
  std::vector<float> *trackcand_parMass;
  std::vector<float> *trackcand_parKE;
  std::vector<float> *trackcand_parEndE;
  std::vector<float> *trackcand_parPx;
  std::vector<float> *trackcand_parPy;
  std::vector<float> *trackcand_parPz;
  std::vector<float> *trackcand_parPhi;
  std::vector<float> *trackcand_parCosPhi;
  std::vector<float> *trackcand_parSinPhi;

  float trackmom_mucand;
  std::vector<float> *trackmom_pcand;

  TTree* _cc1unptree;

  //-------------------------------------------------------------
};


SimpleAna::SimpleAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  //set the fhicl parameters here=================================================================
  _tpcobject_producer        = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::UBXSec");
  _shower_producer           = p.get<std::string>("ShowerProducer",     "pandoraNu::UBXSec");
  _track_producer            = p.get<std::string>("TrackProducer",      "pandoraNu::UBXSec");
  _particle_id_producer           = p.get<std::string>("ParticleIDProducer", "pandoraNupid::UBXSec");   
  _mc_ghost_producer           = p.get<std::string>("MCGhostProducer", "RecoTrueMatcher");

  _pfp_producer                   = p.get<std::string>("PFParticleProducer", "pandoraNu::UBXSec");

  _acpt_producer                  = p.get<std::string>("ACPTProducer", "pandoraCosmicT0Reco");
  _calorimetry_producer           = p.get<std::string>("CalorimetryProducer", "pandoraNucalo");
  _hit_producer                    = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval::UBXSec");
  //_hit_producer                    = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval");

  //_calorimetry_producer           = p.get<std::string>("CalorimetryProducer");

  _trigger_label = p.get<std::string>("TriggerProduct", "triggersim");

  fMinTrk2VtxDist          = p.get      ("MinTrk2VtxDist", 5.);
  fTrackMCSFitLabel        = p.get< std::string > ("TrackMCSFitLabel", "pandoraNuMCSMu" );

    auto const* geom = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fDistToEdgeX             = geom->DetHalfWidth()   - p.get("DistToEdgeX",   10.);    
    fDistToEdgeY             = geom->DetHalfHeight()  - p.get("DistToEdgeY",   20.);    
    fDistToEdgeZ             = geom->DetLength() / 2. - p.get("DistToEdgeZ",   10.);        
 
  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geom->DetHalfHeight(),
                             2.*geom->DetHalfWidth(),
                             geom->DetLength());

  _fiducial_volume.PrintConfig();



  _muon_finder.Configure(p.get<fhicl::ParameterSet>("MuonCandidateFinderSettings"));

  _muon_finder.PrintConfig();

 //Get the tool for MC truth Matching
  const fhicl::ParameterSet& truthParams = p.get<fhicl::ParameterSet>("MCTruthMatching");

  if (truthParams.get<std::string>("tool_type") == "AssociationsTruth")
  {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
  }
  else
  {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::BackTrackerTruth(truthParams));
  }
  



  //===============================================================================================
  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree","");
  _tree->Branch("run",             &_run,              "run/I");
  _tree->Branch("subrun",          &_subrun,           "subrun/I");
  _tree->Branch("event",           &_event,            "event/I");
  _tree->Branch("bnb_correction",  &_bnb_correction,   "bnb_correction/D");
  _tree->Branch("status_ccincl",   &_status_ccincl,    "status_ccincl/O");
  _tree->Branch("status_ccpi0",    &_status_ccpi0,     "status_ccpi0/O");
  _tree->Branch("nu_energy",       &_nu_energy,        "nu_energy/D");
  //declear all the variables for CC1uNP here------------------------
  //art::ServiceHandle<art::TFileService> fs;
  _cc1unptree = fs->make<TTree>("cc1unptree","");
  _cc1unptree->Branch("run",             &_run,              "run/I");
  _cc1unptree->Branch("subrun",          &_subrun,           "subrun/I");
  _cc1unptree->Branch("event",           &_event,            "event/I");
  _cc1unptree->Branch("bnb_correction",  &_bnb_correction,   "bnb_correction/D");
  _cc1unptree->Branch("status_ccincl",   &_status_ccincl,    "status_ccincl/O");
  _cc1unptree->Branch("status_ccpi0",    &_status_ccpi0,     "status_ccpi0/O");
  _cc1unptree->Branch("nu_energy",       &_nu_energy,        "nu_energy/D");
  _cc1unptree->Branch("nu_mode",         &_nu_mode,          "nu_mode/D");
  _cc1unptree->Branch("_lep_mom",        &_lep_mom,          "_lep_mom/D");
  _cc1unptree->Branch("_lep_costheta",        &_lep_costheta,          "_lep_costheta/D");
  _cc1unptree->Branch("_lep_phi",        &_lep_phi,          "_lep_phi/D");
  //_cc1unptree->Branch("fHitNucP4","TLorentzVector",&fHitNucP4);

  _cc1unptree->Branch("_fTruenuvrtxx",   &_fTruenuvrtxx,     "_fTruenuvrtxx/F");
  _cc1unptree->Branch("_fTruenuvrtxy",   &_fTruenuvrtxy,     "_fTruenuvrtxy/F");
  _cc1unptree->Branch("_fTruenuvrtxz",   &_fTruenuvrtxz,     "_fTruenuvrtxz/F");




  _cc1unptree->Branch("trueProtonsTrueMomentum","std::vector<double>",&trueProtonsTrueMomentum);
  _cc1unptree->Branch("trueProtonsTrueTheta","std::vector<double>",&trueProtonsTrueTheta);
  _cc1unptree->Branch("trueProtonsTruePhi","std::vector<double>",&trueProtonsTruePhi);
  _cc1unptree->Branch("trueProtonsEndMomentum","std::vector<double>",&trueProtonsEndMomentum);
  _cc1unptree->Branch("trueProtonsEndProcess","std::vector<std::string>",&trueProtonsEndProcess);
  // variables at G4 stage
  _cc1unptree->Branch("fhg4parpdg","std::vector<int>", &fhg4parpdg);
  _cc1unptree->Branch("fhg4parstatus",  "std::vector<int>",  &fhg4parstatus);
  _cc1unptree->Branch("fhg4parpx",  "std::vector<float>",  &fhg4parpx);
  _cc1unptree->Branch("fhg4parpy",  "std::vector<float>",  &fhg4parpy);
  _cc1unptree->Branch("fhg4parpz",  "std::vector<float>",  &fhg4parpz);
  _cc1unptree->Branch("fhg4partheta",  "std::vector<float>",  &fhg4partheta);
  _cc1unptree->Branch("fhg4parphi",  "std::vector<float>",  &fhg4parphi);
  _cc1unptree->Branch("fhg4parp",  "std::vector<float>",  &fhg4parp);

  _cc1unptree->Branch("OOFVflag", &OOFVflag,   "OOFVflag/O"); 
  //all the labels for cuts  
  _cc1unptree->Branch("ccinclflag",      &_ccinclflag,       "ccinclflag/O");
  _cc1unptree->Branch("ntrk2flag",       &_ntrk2flag,        "ntrk2flag/O");
  _cc1unptree->Branch("noextrkflag",     &_noextrkflag,      "noextrkflag/O");
  _cc1unptree->Branch("upinFVflag",      &_upinFVflag,       "upinFVflag/O"); 
  _cc1unptree->Branch("pdqdxflag",       &_pdqdxflag,        "pdqdxflag/O");   
  //g4 stage variables
  //-------------------------------------------------------------------------------

  _cc1unptree->Branch("pfp_ntrk2flag",   &_pfp_ntrk2flag, "pfp_ntrk2flag/O");
  _cc1unptree->Branch("pfp_upinFVflag", &_pfp_upinFVflag,     "pfp_upinFVflag/O");
  _cc1unptree->Branch("pfp_pdqdxflag", &_pfp_pdqdxflag,       "pfp_pdqdxflag/O");


  _cc1unptree->Branch("npfps",      &_npfps,         "npfps/I");
  _cc1unptree->Branch("npfp_tracks",      &_npfp_tracks,         "npfp_tracks/I");
  _cc1unptree->Branch("npfp_showers",      &_npfp_showers,         "npfp_showers/I");


  _cc1unptree->Branch("pfp_isTrack",    "std::vector<bool>", &pfp_isTrack);
  _cc1unptree->Branch("pfp_isShower",    "std::vector<bool>", &pfp_isShower);
  _cc1unptree->Branch("pfp_isPrimary",   "std::vector<bool>", &pfp_isPrimary);
  _cc1unptree->Branch("pfp_ndaughters",  "std::vector<int>", &pfp_ndaughters);

  _cc1unptree->Branch("pfp_pdg",    "std::vector<int>", &pfp_pdg);
  _cc1unptree->Branch("pfp_origin",    "std::vector<int>", &pfp_origin);
  _cc1unptree->Branch("pfp_parId",    "std::vector<int>", &pfp_parId);
  _cc1unptree->Branch("pfp_status",    "std::vector<int>", &pfp_status);
  _cc1unptree->Branch("pfp_theta",    "std::vector<float>", &pfp_theta);
  _cc1unptree->Branch("pfp_costheta",    "std::vector<float>", &pfp_costheta);
  _cc1unptree->Branch("pfp_phi",    "std::vector<float>", &pfp_phi);
  _cc1unptree->Branch("pfp_mom",    "std::vector<float>", &pfp_mom);
  _cc1unptree->Branch("pfp_startx",    "std::vector<float>", &pfp_startx);
  _cc1unptree->Branch("pfp_starty",    "std::vector<float>", &pfp_starty);
  _cc1unptree->Branch("pfp_startz",    "std::vector<float>", &pfp_startz);
  _cc1unptree->Branch("pfp_endx",    "std::vector<float>", &pfp_endx);
  _cc1unptree->Branch("pfp_endy",    "std::vector<float>", &pfp_endy);
  _cc1unptree->Branch("pfp_endz",    "std::vector<float>", &pfp_endz);
  _cc1unptree->Branch("pfp_endE", "std::vector<float>", &pfp_endE); 
  _cc1unptree->Branch("pfp_KE", "std::vector<float>", &pfp_KE); 
  _cc1unptree->Branch("pfp_Mass", "std::vector<float>", &pfp_Mass); 

  _cc1unptree->Branch("track_pfp_upflag", "std::vector<bool>", &track_pfp_upflag);
  _cc1unptree->Branch("track_pfp_Id", "std::vector<int>", &track_pfp_Id); 
  _cc1unptree->Branch("track_pfp_length", "std::vector<float>", &track_pfp_length); 
  _cc1unptree->Branch("track_pfp_theta", "std::vector<float>", &track_pfp_theta); 
  _cc1unptree->Branch("track_pfp_costheta", "std::vector<float>", &track_pfp_costheta); 
  _cc1unptree->Branch("track_pfp_phi", "std::vector<float>", &track_pfp_phi); 
  _cc1unptree->Branch("track_pfp_startx", "std::vector<float>", &track_pfp_startx); 
  _cc1unptree->Branch("track_pfp_starty", "std::vector<float>", &track_pfp_starty); 
  _cc1unptree->Branch("track_pfp_startz", "std::vector<float>", &track_pfp_startz); 
  _cc1unptree->Branch("track_pfp_endx", "std::vector<float>", &track_pfp_endx); 
  _cc1unptree->Branch("track_pfp_endy", "std::vector<float>", &track_pfp_endy); 
  _cc1unptree->Branch("track_pfp_endz", "std::vector<float>", &track_pfp_endz); 
  _cc1unptree->Branch("track_pfp_Mom", "std::vector<float>", &track_pfp_Mom); 
  _cc1unptree->Branch("track_pfp_trunmeandqdx", "std::vector<float>", &track_pfp_trunmeandqdx); 
  _cc1unptree->Branch("track_pfp_trunmeandqdx_U", "std::vector<float>", &track_pfp_trunmeandqdx_U); 
  _cc1unptree->Branch("track_pfp_trunmeandqdx_V", "std::vector<float>", &track_pfp_trunmeandqdx_V); 

  _cc1unptree->Branch("track_pfp_pida", "std::vector<float>", &track_pfp_pida);
  _cc1unptree->Branch("track_pfp_nhits", "std::vector<int>", &track_pfp_nhits);
  _cc1unptree->Branch("track_pfp_dEdx", "std::vector<std::vector<double>>", &track_pfp_dEdx);
  _cc1unptree->Branch("track_pfp_RR", "std::vector<std::vector<double>>", &track_pfp_RR);

 
  _cc1unptree->Branch("track_pfp_mom_mucand", &track_pfp_mom_mucand, "track_pfp_mom_mucand/F"); 
  _cc1unptree->Branch("track_pfp_mom_pcand", "std::vector<float>", &track_pfp_mom_pcand);


  //all the kinematic variables for cc1unp
  //--------------------------------------------------------------------------
  _cc1unptree->Branch("ntrks",           &_ntrks,             "ntrks/I");

  _cc1unptree->Branch("_nuvtxx_reco",    &_nuvtxx_reco,        "_nuvtxx_reco/F");
  _cc1unptree->Branch("_nuvtxy_reco",    &_nuvtxy_reco,        "_nuvtxy_reco/F");
  _cc1unptree->Branch("_nuvtxz_reco",    &_nuvtxz_reco,        "_nuvtxz_reco/F");

  _cc1unptree->Branch("upflag",   "std::vector<int>", &upflag);

  _cc1unptree->Branch("trackId",   "std::vector<int>", &trackId);
  _cc1unptree->Branch("tracklength","std::vector<float>",&tracklength);
  _cc1unptree->Branch("trackmom",   "std::vector<float>",&trackmom);
  _cc1unptree->Branch("trackcostheta", "std::vector<float>", &trackcostheta);
  _cc1unptree->Branch("trackphi",   "std::vector<float>", &trackphi);

  _cc1unptree->Branch("trackstartx",   "std::vector<double>", &trackstartx);
  _cc1unptree->Branch("trackstarty",   "std::vector<double>", &trackstarty);
  _cc1unptree->Branch("trackstartz",   "std::vector<double>", &trackstartz);

  _cc1unptree->Branch("trackendx",   "std::vector<double>", &trackendx);
  _cc1unptree->Branch("trackendy",   "std::vector<double>", &trackendy);
  _cc1unptree->Branch("trackendz",   "std::vector<double>", &trackendz);
   
  _cc1unptree->Branch("tracktrunmeandqdx",   "std::vector<float>", &tracktrunmeandqdx);
  _cc1unptree->Branch("tracktrunmeandqdx_U",   "std::vector<float>", &tracktrunmeandqdx_U);
  _cc1unptree->Branch("tracktrunmeandqdx_V",   "std::vector<float>", &tracktrunmeandqdx_V);
  _cc1unptree->Branch("trackdEdx","std::vector<std::vector<double>>", &trackdEdx);
  _cc1unptree->Branch("trackRR", "std::vector<std::vector<double>>", &trackRR);
  _cc1unptree->Branch("tracknhits", "std::vector<int>", &tracknhits);
  _cc1unptree->Branch("trackcand_origin", "std::vector<int>", &trackcand_origin);
  _cc1unptree->Branch("trackcand_nuset", "std::vector<int>", &trackcand_nuset);
  _cc1unptree->Branch("trackcand_parPDG", "std::vector<int>", &trackcand_parPDG);
  _cc1unptree->Branch("trackcand_parStatusCode", "std::vector<int>", &trackcand_parStatusCode);
  _cc1unptree->Branch("trackcand_parTheta", "std::vector<float>", &trackcand_parTheta);
  _cc1unptree->Branch("trackcand_parCosTheta", "std::vector<float>", &trackcand_parCosTheta);
  _cc1unptree->Branch("trackcand_parSinTheta", "std::vector<float>", &trackcand_parSinTheta);
  _cc1unptree->Branch("trackcand_parE", "std::vector<float>", &trackcand_parE);
  _cc1unptree->Branch("trackcand_parMass", "std::vector<float>", &trackcand_parMass);
  _cc1unptree->Branch("trackcand_parKE", "std::vector<float>", &trackcand_parKE);
  _cc1unptree->Branch("trackcand_parEndE", "std::vector<float>", &trackcand_parEndE);
  _cc1unptree->Branch("trackcand_parPx", "std::vector<float>", &trackcand_parPx);
  _cc1unptree->Branch("trackcand_parPy", "std::vector<float>", &trackcand_parPy);
  _cc1unptree->Branch("trackcand_parPz", "std::vector<float>", &trackcand_parPz);
  _cc1unptree->Branch("trackcand_parCosPhi", "std::vector<float>", &trackcand_parCosPhi);
  _cc1unptree->Branch("trackcand_parSinPhi", "std::vector<float>", &trackcand_parSinPhi);


  _cc1unptree->Branch("trackmom_mucand", &trackmom_mucand, "trackmom_mucand/F"); 
  _cc1unptree->Branch("trackmom_pcand", "std::vector<float>", &trackmom_pcand); 


  _cc1unptree->Branch("truthtop",         &truthtop,            "truthtop/I");
  _cc1unptree->Branch("truthtop_200thresh",      &truthtop_200thresh,            "truthtop_200thresh/I");
  _cc1unptree->Branch("truthtop_300thresh",      &truthtop_300thresh,            "truthtop_300thresh/I");
  _cc1unptree->Branch("truthtop_400thresh",      &truthtop_400thresh,            "truthtop_400thresh/I");

  _cc1unptree->Branch("truthtop_genie",         &truthtop_genie,            "truthtop_genie/I");
  _cc1unptree->Branch("truthtop_genie_200thresh",      &truthtop_genie_200thresh,            "truthtop_genie_200thresh/I");
  _cc1unptree->Branch("truthtop_genie_300thresh",      &truthtop_genie_300thresh,            "truthtop_genie_300thresh/I");
  _cc1unptree->Branch("truthtop_genie_400thresh",      &truthtop_genie_400thresh,            "truthtop_genie_400thresh/I");
   //
  // Save to tree
  //

  //====================================================
  trueProtonsTrueMomentum=new std::vector<double>;
  trueProtonsTrueTheta=new std::vector<double>;;
  trueProtonsTruePhi=new std::vector<double>;;
  trueProtonsEndMomentum=new std::vector<double>;;
  trueProtonsEndProcess=new std::vector<std::string>;;






  //---------------------------------------------------
  fhg4parpdg=new std::vector<int>;
  fhg4parstatus=new std::vector<int>;
  fhg4parpx=new std::vector<float>;
  fhg4parpy=new std::vector<float>;
  fhg4parpz=new std::vector<float>;
  fhg4partheta=new std::vector<float>;
  fhg4parphi=new std::vector<float>;
  fhg4parp=new std::vector<float>;

  _fg4processname=new std::vector<string>; //Physics process by which the particle was created
  _fg4TrackId=new std::vector<int>;
  _fg4Mother=new std::vector<int>;


  _fg4parpdg=new std::vector<int>; 
  _fg4parstatus=new std::vector<int>; //status code of secondary particles assigned by geant 
  _fg4EngS=new std::vector<float>; //energy of particles at started in GeV; 
  _fg4EngE=new std::vector<float>; //energy of particles at ended in GeV; 
  _fg4px=new std::vector<float>; 
  _fg4py=new std::vector<float>; 
  _fg4pz=new std::vector<float>; 
  _fg4P=new std::vector<float>; 
  _fg4Mass=new std::vector<float>; 

  _fg4StartPointx=new std::vector<float>; 
  _fg4StartPointy=new std::vector<float>; 
  _fg4StartPointz=new std::vector<float>; 
  _fg4EndPointx=new std::vector<float>; 
  _fg4EndPointy=new std::vector<float>; 
  _fg4EndPointz=new std::vector<float>; 


  _fg4StartT=new std::vector<float>;
  _fg4EndT=new std::vector<float>;
  _fg4theta=new std::vector<float>;
  _fg4phi=new std::vector<float>;
  _fg4theta_xz=new std::vector<float>;
  _fg4theta_yz=new std::vector<float>;

  _fg4origin=new std::vector<simb::Origin_t>;
  _fg4MCTruthIndex=new std::vector<int>;
  _fg4NumberDaughters=new std::vector<int>;

  pfp_isTrack=new std::vector<bool>;
  pfp_isShower=new std::vector<bool>;
  pfp_isPrimary=new std::vector<bool>;
  pfp_ndaughters=new std::vector<int>;

  pfp_pdg=new std::vector<int>;
  pfp_origin=new std::vector<int>;
  pfp_status=new std::vector<int>;
  pfp_parId=new std::vector<int>;
  pfp_theta=new std::vector<float>;
  pfp_costheta=new std::vector<float>;
  pfp_phi=new std::vector<float>;
  pfp_mom=new std::vector<float>;
  pfp_startx=new std::vector<float>;
  pfp_starty=new std::vector<float>;
  pfp_startz=new std::vector<float>;
  pfp_endx=new std::vector<float>;
  pfp_endy=new std::vector<float>;
  pfp_endz=new std::vector<float>;
  pfp_endE=new std::vector<float>;
  pfp_KE=new std::vector<float>;
  pfp_Mass=new std::vector<float>;  


  track_pfp_upflag=new std::vector<bool>;
  track_pfp_Id=new std::vector<int>;
  track_pfp_length=new std::vector<float>;
  track_pfp_theta=new std::vector<float>;
  track_pfp_costheta=new std::vector<float>;
  track_pfp_phi=new std::vector<float>;
  track_pfp_startx=new std::vector<float>;
  track_pfp_starty=new std::vector<float>;
  track_pfp_startz=new std::vector<float>;
  track_pfp_endx=new std::vector<float>;
  track_pfp_endy=new std::vector<float>;
  track_pfp_endz=new std::vector<float>;
  track_pfp_Mom=new std::vector<float>;
  track_pfp_trunmeandqdx=new std::vector<float>;
  track_pfp_trunmeandqdx_U=new std::vector<float>;
  track_pfp_trunmeandqdx_V=new std::vector<float>;
 
  track_pfp_pida=new std::vector<float>; 
  track_pfp_nhits=new std::vector<int>;
  track_pfp_dEdx=new std::vector<std::vector<double>>;
  track_pfp_RR=new std::vector<std::vector<double>>;

  track_pfp_mom_pcand=new std::vector<float>;


 
  upflag=new std::vector<int>;
  trackId=new std::vector<int>;
  tracklength = new std::vector<float>;
  trackcostheta = new std::vector<float>;
  trackphi = new std::vector<float>;
  trackmom = new std::vector<float>;

  trackstartx = new std::vector<double>;
  trackstarty = new std::vector<double>;
  trackstartz = new std::vector<double>;
  trackendx = new std::vector<double>;
  trackendy = new std::vector<double>;
  trackendz = new std::vector<double>;

  tracktrunmeandqdx=new std::vector<float>;
  tracktrunmeandqdx_U=new std::vector<float>;
  tracktrunmeandqdx_V=new std::vector<float>;

  trackdEdx=new std::vector<std::vector<double>>;
  trackRR=new std::vector<std::vector<double>>;

  tracknhits=new std::vector<int>;

  trackmom_pcand=new std::vector<float>;

  trackcand_origin=new std::vector<int>;
  trackcand_nuset=new std::vector<int>;

  trackcand_parPDG=new std::vector<int>;
  trackcand_parStatusCode=new std::vector<int>;
  trackcand_parTheta=new std::vector<float>;                  
  trackcand_parCosTheta=new std::vector<float>;
  trackcand_parSinTheta=new std::vector<float>;                  
  trackcand_parE=new std::vector<float>;        
  trackcand_parMass=new std::vector<float>;
  trackcand_parKE=new std::vector<float>;
  trackcand_parEndE=new std::vector<float>;
  trackcand_parPx=new std::vector<float>;
  trackcand_parPy=new std::vector<float>;
  trackcand_parPz=new std::vector<float>;
  trackcand_parPhi=new std::vector<float>;
  trackcand_parCosPhi=new std::vector<float>;
  trackcand_parSinPhi=new std::vector<float>;

  fHitNucP4 = new TLorentzVector(-999,-999,-999,-999);

 
  //-----------------------------------------------------------------
}


// Destructor
 SimpleAna::~ SimpleAna()
{
  delete fHitNucP4;
  delete trueProtonsTrueMomentum;
  delete trueProtonsTrueTheta;
  delete trueProtonsTruePhi;
  delete trueProtonsEndMomentum;
  delete trueProtonsEndProcess;


  delete fhg4parpdg;
  delete fhg4parstatus;
  delete fhg4parpx;
  delete fhg4parpy;
  delete fhg4parpz;
  delete fhg4partheta;
  delete fhg4parphi;
  delete fhg4parp;

  delete _fg4processname; 
  delete _fg4TrackId;
  delete _fg4Mother;


  delete _fg4parpdg; 
  delete _fg4parstatus; 
  delete _fg4EngS;  
  delete _fg4EngE; 
  delete _fg4px; 
  delete _fg4py; 
  delete _fg4pz; 
  delete _fg4P; 
  delete _fg4Mass; 

  delete _fg4StartPointx; 
  delete _fg4StartPointy; 
  delete _fg4StartPointz; 
  delete _fg4EndPointx; 
  delete _fg4EndPointy; 
  delete _fg4EndPointz; 


  delete _fg4StartT;
  delete _fg4EndT;
  delete _fg4theta;
  delete _fg4phi;
  delete _fg4theta_xz;
  delete _fg4theta_yz;

  delete _fg4origin;
  delete _fg4MCTruthIndex;
  delete _fg4NumberDaughters;
  //==================================================
  delete pfp_isTrack;
  delete pfp_isShower; 
  delete pfp_isPrimary;
  delete pfp_ndaughters;

  delete pfp_pdg;
  delete pfp_origin;
  delete pfp_status;
  delete pfp_parId;
  delete pfp_theta;
  delete pfp_costheta;
  delete pfp_phi;
  delete pfp_mom;
  delete pfp_startx;
  delete pfp_starty;
  delete pfp_startz;
  delete pfp_endx;
  delete pfp_endy;
  delete pfp_endz;
  delete pfp_endE;
  delete pfp_KE;
  delete pfp_Mass;

  delete track_pfp_upflag;
  delete track_pfp_Id;
  delete track_pfp_length;
  delete track_pfp_theta;
  delete track_pfp_costheta;
  delete track_pfp_phi;
  delete track_pfp_startx;
  delete track_pfp_starty;
  delete track_pfp_startz;
  delete track_pfp_endx;
  delete track_pfp_endy;
  delete track_pfp_endz;
  delete track_pfp_Mom;
  delete track_pfp_trunmeandqdx;
  delete track_pfp_trunmeandqdx_U;
  delete track_pfp_trunmeandqdx_V;
 
  delete track_pfp_pida; 
  delete track_pfp_nhits;
  delete track_pfp_dEdx;
  delete track_pfp_RR;

  delete track_pfp_mom_pcand;



  //============================================
  delete upflag;
  delete trackId;
  delete tracklength;
  delete trackcostheta;
  delete trackphi;
  delete trackmom;

  delete trackstartx;
  delete trackstarty;
  delete trackstartz;
  delete trackendx;
  delete trackendy;
  delete trackendz;
  
  delete tracktrunmeandqdx; 
  delete tracktrunmeandqdx_U; 
  delete tracktrunmeandqdx_V; 

  delete trackdEdx;
  delete trackRR;
  delete tracknhits;

  delete trackmom_pcand;

  delete trackcand_origin;
  delete trackcand_nuset;
  delete trackcand_parPDG;
  delete trackcand_parStatusCode;
  delete trackcand_parTheta;                  
  delete trackcand_parCosTheta;
  delete trackcand_parSinTheta;                  
  delete trackcand_parE;        
  delete trackcand_parMass;
  delete trackcand_parKE;
  delete trackcand_parEndE;
  delete trackcand_parPx;
  delete trackcand_parPy;
  delete trackcand_parPz;
  delete trackcand_parPhi;
  delete trackcand_parCosPhi;
  delete trackcand_parSinPhi;

  
} 
bool SimpleAna::isContained(double x, double y, double z) const
{
    auto const* geom = lar::providerFrom<geo::Geometry>(); // geometry is needed to go from OpChannel to OpDet 


    double distInX = x - geom->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * geom->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}




int SimpleAna::Topology(int nmuons, int nelectrons, int npions, int npi0, int nprotons)
{
  //// This function return the true topology of the event, numu & anti-numu                                                        
  ////1. CC0Pion0Proton                                                                                                             
  ////2. CC0Pion1Proton                                  
  ////3. CC0Pion2Proton                                                                                                             
  ////4. CC0PionNProton                                                                                                             
  ////5. CC1PionNProton (1 Pion= 1 charged pion || 1 neutral pion)
  ////6. CCNPionNProton    
  ////7. CCnue-antinue   
  ////8. NC                                                                                            
  ////9. OOFV (nu event out of FV)                                                                                                  
  ////10. Cosmic                                                                                                       
  ////11. Other (just in case, let's check!)                                                                        
  /// 12. 2&3&4 -> CC0PinProton (N>0) /// *** add this one!                                                                         
  /// e.g. numu CC inclusive= Topology >0 && Topology < 7 
  
  
  if (nmuons >0 && (nelectrons + npions + npi0 ) == 0 && nprotons ==0 ) return 1;
  if (nmuons >0 && (nelectrons + npions + npi0 ) == 0 && nprotons ==1 ) return 2;
  if (nmuons >0 && (nelectrons + npions + npi0 ) == 0 && nprotons ==2 ) return 3;
  if (nmuons >0 && (nelectrons + npions + npi0 ) == 0 && nprotons >2 ) return 4;
  if (nmuons >0 && (nelectrons ) == 0 && (npions + npi0) == 1 ) return 5;
  if (nmuons >0 && (nelectrons ) == 0 && (npions + npi0) > 1 ) return 6;
  if (nmuons == 0 && nelectrons > 0 ) return 7;
  if (nmuons==0  && nelectrons ==0) return 8;
//  if (OOFVflag) return 9;
//  if (cosmicflag) return 10;  //check with colton how to select cosmic event

  else return 11;


}
void SimpleAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet)
{
  //art::ServiceHandle<cheat::BackTracker> bt;	
  std::map<int,double> trkID_E;	
  for(size_t j = 0; j < track_hits.size(); ++j)
  {	
    art::Ptr<recob::Hit> hit = track_hits[j];
    //const auto& hit = *track_hits[j];
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
    for(size_t k = 0; k < TrackIDs.size(); k++)
    {	
      trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
    }
  }
  //std::cout<<"[SimpleAna] truthMatcher :"<<trkID_E.size()<<std::endl;
  double E_em =0.0;	
  double max_E = -999.0;	
  double total_E = 0.0;	
  int TrackID = -999;	
  double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla	
  //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 	
  //!since we are looking for muons/pions/protons this should be enough 	
  if( !trkID_E.size() ) 
  {
    MCparticle = 0;	
    return; //Ghost track???	
  }	
  for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii)
  {	
    total_E += ii->second;	
    if((ii->second)>max_E)
    {	
      partial_E = ii->second;
      max_E = ii->second;
      TrackID = ii->first;
      if( TrackID < 0 ) E_em += ii->second;
    }	
  }
   	
  //MCparticle = bt->TrackIDToParticle(TrackID);		
  MCparticle = fMCTruthMatching->TrackIDToParticle(TrackID);		
  //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of 	
  //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle	
  //we don't want to track gammas or any other EM activity 	
  if( TrackID < 0 ) return;		
  //Efrac = (partial_E+E_em)/total_E;	
  Efrac = (partial_E)/total_E;		
  //completeness	
  double totenergy =0;	
  for(size_t k = 0; k < all_hits.size(); ++k)
  {	  
    art::Ptr<recob::Hit> hit = all_hits[k];	
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);	
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);	
    for(size_t l = 0; l < TrackIDs.size(); ++l)
    {	
      if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;	
    }	
  } 	
  Ecomplet = partial_E/totenergy;
}
double SimpleAna::GetEnergy(std::vector<art::Ptr<anab::Calorimetry>> calos) {
  double result=-9999;

  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
    std::cout << "old Energy =" << c->KineticEnergy() << "new result = ";  
    std::vector<double> dEdx_v = c->dEdx(); 
    std::vector<double> TrkPtch_v = c->TrkPitchVec(); 
//
    if (dEdx_v.size() == 0){
      return result;
    }
    double totE=0;
    for (unsigned int i(0); i < dEdx_v.size();++i){
      totE += dEdx_v.at(i) * TrkPtch_v.at(i);
    }
    result = totE;
    std::cout << totE << std::endl;
  }
  return result;
}


std::vector<double> SimpleAna::GetRR(std::vector<art::Ptr<anab::Calorimetry>> calos) {
  std::vector<double> result;
  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
   
    std::vector<double> RR_v = c->ResidualRange(); 
    return RR_v;
  }
  return result;
}


std::vector<double> SimpleAna::GetdEdx(std::vector<art::Ptr<anab::Calorimetry>> calos) {
  std::vector<double> result;
  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
   
    std::vector<double> dEdx_v = c->dEdx(); 
//    if (dEdx_v.size() == 0){
//      return 0;
//    }
    return dEdx_v;
  }
  return result;
}



void SimpleAna::analyze(art::Event const & e)
{

  bool isMC=!e.isRealData();
  if(isMC)  fMCTruthMatching->Rebuild(e);
  
  std::cout << "[SimpleAna] Simple UBXSec Validation Module. Starts." << std::endl;

  // Implementation of required member function here.
  art::Handle<std::vector<ubana::SelectionResult>> selection_h;
  e.getByLabel("UBXSec",selection_h);
  if (!selection_h.isValid() || selection_h->empty()) {
    std::cout << "[SimpleAna] SelectionResult handle is not valid or empty." << std::endl;
  }
  //get the selected vertex here and fill selected vertex and selected handle
  std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
  art::fill_ptr_vector(selection_v, selection_h);

  if (selection_v.at(0)->GetSelectionStatus()) {
    std::cout << "[SimpleAna] Event is selected" << std::endl;
  } else {
    std::cout << "[SimpleAna] Event is not selected" << std::endl;
    std::cout << "[SimpleAna] Failure reason " << selection_v.at(0)->GetFailureReason()  << std::endl;
  }


  //get the failure map here check if it past the event selection
  std::map<std::string,bool> failure_map = selection_v.at(0)->GetCutFlowStatus();

  std::cout << "[SimpleAna] Now Printing Cut Flow Status" << std::endl;

  for (auto iter : failure_map) {
    std::cout << "[SimpleAna] Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
  }

  _ccinclflag=selection_v.at(0)->GetSelectionStatus();
  if(!_ccinclflag) return;
  //===================================================================================================================

  std::cout<<"Start Getting the TPC object<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  //e.getByLabel("TPCObjectMaker", tpcobj_h);

 

  
  if (!tpcobj_h.isValid()) {
    std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
  }

  std::cout<<"Start Getting the tracks from PF particles<<<<<<<<<<<<<"<<std::endl;
 

  // Get Tracks 
  art::Handle<std::vector<recob::Track>> track_h;

  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }
  std::vector<art::Ptr<recob::Track>> track_p_v;
  art::fill_ptr_vector(track_p_v, track_h);

  art::FindManyP<recob::OpFlash> opfls_ptr_coll_v(track_h, e, _acpt_producer);
  art::FindManyP<recob::PFParticle> pfp_from_track(track_h, e, _pfp_producer);
  art::FindManyP<anab::Calorimetry> calos_from_track(track_h, e, _calorimetry_producer);






  // Get Tracks
  /*art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }
  */
  std::cout<<"Start Getting the PF Particles<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  // Get PFP
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  if(!pfp_h.isValid()){
       std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
       //throw std::exception();
  }
  if(pfp_h->empty()) {
       std::cout << "[UBXSec] PFP " << _pfp_producer << " is empty." << std::endl;
  }
   
  //===================================================================================== 
                        
  //get tracks from pfp
  /*std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfp_h);
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);
  std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
  */
  /*ubxsec_event->n_pfp = ubxsec_event->n_pfp_primary = 0;
  for (size_t i = 0; i < pfp_h->size(); i++) {
    ubxsec_event->n_pfp++;
    if ((*pfp_h)[i].IsPrimary())
      ubxsec_event->n_pfp_primary++;
  }
  */
  //==================================================================================== 
  std::cout<<"Start Getting the associated TPC object<<<<<<<<<<<<<<<<<<"<<std::endl;

  // if the event is selected, get the associated TPC object here
  art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, e, "UBXSec");
  art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);


  std::cout<<"Start Getting the tracks from associated TPC object<<<<<<<<<<<<<<<<<<"<<std::endl;

  // get the tracks 
  auto const& TPCobjToTrack = *(e.getValidHandle<art::Assns<recob::Track, ubana::TPCObject>>("TPCObjectMaker"));
  for (auto const& TPCobjAndTrack: TPCobjToTrack) {
    std::cout << "TPCObject: " << TPCobjAndTrack.second << " -- Track: " << TPCobjAndTrack.first << std::endl;
  }
  art::FindManyP<recob::Track> tracks_from_tpcobject(tpcobj_h, e, "TPCObjectMaker");
  std::vector<art::Ptr<recob::Track>> tracks = tracks_from_tpcobject.at(tpcobj_candidate.key());

  std::cout<<"Start Getting the pfps from associated TPC object<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  // get the PFParticles
  art::FindManyP<recob::PFParticle> pfps_from_tpcobject(tpcobj_h, e, "TPCObjectMaker");
  std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_tpcobject.at(tpcobj_candidate.key());

  std::cout<<"Start Getting the neutrino vertex from associated TPC object<<<<<<<<<<<<<<<<<<"<<std::endl;

  //get neutrino vertex
  art::FindManyP<recob::Vertex> vertices_from_tpcobject(tpcobj_h, e, "TPCObjectMaker");
  art::Ptr<recob::Vertex> vertex = vertices_from_tpcobject.at(tpcobj_candidate.key()).at(0); 
  //std::vector<art::Ptr<recob::Vertex>> vertex_test = vertices_from_tpcobject.at(tpcobj_candidate.key()); 


  std::cout<<"Start Getting the tracks form PF particles<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);
  //std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
 
  //Get PID information
  art::FindManyP<anab::ParticleID> particleids_from_track(track_h, e, _particle_id_producer);
  if (!particleids_from_track.isValid()) {
    std::cout << "[SimpleAna] anab::ParticleID is not valid." << std::endl;
  }
  std::cout << "[SimpleAna] Numeber of particleids_from_track " << particleids_from_track.size() << std::endl;



  //get the hit of tracks
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;


  if(e.getByLabel(_hit_producer,hitListHandle))
  art::fill_ptr_vector(hitlist, hitListHandle);
  std::cout<<"checkc it the hitlist handle is valid or not "<<hitListHandle.isValid()<<std::endl;

 

  art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> fmhitmc(hitListHandle,e,"crHitRemovalTruthMatch");


  auto const& TrackToHits = *(e.getValidHandle<art::Assns<recob::Hit, recob::Track>>("pandoraNu::UBXSec"));
  for (auto const& TrackAndHit: TrackToHits ) {
    std::cout << "Hit: " << TrackAndHit.first << " -- Track: " << TrackAndHit.second << std::endl;
  }

  art::FindManyP<recob::Hit> fmht (track_h, e, "pandoraNu::UBXSec");
  //art::FindMany<recob::Hit> fmh(track_h, e, "pandoraNu::UBXSec");
  
  std::vector<art::Handle<std::vector<recob::Hit>>> allHitCollections;
  e.getManyByType(allHitCollections);
  for (auto const& hit_h: allHitCollections) {
    std::cout << "Hits from '"
      << hit_h.provenance()->processName() << ":"
      << hit_h.provenance()->moduleLabel() << ":"
      << hit_h.provenance()->productInstanceName()
      << "' have ID " << hit_h.id() << std::endl;
  }
  
  for (std::string label: { 
               "gaushitTruthMatch",
               "crHitRemovalTruthMatch",
               "trajclusterTruthMatch"
    })
  {
    std::cout << "Assn '" << label << "'" << std::endl;
    try {
      auto const& HitToMCParticle = *(e.getValidHandle<art::Assns<recob::Hit, simb::MCParticle, anab::BackTrackerHitMatchingData>>(label));
      for (auto const& HitAndMCParticle: HitToMCParticle ) {
        std::cout << "Hit: " << HitAndMCParticle.first << " -- MCParticle: " << HitAndMCParticle.second << std::endl;
      }
    }
    catch(art::Exception const& e) {
      std::cout << " ==> not present" << std::endl;
    }
  }
  std::cout<<"Start the CC1uNP Selection<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

  art::Handle<std::vector<recob::MCSFitResult>>   MCSFitHandle;
  e.getByLabel(fTrackMCSFitLabel,     MCSFitHandle);
  std::vector<art::Ptr<recob::MCSFitResult> >  mcsfitlist;
  art::fill_ptr_vector(mcsfitlist, MCSFitHandle);

  //==================================================================================================================

  //reco true matching
  
  // Get MCGhosts from the event
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel("RecoTrueMatcher",ghost_h);
  if(!ghost_h.isValid()){
     mf::LogError(__PRETTY_FUNCTION__) << "MCGhost product not found." << std::endl;
      //throw cet::exception();
  }
  std::cout<<"[SimpleAna] Start Looping Over PFParticles<<<<<<<<<"<<std::endl;

  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer); 
  

  _npfps=0;
  _npfp_tracks=0;
  _npfp_showers=0;
  _pfp_ntrk2flag=false;


  pfp_isTrack->clear();
  pfp_isShower->clear();
  pfp_isPrimary->clear();
  pfp_ndaughters->clear();
 
  pfp_pdg->clear();
  pfp_origin->clear();
  pfp_status->clear();
  pfp_parId->clear();
  pfp_theta->clear();
  pfp_costheta->clear();
  pfp_phi->clear();
  pfp_mom->clear();
  pfp_startx->clear();
  pfp_starty->clear();
  pfp_startz->clear();
  pfp_endx->clear();
  pfp_endy->clear();
  pfp_endz->clear();
  pfp_endE->clear();
  pfp_KE->clear();
  pfp_Mass->clear();


  track_pfp_upflag->clear();
  track_pfp_Id->clear();
  track_pfp_length->clear();
  track_pfp_theta->clear();
  track_pfp_costheta->clear();
  track_pfp_phi->clear();
  track_pfp_startx->clear();
  track_pfp_starty->clear();
  track_pfp_startz->clear();
  track_pfp_endx->clear();
  track_pfp_endy->clear();
  track_pfp_endz->clear();
  track_pfp_Mom->clear(); 
  track_pfp_trunmeandqdx->clear();
  track_pfp_trunmeandqdx_U->clear();
  track_pfp_trunmeandqdx_V->clear();

  track_pfp_pida->clear();
  track_pfp_nhits->clear();
  track_pfp_dEdx->clear();
  track_pfp_RR->clear();

  trkf::TrackMomentumCalculator trkm_pfp;
  int TrackpfpCandidate=-999;
  float TrackpfpCandLength=-999.0;


  for(auto pfp : pfps){

    
    const std::vector<size_t> &daughterIDs=pfp->Daughters();


    pfp_isTrack->push_back(lar_pandora::LArPandoraHelper::IsTrack(pfp));
    pfp_isShower->push_back(lar_pandora::LArPandoraHelper::IsShower(pfp));
 
    pfp_isPrimary->push_back(pfp->IsPrimary());

    pfp_ndaughters->push_back(daughterIDs.size());

    //if(daughterIDs.size() !=0) continue;

    if(lar_pandora::LArPandoraHelper::IsTrack(pfp)) {_npfp_tracks=_npfp_tracks+1;}
    if(lar_pandora::LArPandoraHelper::IsShower(pfp)) {_npfp_showers=_npfp_showers+1;}

    if(lar_pandora::LArPandoraHelper::IsShower(pfp)) continue;
    
    _npfps=_npfps+1;

    if(isMC){
    std::vector<art::Ptr<ubana::MCGhost>> mcghosts=mcghost_from_pfp.at(pfp.key());
    std::vector<art::Ptr<simb::MCParticle>> mcpars;

    
    if(mcghosts.size()==0 ||mcghosts.size()>1){
         std::cout<<"[SimpleAna] \t\t mcghosts is either 0 or >1"<<std::endl;
         continue;
    }
    mcpars=mcpar_from_mcghost.at(mcghosts[0].key());
    const auto mc_truth=UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpars[0]->TrackId());
    if(!mc_truth) {
      std::cerr<<"[SimpleAna] Problem with MCTruth pointer.  "<<std::endl;
    }
    pfp_pdg->push_back(mcpars[0]->PdgCode());
    pfp_origin->push_back(mc_truth->Origin());
    pfp_status->push_back(mcpars[0]->StatusCode());
    pfp_parId->push_back(mcpars[0]->TrackId());
    pfp_theta->push_back(mcpars[0]->Momentum().Theta());
    pfp_costheta->push_back(TMath::Cos(mcpars[0]->Momentum().Theta()));
    pfp_phi->push_back(mcpars[0]->Momentum().Phi());
    pfp_mom->push_back(mcpars[0]->P());   
    pfp_startx->push_back(mcpars[0]->Vx());
    pfp_starty->push_back(mcpars[0]->Vy());
    pfp_startz->push_back(mcpars[0]->Vz());
    pfp_endx->push_back(mcpars[0]->EndX());
    pfp_endy->push_back(mcpars[0]->EndY());
    pfp_endz->push_back(mcpars[0]->EndZ());
    pfp_endE->push_back(mcpars[0]->EndE());
    pfp_KE->push_back(mcpars[0]->E()-mcpars[0]->Mass());
    pfp_Mass->push_back(mcpars[0]->Mass());

    //int TrackpfpCandidate=-999;
    //float TrackpfpCandLength=-999.0;
    }//end of if isMC in pfp loop
    std::vector<art::Ptr<recob::Track>> tracks_pfp=tracks_from_pfp.at(pfp.key());
    std::cout<<"[SimpleAna] \t\t  n tracks ass to this pfp: "<<tracks_pfp.size()<<std::endl;
    //loop over all the tracks ass to this pfp
    for(auto track_pfp : tracks_pfp){
        std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track_pfp.key());
        if(pids.size()==0) std::cout<<"[SimpleAna] \t\t Zero Particle ID  "<<std::endl;
        if(pids.size()>1) std::cout<<"[SimpleAna] \t\t Particle ID Vector is bigger than 1. Only 1 saved.   "<<std::endl;
        for(auto pid : pids){
            if(!pid->PlaneID().isValid) continue;
            int planenum=pid->PlaneID().Plane;
            if(planenum<0 || planenum>2) continue;
            /*if(mcpars[0]->PdgCode()==13){ 
                _h_pida_muon->Fill(pid->PIDA());
            } else if(mcpars[0]->PdgCode()==2212)
                _h_pida_proton->Fill(pid->PIDA());
            } 
            */
            track_pfp_pida->push_back(pid->PIDA());
        }
        TVector3  trackPos_pfp=track_pfp->Vertex();
        TVector3  trackEnd_pfp=track_pfp->End();
        track_pfp_Id->push_back(track_pfp->ID());
        track_pfp_length->push_back(track_pfp->Length());
        track_pfp_theta->push_back(track_pfp->Theta());
        track_pfp_costheta->push_back(TMath::Cos(track_pfp->Theta()));
        track_pfp_phi->push_back(track_pfp->Phi());
        track_pfp_startx->push_back(trackPos_pfp.X());
        track_pfp_starty->push_back(trackPos_pfp.Y());
        track_pfp_startz->push_back(trackPos_pfp.Z());
        track_pfp_endx->push_back(trackEnd_pfp.X());
        track_pfp_endy->push_back(trackEnd_pfp.Y());
        track_pfp_endz->push_back(trackEnd_pfp.Z());
        track_pfp_Mom->push_back(track_pfp->VertexMomentum());

        std::vector<art::Ptr<anab::Calorimetry>> calos_track_pfp=calos_from_track.at(track_pfp.key());

        track_pfp_trunmeandqdx->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos_track_pfp));
        track_pfp_trunmeandqdx_U->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos_track_pfp, 0));
        track_pfp_trunmeandqdx_V->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos_track_pfp, 1));  
        track_pfp_upflag->push_back(_muon_finder.MIPConsistency(UBXSecHelper::GetDqDxTruncatedMean(calos_track_pfp), track_pfp->Length()));
        track_pfp_dEdx->push_back(GetdEdx(calos_track_pfp));
        track_pfp_RR->push_back(GetRR(calos_track_pfp));


        std::vector<art::Ptr<recob::Hit>> allKHits_pfp = fmht.at(track_pfp.key());
        track_pfp_nhits->push_back(allKHits_pfp.size()); 

        //if(track_pfp->Length()>TrackpfpCandLength){
        //   TrackpfpCandidate=track_pfp->ID();
        //   TrackpfpCandLength=track_pfp->Length;
        //   track_pfp_mom_mucand=trkm_pfp.GetTrackMomentum(TrackpfpCandLength, 13);
        //}
    }//end of loop over the tracks associated to this pfparticles    
  } //end of loop overr all the pf particles 

  for(size_t ncand_pfp=0; ncand_pfp<track_pfp_Id->size(); ncand_pfp++){
    if(track_pfp_length->at(ncand_pfp)>TrackpfpCandLength){
    TrackpfpCandidate=track_pfp_Id->at(ncand_pfp);
    TrackpfpCandLength=track_pfp_length->at(ncand_pfp);
    track_pfp_mom_mucand=trkm_pfp.GetTrackMomentum(TrackpfpCandLength, 13);
    if(!isContained(track_pfp_endx->at(ncand_pfp), track_pfp_endy->at(ncand_pfp), track_pfp_endz->at(ncand_pfp))){
               track_pfp_mom_mucand=mcsfitlist[track_pfp_Id->at(ncand_pfp)]->bestMomentum();
    }
    }
  }

  _pfp_pdqdxflag=true;
  _pfp_upinFVflag=true;
  track_pfp_mom_pcand->clear();
  for(size_t mcand_pfp=0; mcand_pfp<track_pfp_Id->size(); mcand_pfp++){
    if(track_pfp_Id->at(mcand_pfp)==TrackpfpCandidate) continue;
    track_pfp_mom_pcand->push_back(trkm_pfp.GetTrackMomentum(TrackpfpCandLength, 2212));
    if(track_pfp_upflag->at(mcand_pfp)==1) {_pfp_pdqdxflag=false;}
    if(!isContained(track_pfp_startx->at(mcand_pfp), track_pfp_starty->at(mcand_pfp), track_pfp_startz->at(mcand_pfp)) ||
       !isContained(track_pfp_endx->at(mcand_pfp), track_pfp_endy->at(mcand_pfp), track_pfp_endz->at(mcand_pfp))) 
       {_pfp_upinFVflag=false;}
  }

  _pfp_ntrk2flag=false;
  if(_npfp_tracks>=2) {_pfp_ntrk2flag=true;}
  //--------------------------------------------------------------  
  //==============================================================
 
  
  std::cout<<"Start Getting the MC truth information<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
  

  //
  // MCTruth
  //
   
  bool in_tpcactive = false;
  bool energy_range = false;
  bool right_flavour = false;
  if(isMC){
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (e.getByLabel("generator",mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);
  int iList = 0;

  
  ::geoalgo::Vector truth_nu_vtx (mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz());
  ::art::ServiceHandle<geo::Geometry> geo;
  ::geoalgo::AABox tpcvol(0, (-1.)*(geo->DetHalfHeight()), 0.,
              geo->DetHalfWidth()*2, geo->DetHalfHeight(), geo->DetLength());

  if(tpcvol.Contain(truth_nu_vtx)) in_tpcactive = true;
  else in_tpcactive = false;

  if (mclist[iList]->GetNeutrino().Nu().E() > 0.05 && mclist[iList]->GetNeutrino().Nu().E() < 1.5)
    energy_range = true;
  else 
    energy_range = false;

  if (mclist[iList]->GetNeutrino().CCNC() == 0 && mclist[iList]->GetNeutrino().Nu().PdgCode() == 14)
    right_flavour = true;
  else 
    right_flavour = false;

  
 
  _nu_energy = mclist[iList]->GetNeutrino().Nu().E();

  _lep_mom=mclist[iList]->GetNeutrino().Lepton().P();
  _lep_costheta=mclist[iList]->GetNeutrino().Lepton().Pz()/mclist[iList]->GetNeutrino().Lepton().P();
  _lep_phi   = UBXSecHelper::GetPhi(mclist[iList]->GetNeutrino().Lepton().Px(),
                                                   mclist[iList]->GetNeutrino().Lepton().Py(),
                                                   mclist[iList]->GetNeutrino().Lepton().Pz());

  _nu_ccnc   = mclist[iList]->GetNeutrino().CCNC();
  _nu_mode   = mclist[iList]->GetNeutrino().Mode();
  _nu_pdg    = mclist[iList]->GetNeutrino().Nu().PdgCode();
 
  _fTruenuvrtxx= mclist[iList]->GetNeutrino().Nu().Vx();
  _fTruenuvrtxy= mclist[iList]->GetNeutrino().Nu().Vy(); 
  _fTruenuvrtxz= mclist[iList]->GetNeutrino().Nu().Vz();
 

   std::cout<<"[SimpleAna] in_tpcactive?<<<<<<<<<<<<<<"<<in_tpcactive<<std::endl;
   std::cout<<"[SimpleAna] energy_range?<<<<<<<<<<<<<<"<<energy_range<<std::endl;
   std::cout<<"[SimpleAna] right_flavour?<<<<<<<<<<<<<<"<<right_flavour<<std::endl;

   Int_t nGeniePrimaries=0;
   if (mclist[iList]->NeutrinoSet()) nGeniePrimaries = mclist[iList]->NParticles();

   trueProtonsTrueMomentum->clear();
   trueProtonsTrueTheta->clear();
   trueProtonsTruePhi->clear();
   trueProtonsEndMomentum->clear();
   trueProtonsEndProcess->clear();

   Nmuons_from_genie=0;
   Nelectrons_from_genie=0;
   Npionpms_from_genie=0;
   Npion0s_from_genie=0;
   Nprotons_from_genie=0;
   Nothers_from_genie=0;
   std::string parpri("primary");

   Nprotons_from_genie_200thresh=0;
   Nprotons_from_genie_300thresh=0;
   Nprotons_from_genie_400thresh=0;
   std::cout<<"[SimpleAna] Total number of GENIE parimary particles is "<<nGeniePrimaries<<std::endl;
   for (int igeniepart(0); igeniepart<nGeniePrimaries; igeniepart++){
      simb::MCParticle part = mclist[iList]->GetParticle(igeniepart);
      bool parisPrimary=part.Process()==parpri;
      //std::cout<<"parisPrimary ???  "<<parisPrimary<<std::endl;
      //std::cout<<"mclist origin  is ???"<<mclist[iList]->Origin()<<std::endl;      
      //std::cout<<"mother particle id is ??? "<<part.Mother()<<std::endl;


      if(parisPrimary && part.StatusCode()==1 &&mclist[iList]->Origin()==simb::kBeamNeutrino){
        if(_nu_ccnc==0 && abs(part.PdgCode())==13) {Nmuons_from_genie=Nmuons_from_genie+1;} 
        if(_nu_ccnc==0 && abs(part.PdgCode())==211 ) {Npionpms_from_genie=Npionpms_from_genie+1;} 
        if(_nu_ccnc==0 && abs(part.PdgCode())==111 ) {Npion0s_from_genie=Npion0s_from_genie+1;} 
        if(_nu_ccnc==0 && abs(part.PdgCode())==11 ) {Nelectrons_from_genie=Nelectrons_from_genie+1;} 
        if(_nu_ccnc==0 && abs(part.PdgCode())==2212 ){
          Nprotons_from_genie=Nprotons_from_genie+1;
          trueProtonsTrueMomentum->push_back(part.P());
          trueProtonsTrueTheta->push_back(part.Momentum().Theta());
          trueProtonsTruePhi->push_back(part.Momentum().Phi());
          trueProtonsEndMomentum->push_back(part.EndMomentum().P());
          trueProtonsEndProcess->push_back(part.EndProcess());
        }
          //check the proton's momentum here 
          //if(part.Momentum().Vect().Mag()>0.2){Nprotons_from_genie_200thresh=Nprotons_from_genie_200thresh+1;}
          //if(part.Momentum().Vect().Mag()>0.3){Nprotons_from_genie_300thresh=Nprotons_from_genie_300thresh+1;}
          //if(part.Momentum().Vect().Mag()>0.4){Nprotons_from_genie_400thresh=Nprotons_from_genie_400thresh+1;}
        if(_nu_ccnc==0 && abs(part.PdgCode())==2212 &&part.P()>0.2){Nprotons_from_genie_200thresh=Nprotons_from_genie_200thresh+1;}
        if(_nu_ccnc==0 && abs(part.PdgCode())==2212 &&part.P()>0.3){Nprotons_from_genie_300thresh=Nprotons_from_genie_300thresh+1;}
        if(_nu_ccnc==0 && abs(part.PdgCode())==2212 &&part.P()>0.4){Nprotons_from_genie_400thresh=Nprotons_from_genie_400thresh+1;}
         
        //else {Nothers_from_genie=Nothers_from_genie+1;}
      } //end of if part.StatusCode==1)
   } //end of loop over all the genie particles
    /// Also here we should get things like the true struck neutron momentum - I think we need a GTruth object for this
   art::ValidHandle< std::vector<simb::GTruth> > gtruth = e.getValidHandle< std::vector<simb::GTruth> >("generator");
   if (gtruth->size() <1){
     std::cout << "WARNING NO GTRUTH OBJECT" << std::endl;
   }
   else{
     TLorentzVector v_tmp = gtruth->at(0).fHitNucP4;
     fHitNucP4->SetXYZT(v_tmp.X(), v_tmp.X(), v_tmp.Y(), v_tmp.E());
   }
   //======================================================================================
   
   Int_t TopFlag_genie=Topology(Nmuons_from_genie, Nelectrons_from_genie, Npionpms_from_genie, Npion0s_from_genie, Nprotons_from_genie);
   Int_t TopFlag_genie_200thresh=Topology(Nmuons_from_genie, Nelectrons_from_genie, Npionpms_from_genie, Npion0s_from_genie, Nprotons_from_genie_200thresh);
   Int_t TopFlag_genie_300thresh=Topology(Nmuons_from_genie, Nelectrons_from_genie, Npionpms_from_genie, Npion0s_from_genie, Nprotons_from_genie_300thresh);
   Int_t TopFlag_genie_400thresh=Topology(Nmuons_from_genie, Nelectrons_from_genie, Npionpms_from_genie, Npion0s_from_genie, Nprotons_from_genie_400thresh);


   
  std::cout<<"[SimpleAna] Get the true struck neutron momentum"<<std::endl;
  Int_t nmuons=0;
  Int_t npions=0;
  Int_t npi0=0;
  Int_t nprotons=0;
  Int_t nprotons_200thresh=0;
  Int_t nprotons_300thresh=0;
  Int_t nprotons_400thresh=0;
  Int_t nelectrons=0;
  OOFVflag=false; //nu event outside FV




    
  art::Ptr<simb::MCTruth> mctruth; 
  mctruth = mclist[0];



  std::cout<<"[SimpleAna] Total number of G4 particles is  "<< mclist[iList]->NParticles()<<std::endl;

  art::Handle<std::vector<simb::MCParticle>> mcparticleListHandle;
  std::vector<art::Ptr<simb::MCParticle>> mcparticlelist;

  if (e.getByLabel("largeant", mcparticleListHandle))
  art::fill_ptr_vector(mcparticlelist, mcparticleListHandle);


  fhg4parpdg->clear();
  fhg4parstatus->clear();
  fhg4parpx->clear();
  fhg4parpy->clear();
  fhg4parpz->clear();
  fhg4partheta->clear();
  fhg4parphi->clear();
  fhg4parp->clear();


  _fg4processname->clear();
  _fg4Mother->clear();
  _fg4TrackId->clear();
  
  _fg4parpdg->clear();
  _fg4parstatus->clear();
  _fg4EngS->clear();
  _fg4EngE->clear();
  _fg4Mass->clear();
  _fg4px->clear();
  _fg4py->clear();
  _fg4pz->clear();

  _fg4StartPointx->clear();
  _fg4StartPointy->clear();
  _fg4StartPointz->clear();
 
  _fg4EndPointx->clear();
  _fg4EndPointy->clear();
  _fg4EndPointz->clear();

  _fg4EndT->clear();
  _fg4theta->clear();
  _fg4phi->clear();

  _fg4theta_xz->clear();
  _fg4theta_yz->clear();
  _fg4origin->clear();
  _fg4MCTruthIndex->clear();
  _fg4NumberDaughters->clear();

  std::string pri("primary");

   //check the  topology and classfiy the events into different interactions
  for(size_t p=0; p< mcparticleListHandle->size(); p++){
         
        if(mcparticleListHandle.isValid() && mcparticleListHandle->size()>0)
        {
        simb::MCParticle const & mc_par = mcparticleListHandle->at(p);

        const art::Ptr<simb::MCTruth> mc_truth=fMCTruthMatching->TrackIDToMCTruth(mc_par.TrackId());

        Bool_t isPrimary=mc_par.Process() == pri;
        //std::cout<<"isPrimary"<<mc_par.Process()<<std::endl;

        _fg4processname->push_back(mc_par.Process()); 
        _fg4Mother->push_back(mc_par.Mother());
        _fg4TrackId->push_back(mc_par.TrackId());


        _fg4parpdg->push_back(mc_par.PdgCode());
        _fg4parstatus->push_back(mc_par.StatusCode());
        _fg4EngS->push_back(mc_par.E());
        _fg4EngE->push_back(mc_par.EndE());
        _fg4Mass->push_back(mc_par.Mass());
        _fg4px->push_back(mc_par.Px());
        _fg4py->push_back(mc_par.Py());
        _fg4pz->push_back(mc_par.Pz());
        _fg4P->push_back(mc_par.Momentum().Vect().Mag());
        _fg4StartPointx->push_back(mc_par.Vx());
        _fg4StartPointy->push_back(mc_par.Vy());
        _fg4StartPointz->push_back(mc_par.Vz());
        _fg4EndPointx->push_back(mc_par.EndPosition()[0]);
        _fg4EndPointy->push_back(mc_par.EndPosition()[1]);
        _fg4EndPointz->push_back(mc_par.EndPosition()[2]);
        
        //--------------------------------------------------

        _fg4EndT->push_back(  mc_par.EndT());
        _fg4theta->push_back( mc_par.Momentum().Theta());
        _fg4phi->push_back( mc_par.Momentum().Phi());
        _fg4theta_xz->push_back( std::atan2(mc_par.Px(), mc_par.Pz()));
        _fg4theta_yz->push_back( std::atan2(mc_par.Py(), mc_par.Pz()));

        _fg4origin->push_back(mc_truth->Origin());  //what created particle
        _fg4MCTruthIndex->push_back(mc_truth.key());
        _fg4NumberDaughters->push_back(mc_par.NumberDaughters());

      if( isPrimary && mc_truth->NeutrinoSet() && mc_par.StatusCode()==1 && mc_par.Mother()==0 && mc_truth->Origin()== simb::kBeamNeutrino)  
      //primary tells you if the particle is from Michel electron or decay of other particle
      { 
       
 
        fhg4parpdg->push_back(mc_par.PdgCode());
        fhg4parstatus->push_back(mc_par.StatusCode());
        fhg4parpx->push_back(mc_par.Px());
        fhg4parpy->push_back(mc_par.Py());
        fhg4parpz->push_back(mc_par.Pz());
        fhg4partheta->push_back(mc_par.Momentum().Theta());
        fhg4parphi->push_back(mc_par.Momentum().Phi());
        fhg4parp->push_back(mc_par.Momentum().Vect().Mag());


        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==13) {
          nmuons=nmuons+1;
        }
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==211) {npions=npions+1;}
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==111) {npi0=npi0+1;}
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==11) {nelectrons=nelectrons+1;}
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==2212) {nprotons=nprotons+1;}

        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==2212 && mc_par.Momentum().Vect().Mag() > 0.2 ) {nprotons_200thresh=nprotons_200thresh+1;}
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==2212 && mc_par.Momentum().Vect().Mag() > 0.3 ) {nprotons_300thresh=nprotons_300thresh+1;}
        if(_nu_ccnc==0 && abs(mc_par.PdgCode())==2212 && mc_par.Momentum().Vect().Mag() > 0.4 ) {nprotons_400thresh=nprotons_400thresh+1;}
      }
      }//end of if the mcparticle list handle is valid
   }//end of loop over all the geant 4 particles
   //redefine the OOFV here
   if(!_fiducial_volume.InFV(_fTruenuvrtxx, _fTruenuvrtxy, _fTruenuvrtxz)) {OOFVflag=true;} 
   
   std::cout<<"[SimpleAna] End of loop over all the G4 particles"<<std::endl;

   Int_t TopFlag=Topology(nmuons, nelectrons, npions, npi0, nprotons);
   Int_t TopFlag_200thresh=Topology(nmuons, nelectrons, npions, npi0, nprotons_200thresh);
   Int_t TopFlag_300thresh=Topology(nmuons, nelectrons, npions, npi0, nprotons_300thresh);
   Int_t TopFlag_400thresh=Topology(nmuons, nelectrons, npions, npi0, nprotons_400thresh);
  
   truthtop=TopFlag;
   truthtop_200thresh=TopFlag_200thresh;
   truthtop_300thresh=TopFlag_300thresh;
   truthtop_400thresh=TopFlag_400thresh;

   truthtop_genie=TopFlag_genie;
   truthtop_genie_200thresh=TopFlag_genie_200thresh;
   truthtop_genie_300thresh=TopFlag_genie_300thresh;
   truthtop_genie_400thresh=TopFlag_genie_400thresh;


        
   std::cout<<"Topology of this event is "<<truthtop_genie<<" "<<truthtop_genie_200thresh<<" "<<truthtop_genie_300thresh<<" "<<truthtop_genie_400thresh<<std::endl;
   std::cout<<"Topology of this event is "<<truthtop<<" "<<truthtop_200thresh<<" "<<truthtop_300thresh<<" "<<truthtop_400thresh<<std::endl;




  }//end of if isMC
  //---------------------------------------------------------------------------------
  // EventWeight
  //
  
  double bnb_weight = 1.;

  if(isMC){
  art::Handle<std::vector<evwgh::MCEventWeight>> eventweight_h;
  e.getByLabel("eventweight", eventweight_h);
  if(!eventweight_h.isValid()){
    std::cout << "[SimpleAna] MCEventWeight product not found..." << std::endl;
  }
  std::vector<art::Ptr<evwgh::MCEventWeight>> eventweight_v;
  art::fill_ptr_vector(eventweight_v, eventweight_h);

  if (eventweight_v.size() > 0) {
    art::Ptr<evwgh::MCEventWeight> evt_wgt = eventweight_v.at(0);
    for (auto entry : evt_wgt->fWeight) {
      if (entry.first.find("bnbcorrection") != std::string::npos) {
        bnb_weight *= entry.second.at(0);
        std::cout << "[SimpleAna] BNB Correction Weight: " << bnb_weight << std::endl;
      }
    }
  }

  } //end of if isMC for the bnb correction weight
  //---------------------------------------------------------------------------------
  //get the vertex position
  Double_t xyz[3]={};
  //int Naux=0;
  //for(auto vertex : vertex_test){
  //Naux=Naux+1;
  //if(Naux>0) continue; 
  vertex->XYZ(xyz);
  //}

  _nuvtxx_reco=xyz[0];
  _nuvtxy_reco=xyz[1];
  _nuvtxz_reco=xyz[2];
  TVector3 nuvertexPos(xyz[0],xyz[1],xyz[2]);


  //do the CC1uNP event selection here
  //#1 total number of tracks cut


  upflag->clear();
  trackmom_pcand->clear();


  trackId->clear();
  tracklength->clear();
  trackmom->clear();
  trackcostheta->clear();
  trackphi->clear();


  trackstartx->clear();
  trackstarty->clear();
  trackstartz->clear();
  
 

  trackendx->clear();
  trackendy->clear();
  trackendz->clear();

  tracktrunmeandqdx->clear();
  tracktrunmeandqdx_U->clear();
  tracktrunmeandqdx_V->clear();

  tracknhits->clear();

  trackcand_origin->clear();
  trackcand_nuset->clear();
  trackcand_parPDG->clear();
  trackcand_parStatusCode->clear();
  trackcand_parTheta->clear();                  
  trackcand_parCosTheta->clear();
  trackcand_parSinTheta->clear();                  
  trackcand_parE->clear();        
  trackcand_parMass->clear();
  trackcand_parKE->clear();
  trackcand_parEndE->clear();
  trackcand_parPx->clear();
  trackcand_parPy->clear();
  trackcand_parPz->clear();
  trackcand_parPhi->clear();
  trackcand_parCosPhi->clear();
  trackcand_parSinPhi->clear();


  //Do the vertex track association
  int NumTracksNearVertex=0;

  int TrackCandidate=-999;
  int TrackProtonCandidate=-999;
  double TrackCandLength=0;
  double TrackProtonCandLength=0;


  trkf::TrackMomentumCalculator trkm;

  //std::cout<<"[SimpleAna] Start looping over all the tracks :  "<<std::endl;
  //loop over all the tracks and select muon candidate and proton candidate 
  for(auto track : tracks){ 
    
    TVector3 trackPos=track->Vertex();
    TVector3 trackEnd=track->End();

    double trackToVertexDist = (trackPos - nuvertexPos).Mag();
    if((trackEnd-nuvertexPos).Mag()<trackToVertexDist){
      trackPos = track->End();
      trackEnd = track->Vertex();
      trackToVertexDist=(trackPos-nuvertexPos).Mag();
    }

    std::cout<<"[SimpleAna] distance between vertex and track start is "<<trackToVertexDist<<std::endl;
    std::cout<<"[SimpleAna] distance between vertex and track end is "<<(trackEnd-nuvertexPos).Mag()<<std::endl;
    if(trackToVertexDist<fMinTrk2VtxDist){

       NumTracksNearVertex=NumTracksNearVertex+1;
       //get track length, momentum, angle and truth information of each tracks
  

       trackId->push_back(track->ID());
       tracklength->push_back(track->Length());
       trackmom->push_back(track->VertexMomentum());
       trackcostheta->push_back(TMath::Cos(track->Theta()));
       trackphi->push_back(track->Phi());
       trackstartx->push_back(trackPos.x());
       trackstarty->push_back(trackPos.y());
       trackstartz->push_back(trackPos.z());
       trackendx->push_back(trackEnd.x());
       trackendy->push_back(trackEnd.y());
       trackendz->push_back(trackEnd.z());





       //get dedx and dqdx and trancated mean dQdx for each track associated to the vertex


       std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(track.key());

       tracktrunmeandqdx->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos));
       tracktrunmeandqdx_U->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos, 0));
       tracktrunmeandqdx_V->push_back(UBXSecHelper::GetDqDxTruncatedMean(calos, 1));

       trackdEdx->push_back(GetdEdx(calos));
       trackRR->push_back(GetRR(calos)); 
       //check if this is muon like or proton like event using MIPConsistency
       //if muon like upflag=1, if proton like upflag=0;
       _muon_finder.Reset(); 
       upflag->push_back(_muon_finder.MIPConsistency(UBXSecHelper::GetDqDxTruncatedMean(calos), track->Length()));


       //std::vector<const recob::Hit*> allKHits=fmh.at(track.key());   
       std::vector<art::Ptr<recob::Hit>> allKHits = fmht.at(track.key());
       tracknhits->push_back(allKHits.size());

       //std::cout<<"[SimpleAna] Start getting the MC truth of track "<<track.key()<<"Total number of hits "<<allKHits.size() <<std::endl;
       //try to implement backtracker from Tingjun
       /*std::map<int, double> trk_ide; 
       for(size_t j=0; j<allKHits.size(); ++j){
         art::Ptr<recob::Hit> hit=allKHits[j];
         auto particles=fmhitmc.at(hit.key());
         auto hitmatch=fmhitmc.data(hit.key());
         std::cout<<"particles size = "<<particles.size()<<std::endl;
         for(size_t e=0; e<particles.size(); ++e){
            if(!particles[e]) continue;
            if(!bt->TrackIdToMotherParticle_P(particles[e]->TrackId())) continue;
              size_t trkid=(bt->TrackIdToMotherParticle_P(particles[e]->TrackId()))->TrackId();
            trk_ide[trkId] +=hitmatch[e]->energy;
         }  
       }
       double maxke=-1.;
       double totke=0;
       int Track_id=0;
                          
       for(std::map<int,double>::iterator ii=trk_mu_ide.begin();ii!=trk_mu_ide.end(); ++ii){
           totke += ii->second;
           if((ii->second)>maxke){
                 maxke = ii->second;
                 Track_mu_id=ii->first;
           }
       }

      
       const simb::MCParticle* mparticle=bt->TrackIdToParticle_P(Track_id);
       */
       //=============================================================================================
       if(isMC){
       double tmpEfracm=0;
       double tmpEcompletm=0;
       const simb::MCParticle *mparticle;
       truthMatcher(hitlist, allKHits, mparticle, tmpEfracm, tmpEcompletm);        
       std::cout<<"[SimpleAna] Check if the mparticle is valid or not: "<<mparticle<<std::endl;
       
        
       if(mparticle) 
       {
                       
       const art::Ptr<simb::MCTruth> MCtruth = fMCTruthMatching->ParticleToMCTruth(mparticle);
       //std::cout<<"SimpleAna] MCTruth origin is "<<MCtruth->Origin()<<" MCTruth nuset is "<<MCtruth->NeutrinoSet()<<std::endl; 
       
        
                         trackcand_origin->push_back(MCtruth->Origin());
                         trackcand_nuset->push_back(MCtruth->NeutrinoSet()); 
                         

                         trackcand_parPDG->push_back(mparticle->PdgCode());
                         trackcand_parStatusCode->push_back(mparticle->StatusCode()); 
                         trackcand_parTheta->push_back(mparticle->Momentum().Theta()); 
                         trackcand_parCosTheta->push_back(TMath::Cos(mparticle->Momentum().Theta()));
                         trackcand_parSinTheta->push_back(TMath::Sin(mparticle->Momentum().Theta()));
                         trackcand_parE->push_back(mparticle->E());
                         trackcand_parMass->push_back(mparticle->Mass());
                         trackcand_parKE->push_back((mparticle->E())-(mparticle->Mass()));
                         trackcand_parEndE->push_back(mparticle->EndE());
                         trackcand_parPx->push_back(mparticle->Px());
                         trackcand_parPy->push_back(mparticle->Py());
                         trackcand_parPz->push_back(mparticle->Pz());
                         trackcand_parPhi->push_back(mparticle->Momentum().Phi());
                         trackcand_parCosPhi->push_back(TMath::Cos(mparticle->Momentum().Phi()));
                         trackcand_parSinPhi->push_back(TMath::Cos(mparticle->Momentum().Phi()));
       
       std::cout<<"The PDG Code of the particle is "<<mparticle->PdgCode()<<std::endl;
       }

       }//end of if isMC for the track MCtruth tracking 

       //----------------------------------------------------------------------------------
       //get the track id and use truth matcher to get the truth information of the track
       if(track->Length()>TrackCandLength) {
           TrackCandLength=track->Length();
           TrackCandidate=track->ID();
           trackmom_mucand=trkm.GetTrackMomentum(TrackCandLength, 13);
           if(!isContained(trackPos.x(), trackPos.y(), trackPos.z()) || !isContained(trackEnd.x(), trackEnd.y(), trackEnd.z())){
               trackmom_mucand=mcsfitlist[track->ID()]->bestMomentum();
           }
       }




    }//end of is the track and vertex distance less than fMinTrk2VtxDist

  } //end of loop over all the tracks
  


  //loop over all the tracks within 5 cm from neutrino vertex 
  //make sure all the proton candidates are within the FV
  //also find out the leading proton cand length and ID
  for(auto track : tracks){ 
     TVector3 trackPos=track->Vertex();
     TVector3 trackEnd=track->End();

     double trackToVertexDist = (trackPos - nuvertexPos).Mag();

     if((trackEnd-nuvertexPos).Mag()<trackToVertexDist){
      trackPos = track->End();
      trackEnd = track->Vertex();
      trackToVertexDist=(trackPos-nuvertexPos).Mag();
     }
 
     if(trackToVertexDist<fMinTrk2VtxDist){
       if(track->ID()!=TrackCandidate) continue;
        if(!isContained(track->Vertex().x(), track->Vertex().y(),track->Vertex().z())) continue;
        if(!isContained(track->End().x(), track->End().y(),track->End().z())) continue;
        if(track->Length()>TrackProtonCandLength){
          TrackProtonCandidate=track->ID();
          TrackProtonCandLength=track->Length();
        }
     } //end of if the tracks is within 5cm from the vertex
  } //end of loop over tracks to select proton candidate
  

  std::cout<<"Track Proton Candidate is "<<TrackProtonCandidate<<std::endl;

  _ntrks=NumTracksNearVertex;


  std::cout<<"number of tracks is : "<<_ntrks<<"  number of pfparticles is : "<<_npfps<<std::endl;


  _ntrk2flag=false;
  if(_ntrks>=2) {_ntrk2flag=true;}
  
  _pdqdxflag=true;
  _upinFVflag=true;

  //#3 up in nFV cut only for proton
  for(unsigned int ncand=0; ncand<trackId->size(); ncand++){
     if(trackId->at(ncand)==TrackCandidate) continue;
     trackmom_pcand->push_back(trkm.GetTrackMomentum(tracklength->at(ncand), 2212));
     if((!isContained(trackstartx->at(ncand), trackstarty->at(ncand),trackstartz->at(ncand)))||
              (!isContained(trackendx->at(ncand), trackendy->at(ncand),trackendz->at(ncand)))) 
     {
       _upinFVflag=false;
     }     
     if(upflag->at(ncand)==1) { _pdqdxflag=false;}
  } //end of loop over ncand
  
  
  //===============================================================================================

  //
  // Save to tree
  //




  //if (right_flavour && in_tpcactive && energy_range) 
  {
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
    _bnb_correction = bnb_weight;
    _tree->Fill();
    _cc1unptree->Fill();
  }
}

DEFINE_ART_MODULE(SimpleAna)
