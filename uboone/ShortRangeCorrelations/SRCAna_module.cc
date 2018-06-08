////////////////////////////////////////////////////////////////////////
//
// @file SRCAna_module.cc
//
// @brief Short Range Correlations
//
// @authors Lu Ren (renlu23@nmsu.edu)
//
///////////////////////////////////////////////////////////////////////
#ifndef  SRCAna_Module
#define  SRCAna_Module

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
//include header files for new backtracker test
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>
#include "TTree.h"
#include "TTimeStamp.h"

constexpr double pi = 3.14159265;
constexpr float muon_mass = 0.105658;
constexpr float proton_mass =0.938272;
constexpr int kmaxg4par = 90000;
constexpr int kmaxvertex = 9999;
constexpr int kmaxtrack = 9999;
constexpr int kmaxhit = 9999;

namespace  SRCAna
{
  using namespace std;
  class  SRCAna : public art::EDAnalyzer
  {
  public:
    
    // Standard constructor and destructor for an ART module.
    explicit  SRCAna(fhicl::ParameterSet const& pset);
    virtual ~ SRCAna();
    void beginJob();
    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    bool inFV(double x, double y, double z) const;
    void truthMatcher(vector<art::Ptr<recob::Hit>> all_hits, vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    int Hammer(int nmuons, int nelectrons, int npions, int npi0, int nprotons);
    double GetDistTracktoVtx(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos);
    TVector3 GetStartTrack(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos);
    double GetDqDxTruncatedMean(vector<art::Ptr<anab::Calorimetry>> calos);
    double GetDqDxTruncatedMean_0(vector<art::Ptr<anab::Calorimetry>> calos);
    double GetDqDxTruncatedMean_1(vector<art::Ptr<anab::Calorimetry>> calos);
    double GetDqDxTruncatedMean_alt(vector<art::Ptr<anab::Calorimetry>> calos, int opt);
    bool MIPConsistency(double dqds, double length);
    double GetEnergy(vector<art::Ptr<anab::Calorimetry>> calos);
    vector<double> GetdEdx(vector<art::Ptr<anab::Calorimetry>> calos);
    vector<double> GetRR(vector<art::Ptr<anab::Calorimetry>> calos);
    double GetFlashTrackDist(double flash, double start, double end) const;
    double GetMean(vector<double> dqdx_v);
    double GetMedian(vector<double> dqdx_v);
    double GetVariance(vector<double> dqdx_v);
    double GetSTD(vector<double> dqdx_v);
    void analyze (const art::Event& evt);
    virtual void produces(art::EDProducer*);
    void endSubRun(const art::SubRun &sr); 
    void ClearLocalData();
    TTree*    Tr_ALL_Truth;
    TTree*    Tr_SRC_Truth;
    TTree*    Tr_pottree;
    TTree*    Tr_ALL_Reco;
    TTree*    Tr_SRC_Reco;
    
  private:
    unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
    string _potsum_producer; 
    string _potsum_instance;
    bool isMC;
    bool isData;
    int sr_run;
    int sr_subrun;
    double sr_begintime;
    double sr_endtime;
    double sr_pot;
    //------------------------------------------------------
  
    int Nevt_truth;   
    vector<int> *fhg4parpdg;
    vector<int> *fhg4parstatus;
    vector<float> *fhg4parpx;
    vector<float> *fhg4parpy;
    vector<float> *fhg4parpz;
    vector<float> *fhg4partheta;
    vector<float> *fhg4parphi;
    vector<float> *fhg4parp;
    float vectex[kmaxg4par];
    //===================================================
    int nGEANTparticles;
    vector<int> NuVertexID;
    vector<float> NuVertexx;
    vector<float> NuVertexy;
    vector<float> NuVertexz;
    vector<string> fVertexModuleLabelVec;
    vector<string> fPandoraNuVertexModuleLabelVec;
    vector<string> fVtxTrackAssnsModuleLabelVec;
    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int truth_double_proton; // true 2 proton
    bool truth_vtx_inFV;
    double trueMuonTrueMomentum;
    double trueMuonTrueTheta;
    double trueMuonTruePhi;
    vector<double> *trueProtonsTrueMomentum;
    vector<double> *trueProtonsTrueTheta;
    vector<double> *trueProtonsTruePhi;
    vector<double> *trueProtonsEndMomentum;

 
    TLorentzVector *fHitNucP4;
   
    //----------------------------------------
    int fNRecoTrks=-999;  //total number of tracks including muon, proton and others
    int fNRecoPTrks=-999; //total number of reco proton tracks
    int fNTruePTrks=-999; //total number of true proton tracks   
    float TrunMean_cdd;
    float TrunMean_p_cdd;
    // Other variables that will be shared between different methods.
    detinfo::DetectorProperties const* fDetectorProperties; ///< Pointer to the detector properties
    //declare the fhicl parameters here
    string              fHitsModuleLabel ;    
    string              fTrackModuleLabel;    
    string              fVertexModuleLabel;   
    string              fPandoraNuVertexModuleLabel;
    string              fGenieGenModuleLabel;   
    string              fG4ModuleLabel; 
    string              fOpFlashModuleLabel;    
    string              fCalorimetryModuleLabel; 
    string              fShowerModuleLabel; 
    string              fTrackMCSFitLabel;
    string              fParticleIDModuleLabel;
    float fG4minE;
    double fFlashWidth;      
    double fBeamMin;    
    double fBeamMax;      
    double fdQdx_scale;      
    double fPEThresh;     
    double fMinTrk2VtxDist;     
    double fMinTrackLen;
    double fDistToEdgeX;
    double fDistToEdgeY;
    double fDistToEdgeZ;
    art::EDProducer*           fMyProducerModule;    
    vector<double> _svm_x;

   
    /*
      test variables
     */
    float vtx_truth_reco_lu=-9999.0;
    float vtx_truth_reco_libo=-9999.0;
    float vtx_truth_reco_rac=-9999.0;

    /*
     Reco variables
     */
    float flash_PE = -9999.0;
    float flash_Time = -9999.0;
    bool  flash_Tag;
    bool trackflash_Tag=false;
    int   Reco_vtx_inFV=-1;
    int   Reco_trk_inFV[kmaxtrack];
   

    float Reco_Q2 = -9999.0;
    float Reco_W = -9999.0;
    float Reco_x = -9999.0;
    float Reco_y = -9999.0;
    
    float Reco_E_neutrino = -9999.0;
    float Reco_px_neutrino = -9999.0;
    float Reco_py_neutrino = -9999.0;
    float Reco_pz_neutrino = -9999.0;
    
    int   nvertex;
    float Reco_vtxID_neutrino = -9999.0;
    float Reco_vtxx_neutrino = -9999.0;
    float Reco_vtxy_neutrino = -9999.0;
    float Reco_vtxz_neutrino = -9999.0;

    int   ntrack;
    int   nhit;
    int   Reco_flip_track[kmaxtrack];
    float Reco_px_track[kmaxtrack];
    float Reco_py_track[kmaxtrack];
    float Reco_pz_track[kmaxtrack];
    float Reco_E_track[kmaxtrack];
    float Reco_startx_track[kmaxtrack];
    float Reco_starty_track[kmaxtrack];
    float Reco_startz_track[kmaxtrack];
    float Reco_start_dcosx_track[kmaxtrack];
    float Reco_start_dcosy_track[kmaxtrack];
    float Reco_start_dcosz_track[kmaxtrack];
    float Reco_p_track[kmaxtrack];
    int Reco_nhits_track[kmaxtrack];
    float Reco_TrunMean_dQdx_track[kmaxtrack];
    float Reco_TrunMean_dQdx_scaled_track[kmaxtrack];
    float Reco_TrunMean_dQdx_1st_half_track[kmaxtrack];
    float Reco_TrunMean_dQdx_2nd_half_track[kmaxtrack];
    float Reco_Mean_dQdx_track[kmaxtrack];
    float Reco_Mean_dQdx_1st_half_track[kmaxtrack];
    float Reco_Mean_dQdx_2nd_half_track[kmaxtrack];
    float Reco_range_track[kmaxtrack];
    float Reco_length_track[kmaxtrack];
    bool  Reco_TMD_muon_tag_track[kmaxtrack];
    bool  Reco_PIDA_muon_tag_track[kmaxtrack];
    bool  Reco_PIDA_pion_tag_track[kmaxtrack];
    bool  Reco_PIDA_proton_tag_track[kmaxtrack];
    float Reco_PIDA_track[kmaxtrack];

    vector<double> Reco_dQdx_track;
    vector<double> Reco_dEdx_track;
    vector<double> Reco_ResidualRange_track;

    float Reco_endx_track[kmaxtrack];
    float Reco_endy_track[kmaxtrack];
    float Reco_endz_track[kmaxtrack];
    int   Reco_MC_PDG_track[kmaxtrack];
    int   Reco_MC_origin_track[kmaxtrack];
    float Reco_MC_p_track[kmaxtrack];
    
    // float Reco_PDG_track[kmaxtrack]; 
    float Reco_mass_track[kmaxtrack];
    float Reco_vtx_track_distance[kmaxtrack];
    float Reco_flash_track_distance=-9999;
    
    float Reco_vtx_trk_distance_muon=-9999.0;
    int   Reco_nhits_muon=-1;
    int   Reco_nhits_0_muon=-1;
    int   Reco_nhits_1_muon=-1;
    int   Reco_trk_inFV_muon= 0;
    float Reco_E_muon = -9999.0;
    float Reco_p_muon = -9999.0;
    float Reco_MC_p_muon = -9999.0;
    int   Reco_MC_PDG_muon = -9999.0;
    int   Reco_MC_origin_muon = -9999.0;
    float Reco_px_muon = -9999.0;
    float Reco_py_muon = -9999.0;
    float Reco_pz_muon = -9999.0;
    float Reco_start_vtxx_muon = -9999.0;
    float Reco_start_vtxy_muon = -9999.0;
    float Reco_start_vtxz_muon = -9999.0;
    float Reco_end_vtxx_muon = -9999.0;
    float Reco_end_vtxy_muon = -9999.0;
    float Reco_end_vtxz_muon = -9999.0;
    float Reco_length_muon = -9999.0;
    float Reco_range_muon = -9999.0;
    float Reco_theta_muon = -9999.0;
    float Reco_phi_muon = -9999.0;
    float Reco_p_range_muon = -9999.0;
    float Reco_p_MCS_muon = -9999.0;
    float Reco_TrunMean_dQdx_muon = -9999.0;
    float Reco_TrunMean_dQdx_scaled_muon = -9999.0;
    bool  Reco_TMD_muon_tag_muon;
    bool  Reco_PIDA_muon_tag_muon;
    bool  Reco_PIDA_pion_tag_muon;
    bool  Reco_PIDA_proton_tag_muon;
    float Reco_PIDA_muon=-9999.0;
    float Reco_PID_chi2_muon=-9999.0;
    float Reco_PID_chi2_proton_muon=-9999.0;
    float Reco_PID_chi2_pion_muon=-9999.0;
    float Reco_PID_chi2_kaon_muon=-9999.0;
    float Reco_PID_chi2_muon_muon=-9999.0;
     
    vector<double> Reco_dQdx_muon;
    vector<double> Reco_dEdx_muon;
    vector<double> Reco_ResidualRange_muon;
   vector<double> Reco_dQdx_1_muon;
    vector<double> Reco_dEdx_1_muon;
    vector<double> Reco_ResidualRange_1_muon;
   vector<double> Reco_dQdx_0_muon;
    vector<double> Reco_dEdx_0_muon;
    vector<double> Reco_ResidualRange_0_muon;
  
 /* information from best plane */
    int   Reco_best_plane_muon=-1;
    int   Reco_best_nhits_muon=-1;
    float Reco_best_TrunMean_dQdx_muon = -9999.0;
    float Reco_best_TrunMean_dQdx_scaled_muon = -9999.0;
    bool  Reco_best_TMD_muon_tag_muon;
    float Reco_best_PIDA_muon=-9999.0;
    float Reco_best_PID_chi2_proton_muon=-9999.0;
    vector<double> Reco_best_dQdx_muon;
    vector<double> Reco_best_dEdx_muon;
    vector<double> Reco_best_ResidualRange_muon;

    int   Reco_best_plane_LP=-1;
    int   Reco_best_nhits_LP=-1;
    float Reco_best_TrunMean_dQdx_LP = -9999.0;
    float Reco_best_TrunMean_dQdx_scaled_LP = -9999.0;
    bool  Reco_best_TMD_muon_tag_LP;
    float Reco_best_PIDA_LP=-9999.0;
    float Reco_best_PID_chi2_proton_LP=-9999.0;
    vector<double> Reco_best_dQdx_LP;
    vector<double> Reco_best_dEdx_LP;
    vector<double> Reco_best_ResidualRange_LP;

    int   Reco_best_plane_SP=-1;
    int   Reco_best_nhits_SP=-1;
    float Reco_best_TrunMean_dQdx_SP = -9999.0;
    float Reco_best_TrunMean_dQdx_scaled_SP = -9999.0;
    bool  Reco_best_TMD_muon_tag_SP;
    float Reco_best_PIDA_SP=-9999.0;
    float Reco_best_PID_chi2_proton_SP=-9999.0;
    vector<double> Reco_best_dQdx_SP;
    vector<double> Reco_best_dEdx_SP;
    vector<double> Reco_best_ResidualRange_SP;
     /* end of information from best plane */
    
    float Reco_vtx_trk_distance_LP=-9999.0;
    int   Reco_nhits_LP=-1;
    int   Reco_nhits_0_LP=-1;
    int   Reco_nhits_1_LP=-1;
    int   Reco_trk_inFV_LP= 0;
    float Reco_E_LP = -9999.0;
    float Reco_p_LP = -9999.0;
    float Reco_MC_p_LP = -9999.0;
    int   Reco_MC_PDG_LP = -9999.0;
    int   Reco_MC_origin_LP = -9999.0;
    float Reco_px_LP = -9999.0;
    float Reco_py_LP = -9999.0;
    float Reco_pz_LP = -9999.0;
    float Reco_start_vtxx_LP = -9999.0;
    float Reco_start_vtxy_LP = -9999.0;
    float Reco_start_vtxz_LP = -9999.0;
    float Reco_end_vtxx_LP = -9999.0;
    float Reco_end_vtxy_LP = -9999.0;
    float Reco_end_vtxz_LP = -9999.0;
    float Reco_length_LP = -9999.0;
    float Reco_theta_LP = -9999.0;
    float Reco_phi_LP = -9999.0;
    float Reco_range_LP = -9999.0;
    float Reco_p_range_LP = -9999.0;
    float Reco_p_MCS_LP = -9999.0;
    float Reco_TrunMean_dQdx_LP = -9999.0;
    float Reco_TrunMean_dQdx_scaled_LP = -9999.0;
    bool  Reco_TMD_muon_tag_LP;
    bool  Reco_PIDA_muon_tag_LP;
    bool  Reco_PIDA_pion_tag_LP;
    bool  Reco_PIDA_proton_tag_LP;
    float Reco_PIDA_LP=-9999.0;
     float Reco_PID_chi2_LP=-9999.0;
    float Reco_PID_chi2_proton_LP=-9999.0;
    float Reco_PID_chi2_pion_LP=-9999.0;
    float Reco_PID_chi2_kaon_LP=-9999.0;
    float Reco_PID_chi2_muon_LP=-9999.0;

    vector<double> Reco_dQdx_LP;
    vector<double> Reco_dEdx_LP;
    vector<double> Reco_ResidualRange_LP;
    vector<double> Reco_dQdx_0_LP;
    vector<double> Reco_dEdx_0_LP;
    vector<double> Reco_ResidualRange_0_LP;
    vector<double> Reco_dQdx_1_LP;
    vector<double> Reco_dEdx_1_LP;
    vector<double> Reco_ResidualRange_1_LP;

    float Reco_vtx_trk_distance_SP=-9999.0;
    int   Reco_nhits_0_SP=-1;
    int   Reco_nhits_1_SP=-1;
    int   Reco_nhits_SP=-1;
    int   Reco_trk_inFV_SP= 0;
    float Reco_E_SP = -9999.0;
    float Reco_p_SP = -9999.0;
    float Reco_MC_p_SP = -9999.0;
    int   Reco_MC_PDG_SP = -9999.0;
    int   Reco_MC_origin_SP = -9999.0;
    float Reco_px_SP = -9999.0;
    float Reco_py_SP = -9999.0;
    float Reco_pz_SP = -9999.0;
    float Reco_start_vtxx_SP = -9999.0;
    float Reco_start_vtxy_SP = -9999.0;
    float Reco_start_vtxz_SP = -9999.0;
    float Reco_end_vtxx_SP = -9999.0;
    float Reco_end_vtxy_SP = -9999.0;
    float Reco_end_vtxz_SP = -9999.0;
    float Reco_length_SP = -9999.0;
    float Reco_theta_SP = -9999.0;
    float Reco_phi_SP = -9999.0;
    float Reco_range_SP = -9999.0;
    float Reco_p_range_SP = -9999.0;
    float Reco_p_MCS_SP = -9999.0;
    float Reco_TrunMean_dQdx_SP = -9999.0;
    float Reco_TrunMean_dQdx_scaled_SP = -9999.0;
    bool  Reco_TMD_muon_tag_SP;
    bool  Reco_PIDA_muon_tag_SP;
    bool  Reco_PIDA_pion_tag_SP;
    bool  Reco_PIDA_proton_tag_SP;
    float Reco_PIDA_SP=-9999.0;
    float Reco_PID_chi2_SP=-9999.0;
    float Reco_PID_chi2_proton_SP=-9999.0;
    float Reco_PID_chi2_pion_SP=-9999.0;
    float Reco_PID_chi2_kaon_SP=-9999.0;
    float Reco_PID_chi2_muon_SP=-9999.0;
    vector<double> Reco_dQdx_SP;
    vector<double> Reco_dEdx_SP;
    vector<double> Reco_ResidualRange_SP;
  vector<double> Reco_dQdx_0_SP;
    vector<double> Reco_dEdx_0_SP;
    vector<double> Reco_ResidualRange_0_SP;
  vector<double> Reco_dQdx_1_SP;
    vector<double> Reco_dEdx_1_SP;
    vector<double> Reco_ResidualRange_1_SP;


   
    /*
      Truth common
     */
    int   Truth_CCNC = -999;
    int   Truth_channel = -999;
    int   Truth_inttype = -999;
    int   Truth_nupdg = -999;
    int   Truth_N_muon = -999;
    int   Truth_N_pion = -999;
    int   Truth_N_electron = -999;
    int   Truth_N_neutron = -999;
    int   Truth_N_proton = -999; 
    int   Truth_vtx_inFV = -999;
    
    float Truth_Q2 = -9999.0;
    float Truth_W = -9999.0;
    float Truth_x = -9999.0;
    float Truth_y = -9999.0;
    
    float Truth_E_neutrino = -9999.0;
    float Truth_px_neutrino = -9999.0;
    float Truth_py_neutrino = -9999.0;
    float Truth_pz_neutrino = -9999.0;
    float Truth_vtxx_neutrino = -9999.0;
    float Truth_vtxy_neutrino = -9999.0;
    float Truth_vtxz_neutrino = -9999.0;

    float Truth_E_muon = -9999.0;
    float Truth_KE_muon = -9999.0;
    float Truth_p_muon = -9999.0;
    float Truth_px_muon = -9999.0;
    float Truth_py_muon = -9999.0;
    float Truth_pz_muon = -9999.0;
    float Truth_start_vtxx_muon = -9999.0;
    float Truth_start_vtxy_muon = -9999.0;
    float Truth_start_vtxz_muon = -9999.0;
    float Truth_end_vtxx_muon = -9999.0;
    float Truth_end_vtxy_muon = -9999.0;
    float Truth_end_vtxz_muon = -9999.0;

     /*
      Truth SRC only
     */

    float Truth_E_proton_leading = -9999.0;
    float Truth_E_proton_recoil = -9999.0;
    float Truth_KE_proton_leading = -9999.0;
    float Truth_KE_proton_recoil = -9999.0;
    float Truth_px_proton_leading = -9999.0;
    float Truth_px_proton_recoil = -9999.0;
    float Truth_py_proton_leading = -9999.0;
    float Truth_py_proton_recoil = -9999.0;
    float Truth_pz_proton_leading = -9999.0;
    float Truth_pz_proton_recoil = -9999.0;
    float Truth_start_vtxx_proton_leading = -9999.0;
    float Truth_start_vtxx_proton_recoil = -9999.0;
    float Truth_start_vtxy_proton_leading = -9999.0;
    float Truth_start_vtxy_proton_recoil = -9999.0;
    float Truth_start_vtxz_proton_leading = -9999.0;
    float Truth_start_vtxz_proton_recoil = -9999.0;
    float Truth_end_vtxx_proton_leading = -9999.0;
    float Truth_end_vtxx_proton_recoil = -9999.0;
    float Truth_end_vtxy_proton_leading = -9999.0;
    float Truth_end_vtxy_proton_recoil = -9999.0;
    float Truth_end_vtxz_proton_leading = -9999.0;
    float Truth_end_vtxz_proton_recoil = -9999.0;
    float Truth_theta_wrt_q_proton_leading = -9999.0;
    float Truth_theta_wrt_q_proton_recoil = -9999.0;
    float Truth_theta_wrt_beam_proton_leading = -9999.0;
    float Truth_theta_wrt_beam_proton_recoil = -9999.0;
    float Truth_alpha_proton_leading = -9999.0;
    float Truth_alpha_proton_recoil = -9999.0;
    
  }; // class  SRCAna
 
  // Constructor
  SRCAna:: SRCAna(fhicl::ParameterSet const& parameterSet) : EDAnalyzer(parameterSet)
  {
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }
  // Destructor
  SRCAna::~ SRCAna()
  {
  }
  //-----------------------------------------------------------------------
  void  SRCAna::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    Tr_ALL_Truth=tfs->make<TTree>("Tr_ALL_Truth", "MC Holder"); 
    Tr_ALL_Truth->Branch("Run",&fRun,"fRun/I");
    Tr_ALL_Truth->Branch("SubRun",&fSubRun,"fSubRun/I");
    Tr_ALL_Truth->Branch("Event",&fEvent,"fEvent/I");
    Tr_ALL_Truth->Branch("Truth_CCNC",&Truth_CCNC,"Truth_CCNC/I");
    Tr_ALL_Truth->Branch("Truth_channel",&Truth_channel,"Truth_channel/I");
    Tr_ALL_Truth->Branch("Truth_inttype",&Truth_inttype,"Truth_inttype/I");
    Tr_ALL_Truth->Branch("Truth_nupdg",&Truth_nupdg,"Truth_nupdg/I");
    Tr_ALL_Truth->Branch("Truth_N_muon",&Truth_N_muon,"Truth_N_muon/I");
    Tr_ALL_Truth->Branch("Truth_N_pion",&Truth_N_pion,"Truth_N_pion/I");
    Tr_ALL_Truth->Branch("Truth_N_electron",&Truth_N_electron,"Truth_N_electron/I");
    Tr_ALL_Truth->Branch("Truth_N_neutron",&Truth_N_neutron,"Truth_N_neutron/I");
    Tr_ALL_Truth->Branch("Truth_N_proton",&Truth_N_proton,"Truth_N_proton/I");
    Tr_ALL_Truth->Branch("Truth_Q2",&Truth_Q2,"Truth_Q2/F");
    Tr_ALL_Truth->Branch("Truth_W",&Truth_W,"Truth_W/F");
    Tr_ALL_Truth->Branch("Truth_x",&Truth_x,"Truth_x/F");
    Tr_ALL_Truth->Branch("Truth_y",&Truth_y,"Truth_y/F");
    // Tr_ALL_Truth->Branch("fHitNucP4","TLorentzVector",&fHitNucP4);
    Tr_ALL_Truth->Branch("Truth_vtx_inFV", &truth_vtx_inFV, "Truth_vtx_inFV/O");
    Tr_ALL_Truth->Branch("Truth_E_muon", &Truth_E_muon , "Truth_E_muon/F");
    Tr_ALL_Truth->Branch("Truth_KE_muon", &Truth_KE_muon , "Truth_KE_muon/F");
    Tr_ALL_Truth->Branch("Truth_E_neutrino", &Truth_E_neutrino , "Truth_E_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_px_muon", &Truth_px_muon , "Truth_px_muon/F");
    Tr_ALL_Truth->Branch("Truth_px_neutrino", &Truth_px_neutrino , "Truth_px_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_py_muon", &Truth_py_muon , "Truth_py_muon/F");
    Tr_ALL_Truth->Branch("Truth_py_neutrino", &Truth_py_neutrino, "Truth_py_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_pz_muon", &Truth_pz_muon , "Truth_pz_muon/F");
    Tr_ALL_Truth->Branch("Truth_pz_neutrino", &Truth_pz_neutrino , "Truth_pz_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_vtxx_neutrino", &Truth_vtxx_neutrino , "Truth_vtxx_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_vtxy_neutrino", &Truth_vtxy_neutrino , "Truth_vtxy_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_vtxz_neutrino", &Truth_vtxz_neutrino , "Truth_vtxz_neutrino/F");
    Tr_ALL_Truth->Branch("Truth_start_vtxx_muon", &Truth_start_vtxx_muon , "Truth_start_vtxx_muon/F");
    Tr_ALL_Truth->Branch("Truth_start_vtxy_muon", &Truth_start_vtxy_muon, "Truth_start_vtxy_muon/F");
    Tr_ALL_Truth->Branch("Truth_start_vtxz_muon", &Truth_start_vtxz_muon , "Truth_start_vtxz_muon/F");
    Tr_ALL_Truth->Branch("Truth_end_vtxx_muon", &Truth_end_vtxx_muon , "Truth_end_vtxx_muon/F");
    Tr_ALL_Truth->Branch("Truth_end_vtxy_muon", &Truth_end_vtxy_muon , "Truth_end_vtxy_muon/F");
    Tr_ALL_Truth->Branch("Truth_end_vtxz_muon", &Truth_end_vtxz_muon , "Truth_end_vtxz_muon/F");
    Tr_ALL_Truth->Branch("Truth_p_muon", &Truth_p_muon , "Truth_p_muon/F");
       
    Tr_SRC_Truth=tfs->make<TTree>("Tr_SRC_Truth", "MC Holder"); 
   Tr_SRC_Truth->Branch("Run",&fRun,"fRun/I");
    Tr_SRC_Truth->Branch("SubRun",&fSubRun,"fSubRun/I");
    Tr_SRC_Truth->Branch("Event",&fEvent,"fEvent/I");
    Tr_SRC_Truth->Branch("Truth_CCNC",&Truth_CCNC,"Truth_CCNC/I");
    Tr_SRC_Truth->Branch("Truth_channel",&Truth_channel,"Truth_channel/I");
    Tr_SRC_Truth->Branch("Truth_inttype",&Truth_inttype,"Truth_inttype/I");
    Tr_SRC_Truth->Branch("Truth_nupdg",&Truth_nupdg,"Truth_nupdg/I");
    Tr_SRC_Truth->Branch("Truth_N_muon",&Truth_N_muon,"Truth_N_muon/I");
    Tr_SRC_Truth->Branch("Truth_N_pion",&Truth_N_pion,"Truth_N_pion/I");
    Tr_SRC_Truth->Branch("Truth_N_electron",&Truth_N_electron,"Truth_N_electron/I");
    Tr_SRC_Truth->Branch("Truth_N_neutron",&Truth_N_neutron,"Truth_N_neutron/I");
    Tr_SRC_Truth->Branch("Truth_N_proton",&Truth_N_proton,"Truth_N_proton/I");   
    Tr_SRC_Truth->Branch("Truth_Q2",&Truth_Q2,"Truth_Q2/F");
    Tr_SRC_Truth->Branch("Truth_W",&Truth_W,"Truth_W/F");
    Tr_SRC_Truth->Branch("Truth_x",&Truth_x,"Truth_x/F");
    Tr_SRC_Truth->Branch("Truth_y",&Truth_y,"Truth_y/F");
    // Tr_SRC_Truth->Branch("fHitNucP4","TLorentzVector",&fHitNucP4);
    Tr_SRC_Truth->Branch("Truth_vtx_inFV", &truth_vtx_inFV, "Truth_vtx_inFV/I");
    Tr_SRC_Truth->Branch("Truth_E_muon", &Truth_E_muon , "Truth_E_muon/F");
    Tr_SRC_Truth->Branch("Truth_KE_muon", &Truth_KE_muon , "Truth_KE_muon/F");
    Tr_SRC_Truth->Branch("Truth_E_neutrino", &Truth_E_neutrino , "Truth_E_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_px_muon", &Truth_px_muon , "Truth_px_muon/F");
    Tr_SRC_Truth->Branch("Truth_px_neutrino", &Truth_px_neutrino , "Truth_px_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_py_muon", &Truth_py_muon , "Truth_py_muon/F");
    Tr_SRC_Truth->Branch("Truth_p_muon", &Truth_p_muon , "Truth_p_muon/F");
    Tr_SRC_Truth->Branch("Truth_py_neutrino", &Truth_py_neutrino, "Truth_py_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_pz_muon", &Truth_pz_muon , "Truth_pz_muon/F");
    Tr_SRC_Truth->Branch("Truth_pz_neutrino", &Truth_pz_neutrino , "Truth_pz_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_vtxx_neutrino", &Truth_vtxx_neutrino , "Truth_vtxx_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_vtxy_neutrino", &Truth_vtxy_neutrino , "Truth_vtxy_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_vtxz_neutrino", &Truth_vtxz_neutrino , "Truth_vtxz_neutrino/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxx_muon", &Truth_start_vtxx_muon , "Truth_start_vtxx_muon/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxy_muon", &Truth_start_vtxy_muon, "Truth_start_vtxy_muon/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxz_muon", &Truth_start_vtxz_muon , "Truth_start_vtxz_muon/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxx_muon", &Truth_end_vtxx_muon , "Truth_end_vtxx_muon/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxy_muon", &Truth_end_vtxy_muon , "Truth_end_vtxy_muon/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxz_muon", &Truth_end_vtxz_muon , "Truth_end_vtxz_muon/F");
    Tr_SRC_Truth->Branch("Truth_KE_proton_leading", &Truth_KE_proton_leading , "Truth_KE_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_KE_proton_recoil", &Truth_KE_proton_recoil, "Truth_KE_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_E_proton_leading", &Truth_E_proton_leading , "Truth_E_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_E_proton_recoil", &Truth_E_proton_recoil, "Truth_E_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_px_proton_leading", &Truth_px_proton_leading , "Truth_px_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_px_proton_recoil", &Truth_px_proton_recoil , "Truth_px_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_py_proton_leading", &Truth_py_proton_leading , "Truth_py_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_py_proton_recoil", &Truth_py_proton_recoil , "Truth_py_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_pz_proton_leading", &Truth_pz_proton_leading , "Truth_pz_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_pz_proton_recoil", &Truth_pz_proton_recoil , "Truth_pz_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxx_proton_leading", &Truth_start_vtxx_proton_leading , "Truth_start_vtxx_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxx_proton_recoil", &Truth_start_vtxx_proton_recoil , "Truth_start_vtxx_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxy_proton_leading", &Truth_start_vtxy_proton_leading , "Truth_start_vtxy_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxy_proton_recoil", &Truth_start_vtxy_proton_recoil , "Truth_start_vtxy_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxz_proton_leading", &Truth_start_vtxz_proton_leading , "Truth_start_vtxz_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_start_vtxz_proton_recoil", &Truth_start_vtxz_proton_recoil , "Truth_start_vtxz_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxx_proton_leading", &Truth_end_vtxx_proton_leading , "Truth_end_vtxx_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxx_proton_recoil", &Truth_end_vtxx_proton_recoil , "Truth_end_vtxx_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxy_proton_leading", &Truth_end_vtxy_proton_leading, "Truth_end_vtxy_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxy_proton_recoil", &Truth_end_vtxy_proton_recoil , "Truth_end_vtxy_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxz_proton_leading", &Truth_end_vtxz_proton_leading , "Truth_end_vtxz_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_end_vtxz_proton_recoil", &Truth_end_vtxz_proton_recoil , "Truth_end_vtxz_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_theta_wrt_q_proton_leading", &Truth_theta_wrt_q_proton_leading , "Truth_theta_wrt_q_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_theta_wrt_q_proton_recoil", &Truth_theta_wrt_q_proton_recoil , "Truth_theta_wrt_q_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_theta_wrt_beam_proton_leading", &Truth_theta_wrt_beam_proton_leading , "Truth_theta_wrt_beam_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_theta_wrt_beam_proton_recoil", &Truth_theta_wrt_beam_proton_recoil , "Truth_theta_wrt_beam_proton_recoil/F");
    Tr_SRC_Truth->Branch("Truth_alpha_proton_leading", &Truth_alpha_proton_leading , "Truth_alpha_proton_leading/F");
    Tr_SRC_Truth->Branch("Truth_alpha_proton_recoil", &Truth_alpha_proton_recoil , "Truth_alpha_proton_recoil/F");
  
     //========================================================================================= 
    Tr_ALL_Reco=tfs->make<TTree>("Tr_ALL_Reco","Data Holder");  
    Tr_ALL_Reco->Branch("Run",&fRun,"fRun/I");
    Tr_ALL_Reco->Branch("SubRun",&fSubRun,"fSubRun/I");
    Tr_ALL_Reco->Branch("Event",&fEvent,"fEvent/I");
    Tr_ALL_Reco->Branch("Truth_CCNC",&Truth_CCNC,"Truth_CCNC/I");
    Tr_ALL_Reco->Branch("Truth_channel",&Truth_channel,"Truth_channel/I");
    Tr_ALL_Reco->Branch("Truth_inttype",&Truth_inttype,"Truth_inttype/I");
    Tr_ALL_Reco->Branch("Truth_nupdg",&Truth_nupdg,"Truth_nupdg/I");
    
    Tr_ALL_Reco->Branch("flash_PE", &flash_PE , "flash_PE/F");
    Tr_ALL_Reco->Branch("flash_Time", &flash_Time , "flash_Time/F");
    Tr_ALL_Reco->Branch("flash_Tag", &flash_Tag , "flash_Tag/O");
    Tr_ALL_Reco->Branch("trackflash_Tag", &trackflash_Tag , "trackflash_Tag/O");
    
    Tr_ALL_Reco->Branch("Reco_pz_neutrino", &Reco_pz_neutrino , "Reco_pz_neutrino/F");
    Tr_ALL_Reco->Branch("Reco_E_neutrino", &Reco_E_neutrino , "Reco_E_neutrino/F");
    Tr_ALL_Reco->Branch("Reco_px_neutrino", &Reco_px_neutrino , "Reco_px_neutrino/F");
    Tr_ALL_Reco->Branch("Reco_py_neutrino", &Reco_py_neutrino, "Reco_py_neutrino/F");
    
    
    Tr_ALL_Reco->Branch("nvertex", &nvertex,"nvertex/I");
    Tr_ALL_Reco->Branch("Reco_vtxID_neutrino", &Reco_vtxID_neutrino,"Reco_vtxID_neutrino/I");
    Tr_ALL_Reco->Branch("Reco_vtxx_neutrino", &Reco_vtxx_neutrino,"Reco_vtxx_neutrino/F");
    Tr_ALL_Reco->Branch("Reco_vtxy_neutrino", &Reco_vtxy_neutrino,"Reco_vtxy_neutrino/F");
    Tr_ALL_Reco->Branch("Reco_vtxz_neutrino", &Reco_vtxz_neutrino,"Reco_vtxz_neutrino/F");
    
    Tr_ALL_Reco->Branch("ntrack", &ntrack,"ntrack/I");
    Tr_ALL_Reco->Branch("Reco_vtx_inFV", &Reco_vtx_inFV, "Reco_vtx_inFV/I");
    Tr_ALL_Reco->Branch("Reco_trk_inFV", &Reco_trk_inFV, "Reco_trk_inFV[ntrack]/I");
    
    Tr_ALL_Reco->Branch("Reco_px_track", &Reco_px_track,"Reco_px_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_py_track", &Reco_py_track,"Reco_py_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_pz_track", &Reco_pz_track,"Reco_pz_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_p_track", &Reco_p_track,"Reco_p_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_E_track", &Reco_E_track,"Reco_E_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_startx_track", &Reco_startx_track,"Reco_startx_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_starty_track", &Reco_starty_track,"Reco_starty_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_startz_track", &Reco_startz_track,"Reco_startz_track[ntrack]/F");

    Tr_ALL_Reco->Branch("Reco_start_dcosx_track", &Reco_start_dcosx_track,"Reco_start_dcosx_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_start_dcosy_track", &Reco_start_dcosy_track,"Reco_start_dcosy_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_start_dcosz_track", &Reco_start_dcosz_track,"Reco_start_dcosz_track[ntrack]/F");
  

    Tr_ALL_Reco->Branch("Reco_endx_track", &Reco_endx_track,"Reco_endx_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_endy_track", &Reco_endy_track,"Reco_endy_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_endz_track", &Reco_endz_track,"Reco_endz_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_mass_track", &Reco_mass_track,"Reco_mass_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_MC_PDG_track", &Reco_MC_PDG_track,"Reco_MC_PDG_track[ntrack]/I");
    Tr_ALL_Reco->Branch("Reco_MC_origin_track", &Reco_MC_origin_track,"Reco_MC_origin_track[ntrack]/I");

    Tr_ALL_Reco->Branch("Reco_MC_p_track", &Reco_MC_p_track,"Reco_MC_p_track[ntrack]/F");

    Tr_ALL_Reco->Branch("Reco_vtx_track_distance", &Reco_vtx_track_distance,"Reco_vtx_track_distance[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_flip_track", &Reco_flip_track,"Reco_flip_track[ntrack]/I");
    
     Tr_ALL_Reco->Branch("Reco_flash_track_distance", &Reco_flash_track_distance, "Reco_flash_track_distance/F");
    Tr_ALL_Reco->Branch("Truth_p_muon", &Truth_p_muon , "Truth_p_muon/F");
    Tr_ALL_Reco->Branch("Reco_nhits_track", &Reco_nhits_track,"Reco_nhits_track[ntrack]/I");
     Tr_ALL_Reco->Branch("nhit", &Reco_nhits_track,"Reco_nhits_track/I");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_track", &Reco_TrunMean_dQdx_track,"Reco_TrunMean_dQdx_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_scaled_track", &Reco_TrunMean_dQdx_scaled_track,"Reco_TrunMean_dQdx_scaled_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_1st_half_track", &Reco_TrunMean_dQdx_1st_half_track,"Reco_TrunMean_dQdx_1st_half_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_2nd_half_track", &Reco_TrunMean_dQdx_2nd_half_track,"Reco_TrunMean_dQdx_2nd_half_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_Mean_dQdx_track", &Reco_Mean_dQdx_track,"Reco_Mean_dQdx_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_Mean_dQdx_1st_half_track", &Reco_Mean_dQdx_1st_half_track,"Reco_Mean_dQdx_1st_half_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_Mean_dQdx_2nd_half_track", &Reco_Mean_dQdx_2nd_half_track,"Reco_Mean_dQdx_2nd_half_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_range_track", &Reco_range_track,"Reco_range_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_length_track", &Reco_length_track,"Reco_length_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_TMD_muon_tag_track", &Reco_TMD_muon_tag_track,"Reco_TMD_muon_tag_track[ntrack]/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_muon_tag_track", &Reco_PIDA_muon_tag_track,"Reco_PIDA_muon_tag_track[ntrack]/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_pion_tag_track", &Reco_PIDA_pion_tag_track,"Reco_PIDA_pion_tag_track[ntrack]/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_proton_tag_track", &Reco_PIDA_proton_tag_track,"Reco_PIDA_proton_tag_track[ntrack]/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_track", &Reco_PIDA_track,"Reco_PIDA_track[ntrack]/F");
    Tr_ALL_Reco->Branch("Reco_dEdx_track", "vector<double>", &Reco_dEdx_track);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_track", &Reco_ResidualRange_track,"vector<double> Reco_ResidualRange_track[ntrack]/F");
  
    /* muon block */
    Tr_ALL_Reco->Branch("Reco_E_muon", &Reco_E_muon , "Reco_E_muon/F");
    Tr_ALL_Reco->Branch("Reco_px_muon", &Reco_px_muon , "Reco_px_muon/F");
    Tr_ALL_Reco->Branch("Reco_py_muon", &Reco_py_muon , "Reco_py_muon/F");
    Tr_ALL_Reco->Branch("Reco_pz_muon", &Reco_pz_muon , "Reco_pz_muon/F");
    Tr_ALL_Reco->Branch("Reco_trk_inFV_muon", &Reco_trk_inFV_muon, "Reco_trk_inFV_muon/I");
    Tr_ALL_Reco->Branch("Reco_nhits_muon", &Reco_nhits_muon, "Reco_nhits_muon/I");
    Tr_ALL_Reco->Branch("Reco_nhits_0_muon", &Reco_nhits_0_muon, "Reco_nhits_0_muon/I");
    Tr_ALL_Reco->Branch("Reco_nhits_1_muon", &Reco_nhits_1_muon, "Reco_nhits_1_muon/I");
    Tr_ALL_Reco->Branch("Reco_p_muon", &Reco_p_muon , "Reco_p_muon/F");
    Tr_ALL_Reco->Branch("Reco_p_range_muon", &Reco_p_range_muon , "Reco_p_range_muon/F");
    Tr_ALL_Reco->Branch("Reco_p_MCS_muon", &Reco_p_MCS_muon , "Reco_p_MCS_muon/F");
    Tr_ALL_Reco->Branch("Reco_range_muon", &Reco_range_muon , "Reco_range_muon/F");
    Tr_ALL_Reco->Branch("Reco_length_muon", &Reco_length_muon , "Reco_length_muon/F"); 
    Tr_ALL_Reco->Branch("Reco_theta_muon", &Reco_theta_muon , "Reco_theta_muon/F"); 
    Tr_ALL_Reco->Branch("Reco_phi_muon", &Reco_phi_muon , "Reco_phi_muon/F"); 
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_muon", &Reco_TrunMean_dQdx_muon , "Reco_TrunMean_dQdx_muon/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_scaled_muon", &Reco_TrunMean_dQdx_scaled_muon , "Reco_TrunMean_dQdx_scaled_muon/F");   
    Tr_ALL_Reco->Branch("Reco_MC_p_muon", &Reco_MC_p_muon , "Reco_MC_p_muon/F");
    Tr_ALL_Reco->Branch("Reco_MC_PDG_muon", &Reco_MC_PDG_muon , "Reco_MC_PDG_muon/I");
    Tr_ALL_Reco->Branch("Reco_MC_origin_muon", &Reco_MC_origin_muon , "Reco_MC_origin_muon/I");
    Tr_ALL_Reco->Branch("Reco_TMD_muon_tag_muon", &Reco_TMD_muon_tag_muon , "Reco_TMD_muon_tag_muon/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_muon_tag_muon", &Reco_PIDA_muon_tag_muon,"Reco_PIDA_muon_tag_muon/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_pion_tag_muon", &Reco_PIDA_pion_tag_muon,"Reco_PIDA_pion_tag_muon/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_proton_tag_muon", &Reco_PIDA_proton_tag_muon,"Reco_PIDA_proton_tag_muon/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_muon", &Reco_PIDA_muon,"Reco_PIDA_muon/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_muon", &Reco_PID_chi2_muon,"Reco_PID_chi2_muon/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_proton_muon", &Reco_PID_chi2_proton_muon,"Reco_PID_chi2_proton_muon/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_pion_muon", &Reco_PID_chi2_pion_muon,"Reco_PID_chi2_pion_muon/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_kaon_muon", &Reco_PID_chi2_kaon_muon,"Reco_PID_chi2_kaon_muon/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_muon_muon", &Reco_PID_chi2_muon_muon,"Reco_PID_chi2_muon_muon/F");
    Tr_ALL_Reco->Branch("Reco_dQdx_muon", "vector<double>", &Reco_dQdx_muon);
    Tr_ALL_Reco->Branch("Reco_dEdx_muon", "vector<double>", &Reco_dEdx_muon);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_muon", "vector<double>", &Reco_ResidualRange_muon);
    Tr_ALL_Reco->Branch("Reco_dQdx_0_muon", "vector<double>", &Reco_dQdx_0_muon);
    Tr_ALL_Reco->Branch("Reco_dEdx_0_muon", "vector<double>", &Reco_dEdx_0_muon);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_0_muon", "vector<double>", &Reco_ResidualRange_0_muon);
    Tr_ALL_Reco->Branch("Reco_dQdx_1_muon", "vector<double>", &Reco_dQdx_1_muon);
    Tr_ALL_Reco->Branch("Reco_dEdx_1_muon", "vector<double>", &Reco_dEdx_1_muon);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_1_muon", "vector<double>", &Reco_ResidualRange_1_muon);
    Tr_ALL_Reco->Branch("Reco_start_vtxx_muon", &Reco_start_vtxx_muon , "Reco_start_vtxx_muon/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxy_muon", &Reco_start_vtxy_muon, "Reco_start_vtxy_muon/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxz_muon", &Reco_start_vtxz_muon , "Reco_start_vtxz_muon/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxx_muon", &Reco_end_vtxx_muon , "Reco_end_vtxx_muon/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxy_muon", &Reco_end_vtxy_muon , "Reco_end_vtxy_muon/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxz_muon", &Reco_end_vtxz_muon , "Reco_end_vtxz_muon/F");
    Tr_ALL_Reco->Branch("Reco_vtx_trk_distance_muon", &Reco_vtx_trk_distance_muon , "Reco_vtx_trk_distance_muon/F");
    
    Tr_ALL_Reco->Branch("Reco_best_nhits_muon", &Reco_best_nhits_muon, "Reco_best_nhits_muon/I");
    Tr_ALL_Reco->Branch("Reco_best_plane_muon", &Reco_best_plane_muon, "Reco_best_plane_muon/I");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_muon", &Reco_best_TrunMean_dQdx_muon , "Reco_best_TrunMean_dQdx_muon/F");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_scaled_muon", &Reco_best_TrunMean_dQdx_scaled_muon , "Reco_best_TrunMean_dQdx_scaled_muon/F");   
    Tr_ALL_Reco->Branch("Reco_best_PIDA_muon", &Reco_best_PIDA_muon,"Reco_best_PIDA_muon/F");
    Tr_ALL_Reco->Branch("Reco_best_TMD_muon_tag_muon", &Reco_best_TMD_muon_tag_muon , "Reco_best_TMD_muon_tag_muon/O");
    Tr_ALL_Reco->Branch("Reco_best_PID_chi2_proton_muon", &Reco_best_PID_chi2_proton_muon,"Reco_best_PID_chi2_proton_muon/F");
    Tr_ALL_Reco->Branch("Reco_best_dQdx_muon", "vector<double>", &Reco_best_dQdx_muon);
    Tr_ALL_Reco->Branch("Reco_best_dEdx_muon", "vector<double>", &Reco_best_dEdx_muon);
    Tr_ALL_Reco->Branch("Reco_best_ResidualRange_muon", "vector<double>", &Reco_best_ResidualRange_muon);
    
    Tr_ALL_Reco->Branch("Reco_best_nhits_LP", &Reco_best_nhits_LP, "Reco_best_nhits_LP/I");
    Tr_ALL_Reco->Branch("Reco_best_plane_LP", &Reco_best_plane_LP, "Reco_best_plane_LP/I");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_LP", &Reco_best_TrunMean_dQdx_LP , "Reco_best_TrunMean_dQdx_LP/F");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_scaled_LP", &Reco_best_TrunMean_dQdx_scaled_LP , "Reco_best_TrunMean_dQdx_scaled_LP/F");   
    Tr_ALL_Reco->Branch("Reco_best_PIDA_LP", &Reco_best_PIDA_LP,"Reco_best_PIDA_LP/F");
    Tr_ALL_Reco->Branch("Reco_best_TMD_muon_tag_LP", &Reco_best_TMD_muon_tag_LP , "Reco_best_TMD_muon_tag_LP/O");
    Tr_ALL_Reco->Branch("Reco_best_PID_chi2_proton_LP", &Reco_best_PID_chi2_proton_LP,"Reco_best_PID_chi2_proton_LP/F");
    Tr_ALL_Reco->Branch("Reco_best_dQdx_LP", "vector<double>", &Reco_best_dQdx_LP);
    Tr_ALL_Reco->Branch("Reco_best_dEdx_LP", "vector<double>", &Reco_best_dEdx_LP);
    Tr_ALL_Reco->Branch("Reco_best_ResidualRange_LP", "vector<double>", &Reco_best_ResidualRange_LP);

    Tr_ALL_Reco->Branch("Reco_best_nhits_SP", &Reco_best_nhits_SP, "Reco_best_nhits_SP/I");
    Tr_ALL_Reco->Branch("Reco_best_plane_SP", &Reco_best_plane_SP, "Reco_best_plane_SP/I");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_SP", &Reco_best_TrunMean_dQdx_SP , "Reco_best_TrunMean_dQdx_SP/F");
    Tr_ALL_Reco->Branch("Reco_best_TrunMean_dQdx_scaled_SP", &Reco_best_TrunMean_dQdx_scaled_SP , "Reco_best_TrunMean_dQdx_scaled_SP/F");   
    Tr_ALL_Reco->Branch("Reco_best_PIDA_SP", &Reco_best_PIDA_SP,"Reco_best_PIDA_SP/F");
    Tr_ALL_Reco->Branch("Reco_best_TMD_muon_tag_SP", &Reco_best_TMD_muon_tag_SP , "Reco_best_TMD_muon_tag_SP/O");
    Tr_ALL_Reco->Branch("Reco_best_PID_chi2_proton_SP", &Reco_best_PID_chi2_proton_SP,"Reco_best_PID_chi2_proton_SP/F");
    Tr_ALL_Reco->Branch("Reco_best_dQdx_SP", "vector<double>", &Reco_best_dQdx_SP);
    Tr_ALL_Reco->Branch("Reco_best_dEdx_SP", "vector<double>", &Reco_best_dEdx_SP);
    Tr_ALL_Reco->Branch("Reco_best_ResidualRange_SP", "vector<double>", &Reco_best_ResidualRange_SP);

  
   
   /* Long Proton block */
    Tr_ALL_Reco->Branch("Reco_E_LP", &Reco_E_LP , "Reco_E_LP/F");
    Tr_ALL_Reco->Branch("Reco_px_LP", &Reco_px_LP , "Reco_px_LP/F");
    Tr_ALL_Reco->Branch("Reco_py_LP", &Reco_py_LP , "Reco_py_LP/F");
    Tr_ALL_Reco->Branch("Reco_pz_LP", &Reco_pz_LP , "Reco_pz_LP/F");
    Tr_ALL_Reco->Branch("Reco_trk_inFV_LP", &Reco_trk_inFV_LP, "Reco_trk_inFV_LP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_LP", &Reco_nhits_LP, "Reco_nhits_LP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_0_LP", &Reco_nhits_0_LP, "Reco_nhits_0_LP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_1_LP", &Reco_nhits_1_LP, "Reco_nhits_1_LP/I");
    Tr_ALL_Reco->Branch("Reco_p_LP", &Reco_p_LP , "Reco_p_LP/F");
    Tr_ALL_Reco->Branch("Reco_p_range_LP", &Reco_p_range_LP , "Reco_p_range_LP/F");
    Tr_ALL_Reco->Branch("Reco_p_MCS_LP", &Reco_p_MCS_LP , "Reco_p_MCS_LP/F");
    Tr_ALL_Reco->Branch("Reco_range_LP", &Reco_range_LP , "Reco_range_LP/F");
    Tr_ALL_Reco->Branch("Reco_length_LP", &Reco_length_LP , "Reco_length_LP/F"); 
    Tr_ALL_Reco->Branch("Reco_theta_LP", &Reco_theta_LP , "Reco_theta_LP/F"); 
    Tr_ALL_Reco->Branch("Reco_phi_LP", &Reco_phi_LP , "Reco_phi_LP/F"); 
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_LP", &Reco_TrunMean_dQdx_LP , "Reco_TrunMean_dQdx_LP/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_scaled_LP", &Reco_TrunMean_dQdx_scaled_LP , "Reco_TrunMean_dQdx_scaled_LP/F");   
    Tr_ALL_Reco->Branch("Reco_MC_p_LP", &Reco_MC_p_LP , "Reco_MC_p_LP/F");
    Tr_ALL_Reco->Branch("Reco_MC_PDG_LP", &Reco_MC_PDG_LP , "Reco_MC_PDG_LP/I");
    Tr_ALL_Reco->Branch("Reco_MC_origin_LP", &Reco_MC_origin_LP , "Reco_MC_origin_LP/I");
    Tr_ALL_Reco->Branch("Reco_TMD_muon_tag_LP", &Reco_TMD_muon_tag_LP , "Reco_TMD_muon_tag_LP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_muon_tag_LP", &Reco_PIDA_muon_tag_LP,"Reco_PIDA_muon_tag_LP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_pion_tag_LP", &Reco_PIDA_pion_tag_LP,"Reco_PIDA_pion_tag_LP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_proton_tag_LP", &Reco_PIDA_proton_tag_LP,"Reco_PIDA_proton_tag_LP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_LP", &Reco_PIDA_LP,"Reco_PIDA_LP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_LP", &Reco_PID_chi2_LP,"Reco_PID_chi2_LP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_proton_LP", &Reco_PID_chi2_proton_LP,"Reco_PID_chi2_proton_LP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_pion_LP", &Reco_PID_chi2_pion_LP,"Reco_PID_chi2_pion_LP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_kaon_LP", &Reco_PID_chi2_kaon_LP,"Reco_PID_chi2_kaon_LP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_muon_LP", &Reco_PID_chi2_muon_LP,"Reco_PID_chi2_muon_LP/F");
    Tr_ALL_Reco->Branch("Reco_dQdx_LP", "vector<double>", &Reco_dQdx_LP);
    Tr_ALL_Reco->Branch("Reco_dEdx_LP", "vector<double>", &Reco_dEdx_LP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_LP", "vector<double>", &Reco_ResidualRange_LP);
    Tr_ALL_Reco->Branch("Reco_dQdx_0_LP", "vector<double>", &Reco_dQdx_0_LP);
    Tr_ALL_Reco->Branch("Reco_dEdx_0_LP", "vector<double>", &Reco_dEdx_0_LP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_0_LP", "vector<double>", &Reco_ResidualRange_0_LP);
    Tr_ALL_Reco->Branch("Reco_dQdx_1_LP", "vector<double>", &Reco_dQdx_1_LP);
    Tr_ALL_Reco->Branch("Reco_dEdx_1_LP", "vector<double>", &Reco_dEdx_1_LP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_1_LP", "vector<double>", &Reco_ResidualRange_1_LP);
    Tr_ALL_Reco->Branch("Reco_start_vtxx_LP", &Reco_start_vtxx_LP , "Reco_start_vtxx_LP/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxy_LP", &Reco_start_vtxy_LP, "Reco_start_vtxy_LP/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxz_LP", &Reco_start_vtxz_LP , "Reco_start_vtxz_LP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxx_LP", &Reco_end_vtxx_LP , "Reco_end_vtxx_LP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxy_LP", &Reco_end_vtxy_LP , "Reco_end_vtxy_LP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxz_LP", &Reco_end_vtxz_LP , "Reco_end_vtxz_LP/F");
    Tr_ALL_Reco->Branch("Reco_vtx_trk_distance_LP", &Reco_vtx_trk_distance_LP , "Reco_vtx_trk_distance_LP/F");
    
   
   /* Short Proton block */
    Tr_ALL_Reco->Branch("Reco_E_SP", &Reco_E_SP , "Reco_E_SP/F");
    Tr_ALL_Reco->Branch("Reco_px_SP", &Reco_px_SP , "Reco_px_SP/F");
    Tr_ALL_Reco->Branch("Reco_py_SP", &Reco_py_SP , "Reco_py_SP/F");
    Tr_ALL_Reco->Branch("Reco_pz_SP", &Reco_pz_SP , "Reco_pz_SP/F");
    Tr_ALL_Reco->Branch("Reco_trk_inFV_SP", &Reco_trk_inFV_SP, "Reco_trk_inFV_SP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_SP", &Reco_nhits_SP, "Reco_nhits_SP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_0_SP", &Reco_nhits_0_SP, "Reco_nhits_0_SP/I");
    Tr_ALL_Reco->Branch("Reco_nhits_1_SP", &Reco_nhits_1_SP, "Reco_nhits_1_SP/I");
    Tr_ALL_Reco->Branch("Reco_p_SP", &Reco_p_SP , "Reco_p_SP/F");
    Tr_ALL_Reco->Branch("Reco_p_range_SP", &Reco_p_range_SP , "Reco_p_range_SP/F");
    Tr_ALL_Reco->Branch("Reco_p_MCS_SP", &Reco_p_MCS_SP , "Reco_p_MCS_SP/F");
    Tr_ALL_Reco->Branch("Reco_range_SP", &Reco_range_SP , "Reco_range_SP/F");
    Tr_ALL_Reco->Branch("Reco_length_SP", &Reco_length_SP , "Reco_length_SP/F"); 
    Tr_ALL_Reco->Branch("Reco_theta_SP", &Reco_theta_SP , "Reco_theta_SP/F"); 
    Tr_ALL_Reco->Branch("Reco_phi_SP", &Reco_phi_SP , "Reco_phi_SP/F"); 
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_SP", &Reco_TrunMean_dQdx_SP , "Reco_TrunMean_dQdx_SP/F");
    Tr_ALL_Reco->Branch("Reco_TrunMean_dQdx_scaled_SP", &Reco_TrunMean_dQdx_scaled_SP , "Reco_TrunMean_dQdx_scaled_SP/F");   
    Tr_ALL_Reco->Branch("Reco_MC_p_SP", &Reco_MC_p_SP , "Reco_MC_p_SP/F");
    Tr_ALL_Reco->Branch("Reco_MC_PDG_SP", &Reco_MC_PDG_SP , "Reco_MC_PDG_SP/I");
    Tr_ALL_Reco->Branch("Reco_MC_origin_SP", &Reco_MC_origin_SP , "Reco_MC_origin_SP/I");
    Tr_ALL_Reco->Branch("Reco_TMD_muon_tag_SP", &Reco_TMD_muon_tag_SP , "Reco_TMD_muon_tag_SP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_muon_tag_SP", &Reco_PIDA_muon_tag_SP,"Reco_PIDA_muon_tag_SP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_pion_tag_SP", &Reco_PIDA_pion_tag_SP,"Reco_PIDA_pion_tag_SP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_proton_tag_SP", &Reco_PIDA_proton_tag_SP,"Reco_PIDA_proton_tag_SP/O");
    Tr_ALL_Reco->Branch("Reco_PIDA_SP", &Reco_PIDA_SP,"Reco_PIDA_SP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_SP", &Reco_PID_chi2_SP,"Reco_PID_chi2_SP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_proton_SP", &Reco_PID_chi2_proton_SP,"Reco_PID_chi2_proton_SP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_pion_SP", &Reco_PID_chi2_pion_SP,"Reco_PID_chi2_pion_SP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_kaon_SP", &Reco_PID_chi2_kaon_SP,"Reco_PID_chi2_kaon_SP/F");
    Tr_ALL_Reco->Branch("Reco_PID_chi2_muon_SP", &Reco_PID_chi2_muon_SP,"Reco_PID_chi2_muon_SP/F");
    Tr_ALL_Reco->Branch("Reco_dQdx_SP", "vector<double>", &Reco_dQdx_SP);
    Tr_ALL_Reco->Branch("Reco_dEdx_SP", "vector<double>", &Reco_dEdx_SP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_SP", "vector<double>", &Reco_ResidualRange_SP);
    Tr_ALL_Reco->Branch("Reco_dQdx_0_SP", "vector<double>", &Reco_dQdx_0_SP);
    Tr_ALL_Reco->Branch("Reco_dEdx_0_SP", "vector<double>", &Reco_dEdx_0_SP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_0_SP", "vector<double>", &Reco_ResidualRange_0_SP);
    Tr_ALL_Reco->Branch("Reco_dQdx_1_SP", "vector<double>", &Reco_dQdx_1_SP);
    Tr_ALL_Reco->Branch("Reco_dEdx_1_SP", "vector<double>", &Reco_dEdx_1_SP);
    Tr_ALL_Reco->Branch("Reco_ResidualRange_1_SP", "vector<double>", &Reco_ResidualRange_1_SP);
    Tr_ALL_Reco->Branch("Reco_start_vtxx_SP", &Reco_start_vtxx_SP , "Reco_start_vtxx_SP/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxy_SP", &Reco_start_vtxy_SP, "Reco_start_vtxy_SP/F");
    Tr_ALL_Reco->Branch("Reco_start_vtxz_SP", &Reco_start_vtxz_SP , "Reco_start_vtxz_SP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxx_SP", &Reco_end_vtxx_SP , "Reco_end_vtxx_SP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxy_SP", &Reco_end_vtxy_SP , "Reco_end_vtxy_SP/F");
    Tr_ALL_Reco->Branch("Reco_end_vtxz_SP", &Reco_end_vtxz_SP , "Reco_end_vtxz_SP/F");
    Tr_ALL_Reco->Branch("Reco_vtx_trk_distance_SP", &Reco_vtx_trk_distance_SP , "Reco_vtx_trk_distance_SP/F");
    
 
   
  
  
    /*
      this block is fir testing only 
    */
    Tr_ALL_Reco->Branch("Truth_vtx_inFV", &truth_vtx_inFV, "Truth_vtx_inFV/O");
    Tr_ALL_Reco->Branch("Truth_vtxx_neutrino", &Truth_vtxx_neutrino , "Truth_vtxx_neutrino/F");
    Tr_ALL_Reco->Branch("Truth_vtxy_neutrino", &Truth_vtxy_neutrino , "Truth_vtxy_neutrino/F");
    Tr_ALL_Reco->Branch("Truth_vtxz_neutrino", &Truth_vtxz_neutrino , "Truth_vtxz_neutrino/F");
    Tr_ALL_Reco->Branch("vtx_truth_reco_lu", &vtx_truth_reco_lu,"vtx_truth_reco_lu/F");
    Tr_ALL_Reco->Branch("vtx_truth_reco_libo", &vtx_truth_reco_libo,"vtx_truth_reco_libo/F");
    Tr_ALL_Reco->Branch("vtx_truth_reco_rac", &vtx_truth_reco_rac,"vtx_truth_reco_rac/F");
 
    //=================================================================   
    Tr_SRC_Reco=tfs->make<TTree>("Tr_SRC_Reco","Data Holder");  
    Tr_SRC_Reco->Branch("flash_PE", &flash_PE , "flash_PE/F");
    Tr_SRC_Reco->Branch("flash_Time", &flash_Time , "flash_Time/F");
    Tr_SRC_Reco->Branch("flash_Tag", &flash_Tag , "flash_Tag/O");
   
     //=================================================================
 
    Tr_pottree =tfs->make<TTree>("Tr_pottree","");
    Tr_pottree->Branch("run",                &sr_run,                "run/I");
    Tr_pottree->Branch("subrun",             &sr_subrun,             "subrun/I");
    Tr_pottree->Branch("begintime",          &sr_begintime,          "begintime/D");
    Tr_pottree->Branch("endtime",            &sr_endtime,            "endtime/D");
    Tr_pottree->Branch("pot", &sr_pot, "pot/D");

  }
  void  SRCAna::beginRun(const art::Run& /*run*/)
  {
  }
  
  //-----------------------------------------------------------------------
  void  SRCAna::reconfigure(fhicl::ParameterSet const& pset)
  {
    fVertexModuleLabelVec        = pset.get< vector<string> >("VertexModuleLabelVec",       vector<string>() ={"vertex3d"});
    fPandoraNuVertexModuleLabelVec  = pset.get< vector<string> >("PandoraNuVertexModuleLabelVec", vector<string>()={"pandoraNu"}); 
    fVtxTrackAssnsModuleLabelVec = pset.get< vector<string> >("VtxTrackAssnModuleLabelVec", vector<string>() ={"neutrinoID"});
    if (fVertexModuleLabelVec.size() != fVtxTrackAssnsModuleLabelVec.size())
      {
        mf::LogError("TPCNeutrinoIDFilter") << "Mismatch between string vector lengths input from fhicl!" << endl;
      }
    auto const* geom = lar::providerFrom<geo::Geometry>(); // geometry is needed to go from OpChannel to OpDet 
    //fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fOpFlashModuleLabel      = pset.get<string> ("OpFlashModuleLabel", "simpleFlashBeam");
    fHitsModuleLabel         = pset.get<string> ("HitsModuleLabel", "gaushit");    
    fTrackModuleLabel        = pset.get<string> ("TrackModuleLabel",  "pandoraNu");    
    fVertexModuleLabel       = pset.get<string> ("VertexModuleLabel",  "pandoraNu");    
    fPandoraNuVertexModuleLabel=pset.get<string> ("PandoraNuVertexModuleLabel", "pandoraNu");
    fGenieGenModuleLabel     = pset.get<string> ("GenieGenModuleLabel", "generator");        
    fG4ModuleLabel           = pset.get<string> ("G4ModuleLabel", "largeant");
    fCalorimetryModuleLabel  = pset.get<string> ("CalorimetryModuleLabel", "pandoraNucalo");
    fShowerModuleLabel       = pset.get<string> ("ShowerModuleLabel", "showerrecopandora");
    fTrackMCSFitLabel        = pset.get<string> ("TrackMCSFitLabel", "pandoraNuMCSMu" );
    fParticleIDModuleLabel    = pset.get<string> ("ParticleIDModuleLabel");
    fDistToEdgeX             = geom->DetHalfWidth()   - pset.get("DistToEdgeX",   10.);    
    fDistToEdgeY             = geom->DetHalfHeight()  - pset.get("DistToEdgeY",   20.);    
    fDistToEdgeZ             = geom->DetLength() / 2. - pset.get("DistToEdgeZ",   10.);        
    fFlashWidth              = pset.get      ("FlashWidth", 80.);    
     fdQdx_scale              = pset.get("dQdxScale",196.98);
     // fBeamMin                 = pset.get      ("BeamMin", 3.65);   //extbnb offbeam data
     //fBeamMax                 = pset.get      ("BeamMax", 5.25);   //extbnb offbeam data
    // fBeamMin                 = pset.get      ("BeamMin", 3.3);   //bnb onbeam data
    //fBeamMax                 = pset.get      ("BeamMax", 4.9);   //bnb onbeam data
   fBeamMin                 = pset.get      ("BeamMin", 3.2);   //BNB+COSMIC MC
    fBeamMax                 = pset.get      ("BeamMax", 4.8);   //BNB+COSMIC MC
    fPEThresh                = pset.get      ("PEThresh", 50.);    
    fMinTrk2VtxDist          = pset.get      ("MinTrk2VtxDist", 5.);    
    fMinTrackLen             = pset.get      ("MinTrackLen", 75.);
    fG4minE                  = pset.get      ("G4minE",0.01); 
    _svm_x = pset.get<vector<double>> ("SVM_X");

    //Get the tool for MC truth Matching
    const fhicl::ParameterSet& truthParams = pset.get<fhicl::ParameterSet>("MCTruthMatching");    
    if (truthParams.get<string>("tool_type") == "AssociationsTruth")
      {
        fMCTruthMatching = unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
      }
    else
      {
        fMCTruthMatching = unique_ptr<truth::IMCTruthMatching>(new truth::BackTrackerTruth(truthParams));
      }
    _potsum_producer = pset.get<string>("POTSummaryProducer");
    _potsum_instance = pset.get<string>("POTSummaryInstance");
    
    fhg4parpdg = new vector<int>;
    fhg4parstatus= new vector<int>;
    fhg4parpx= new vector<float>;
    fhg4parpy= new vector<float>;
    fhg4parpz= new vector<float>;
    fhg4parp = new vector<float>;
    fhg4partheta= new vector<float>;
    fhg4parphi= new vector<float>;
    fHitNucP4 = new TLorentzVector(-999,-999,-999,-999);
    return;
  }
  bool SRCAna::inFV(double x, double y, double z) const
  {
    auto const* geom = lar::providerFrom<geo::Geometry>();
    double distInX = x - geom->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * geom->DetLength();
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    return false;
  }
 
vector<double> SRCAna::GetRR(vector<art::Ptr<anab::Calorimetry>> calos) {
  vector<double> result;
  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
   
    vector<double> RR_v = c->ResidualRange(); 
    return RR_v;
  }
  return result;
}
vector<double> SRCAna::GetdEdx(vector<art::Ptr<anab::Calorimetry>> calos) {
  vector<double> result;
  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
   
    vector<double> dEdx_v = c->dEdx(); 
//    if (dEdx_v.size() == 0){
//      return 0;
//    }
    return dEdx_v;
  }
  return result;
}
double SRCAna::GetEnergy(vector<art::Ptr<anab::Calorimetry>> calos) {
  double result=-9999;
  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
    cout << "old Energy =" << c->KineticEnergy() << "new result = ";  
    vector<double> dEdx_v = c->dEdx(); 
    vector<double> TrkPtch_v = c->TrkPitchVec(); 
//
    if (dEdx_v.size() == 0){
      return result;
    }
    double totE=0;
    for (unsigned int i(0); i < dEdx_v.size();++i){
      totE += dEdx_v.at(i) * TrkPtch_v.at(i);
    }
    result = totE;
    cout << totE << endl;
  }
  return result;
}
double SRCAna::GetMean(vector<double> dqdx_v) {
  double mean = -9999;
  size_t size = dqdx_v.size();
  if (size == 0)
    return mean;
  double sum = 0;
  for (auto v : dqdx_v) {
    sum += v;
  }
  mean = sum/(double)size;
  return mean;
}
double SRCAna::GetMedian(vector<double> dqdx_v) {
  double median = -9999;
  size_t size = dqdx_v.size();
  if (size == 0)
    return median;
  sort(dqdx_v.begin(), dqdx_v.end());
  if (size % 2 == 0){
    median = (dqdx_v[size/2 - 1] + dqdx_v[size/2]) / 2;
  }
  else{
    median = dqdx_v[size/2];
  }
  return median;
}
double SRCAna::GetVariance(vector<double> dqdx_v) {
  double variance = -1;
  double sum = 0;
  double sum2 = 0;
  size_t size = dqdx_v.size();
  if (size == 0)
    return variance;
  for (auto value : dqdx_v) {
    sum  += value;
    sum2 += value*value;
  }  
  variance = sum2/(double)size - (sum/(double)size)*(sum/(double)size);
  return variance;
}
double SRCAna::GetSTD(vector<double> dqdx_v) {
  if (dqdx_v.size() == 0)
    return -9999;
  double variance = GetVariance(dqdx_v);
  if (variance > 0)
    return sqrt(variance);
  else 
    return -9999;
}

double SRCAna::GetDqDxTruncatedMean(vector<art::Ptr<anab::Calorimetry>> calos) {
  double result = -9999;
  double n = 1.;
  for (auto c : calos) {
   
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
    cout<<"C? "<<c<<" , plane num: "<<planenum<<endl;
    vector<double> dqdx_v = c->dQdx(); 
    if (dqdx_v.size() == 0)
      return result;
    // for (auto q : dqdx_v) {
    //   cout << "dqdx before trim: " << q << endl;
    // }
    
    double median = GetMedian(dqdx_v);
    double std    = GetSTD(dqdx_v);
    
    //      cout << "median " << median << endl;
    // cout << "std    " << std << endl;
    
    vector<double> dqdx_v_trimmed;
    dqdx_v_trimmed.clear();
    
    for (auto q : dqdx_v) {
      if (q > median - n * std && 
          q < median + n * std) {
         dqdx_v_trimmed.emplace_back(q);
	 //	 cout << "dqdx after trim: " << q <<endl;
      }
    }
    
    result = GetMean(dqdx_v_trimmed);
    //cout<<dqdx_v.size()<<" "<<dqdx_v_trimmed.size()<<endl;
   
  }
    
  return result;
}
  /////
  double SRCAna::GetDqDxTruncatedMean_0(vector<art::Ptr<anab::Calorimetry>> calos) 
  {
    double result = -9999;
    double n = 1.;
    for (auto c : calos) 
      {   
	if (!c) continue;
	if (!c->PlaneID().isValid) continue;
	int planenum = c->PlaneID().Plane;
	if (planenum != 0) continue;
	cout<<"C? "<<c<<" , plane num: "<<planenum<<endl;
	vector<double> dqdx_v = c->dQdx(); 
	if (dqdx_v.size() == 0)
	  return result;
	double median = GetMedian(dqdx_v);
	double std    = GetSTD(dqdx_v);
	vector<double> dqdx_v_trimmed;
	dqdx_v_trimmed.clear();
	for (auto q : dqdx_v) 
	  {
	    if (q > median - n * std && q < median + n * std) 
	      {
		dqdx_v_trimmed.emplace_back(q);
	      }
	  }
	result = GetMean(dqdx_v_trimmed);
      }   
    return result;
  }
  
  ////
 double SRCAna::GetDqDxTruncatedMean_1(vector<art::Ptr<anab::Calorimetry>> calos) 
  {
    double result = -9999;
    double n = 1.;
    for (auto c : calos) 
      {   
	if (!c) continue;
	if (!c->PlaneID().isValid) continue;
	int planenum = c->PlaneID().Plane;
	if (planenum != 1) continue;
	vector<double> dqdx_v = c->dQdx(); 
	if (dqdx_v.size() == 0)
	  return result;
	double median = GetMedian(dqdx_v);
	double std    = GetSTD(dqdx_v);
	vector<double> dqdx_v_trimmed;
	dqdx_v_trimmed.clear();
	for (auto q : dqdx_v) 
	  {
	    if (q > median - n * std && q < median + n * std) 
	      {
		dqdx_v_trimmed.emplace_back(q);
	      }
	  }
	result = GetMean(dqdx_v_trimmed);
      }   
    return result;
  }

  double SRCAna::GetDqDxTruncatedMean_alt(vector<art::Ptr<anab::Calorimetry>> calos, int opt = 0) 
  {
    double result_0 = -9999;//truncated mean
    double result_1 = -9999;//truncated mean using 1st half of the track
    double result_2 = -9999;//truncated mean using 2nd half of the track
    double result_3 = -9999;//mean
    double result_4 = -9999;//mean using 1st half of the track
    double result_5 = -9999;//mean using 2nd half of the track 
    double n = 1.;
    
    vector<double> dqdx_v = calos[2]->dQdx(); 
    if (dqdx_v.size() == 0)
      return result_0;
    size_t count = dqdx_v.size();
    size_t half_count = 0;
    if(count%2==0)
      half_count = count/2;
    else
      half_count = (count+1)/2;
    vector<double>::const_iterator first=dqdx_v.begin();
    vector<double>::const_iterator half=dqdx_v.begin()+ half_count;
    vector<double>::const_iterator end=dqdx_v.begin()+ count;
    
    vector<double> dqdx_v_1(first, half);
    vector<double> dqdx_v_2(half, end);	
    
    
    if(opt ==0)
      {
	double median = GetMedian(dqdx_v);
	double std    = GetSTD(dqdx_v);
	vector<double> dqdx_v_trimmed;
	dqdx_v_trimmed.clear();
	for (size_t q=0;q<count;q++) 
	  {
	    if ( dqdx_v[q] > median - n * std && dqdx_v[q] < median + n * std) 
	      {
		dqdx_v_trimmed.emplace_back( dqdx_v[q]);
	      }
	  }
	result_0 = GetMean(dqdx_v_trimmed);
	return result_0;
      }
    else if(opt==1)
      {
	double median = GetMedian(dqdx_v_1);
	double std    = GetSTD(dqdx_v_1);
	vector<double> dqdx_v_trimmed;
	dqdx_v_trimmed.clear();
	for (size_t q=0;q<dqdx_v_1.size();q++) 
	  {
	    if ( dqdx_v_1[q] > median - n * std && dqdx_v_1[q] < median + n * std) 
	      {
		dqdx_v_trimmed.emplace_back( dqdx_v_1[q]);
	      }
	  }
	result_1 = GetMean(dqdx_v_trimmed);
	return result_1;
      }
    else if(opt==2)
      {
	double median = GetMedian(dqdx_v_2);
	double std    = GetSTD(dqdx_v_2);
	vector<double> dqdx_v_trimmed;
	dqdx_v_trimmed.clear();
	for (size_t q=0;q<dqdx_v_2.size();q++) 
	  {
	    if ( dqdx_v_2[q] > median - n * std && dqdx_v_2[q] < median + n * std) 
	      {
		dqdx_v_trimmed.emplace_back( dqdx_v_2[q]);
	      }
	  }
	result_2 = GetMean(dqdx_v_trimmed);
	return result_2;
      }
    else if(opt==3)
      {
	result_3 = GetMean(dqdx_v);
	return result_3;
      }
    else if(opt==4)
      {
	result_4 = GetMean(dqdx_v_1);
	return result_4;
      }
    else if(opt==5)
      {
	result_5 = GetMean(dqdx_v_2);
	return result_5;
      }
    else
      return -9999;
      
  }
  bool SRCAna::MIPConsistency(double dqds, double length) {
    if (length > 1000)
      return true;
    if (length < 0) {
      //     cout << "[MuonCandidateFinder] Track length < 0?!" << endl;
      return false;
    }
    int l = round(length);
    double dqds_cut = _svm_x.at(l);
    //    cout << "[MuonCandidateFinder] Track length is " << length << ", dqds_cut is " << dqds_cut << ", dqds value is " << dqds*fdQdx_scale << endl;
 
    //if (dqds*242.77 <= dqds_cut) //data
    //if (dqds*196.98 <= dqds_cut) //MC
    if (dqds * fdQdx_scale <= dqds_cut) //Use a fcl variable to change from data to MC
      return true;
  
    return false;
}

double SRCAna::GetFlashTrackDist(double flash, double start, double end) const
{
    if (end >= start) {
        if (flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if (flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}
void SRCAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet)
{       
  std::map<int,double> trkID_E;        
  for(size_t j = 0; j < track_hits.size(); ++j)
    {        
      art::Ptr<recob::Hit> hit = track_hits[j];
      std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
      for(size_t k = 0; k < TrackIDs.size(); k++)
	{        
	  trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
	}
    }
  double E_em =0.0;        
  double max_E = -999.0;        
  double total_E = 0.0;        
  int TrackID = -999;        
  double partial_E =0.0;          
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
  MCparticle = fMCTruthMatching->TrackIDToParticle(TrackID);                      
  if( TrackID < 0 ) return;                      
  Efrac = (partial_E)/total_E;                
  //completeness        
  double totenergy =0;        
  for(size_t k = 0; k < all_hits.size(); ++k)
  {          
    art::Ptr<recob::Hit> hit = all_hits[k];              
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);        
    for(size_t l = 0; l < TrackIDs.size(); ++l)
      {        
	if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;        
      }        
  }         
  Ecomplet = partial_E/totenergy;
}

  ///
void SRCAna::endSubRun(const art::SubRun& sr){
  // if (_debug) cout << "[CC1uNPSelAna::endSubRun] Starts" << endl;
  // Saving run and subrun number on file so that we can run Zarko's script easily
  //_run_subrun_list_file << sr.run() << " " << sr.subRun() << endl;
  
  sr_run       = sr.run();
  sr_subrun    = sr.subRun();
  sr_begintime = sr.beginTime().value();
  sr_endtime   = sr.endTime().value();
  
  art::Handle<sumdata::POTSummary> potsum_h;
  
  // MC
  // if (isMC) {
  //    if (_debug) cout << "CC1uNPSelAna::endSubRun] Getting POT for MC" << endl;
  //    if(sr.getByLabel(_potsum_producer, potsum_h)) {
  //       if (_debug) cout << "CC1uNPSelAna POT are valid, POT = " << potsum_h->totpot  << endl;
  //       sr_pot = (double)potsum_h->totpot/1.0e16;
  //    }
  //    else
  //    sr_pot = 0.;
  // }
  /*
   // Data
   if (isData) {
     if (_debug) cout << "[CC1uNPSelAna::endSubRun] Getting POT for DATA, producer " << _potsum_producer << ", instance " << _potsum_instance << endl;
     if (sr.getByLabel(_potsum_producer, _potsum_instance, potsum_h)){
        if (_debug) cout << "[CC1uNPSelAna::endSubRun] POT are valid" << endl;
          _sr_pot = potsum_h->totpot;
     }
   else
   _sr_pot = 0;
  }
  */
  Tr_pottree->Fill();
  
  //  if (_debug) cout << "[CC1uNPSelAna::endSubRun] Ends" << endl;
}

void SRCAna::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();
}
void SRCAna::ClearLocalData(){
  fill(vectex,       vectex+sizeof(vectex)/sizeof(vectex[0]), -9999.);
}

  int SRCAna::Hammer(int nmuons, int nelectrons, int npions, int npi0, int nprotons)
  {
    if (nmuons ==1 && (nelectrons + npions + npi0 ) == 0 && nprotons ==2 ) return 1;
    else return -1;
  }
  double SRCAna::GetDistTracktoVtx(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos)
  { 
    TVector3 trk_pos = InputTrackPtr->Vertex();
    TVector3 trk_end = InputTrackPtr->End();
    double diststart = (trk_pos - InputvertexPos).Mag();
    double distend = (trk_pos - InputvertexPos).Mag(); 
    if(diststart> distend ) return distend;
    else return diststart;    
  }
TVector3 SRCAna::GetStartTrack(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos){

  TVector3 trackPos = InputTrackPtr->Vertex();
  TVector3 trackEnd = InputTrackPtr->End();

  // Take the closer end---------------------------------                                                                                                 
  double diststart = (trackPos - InputvertexPos).Mag();
  double distend = (trackEnd - InputvertexPos).Mag();

  if(diststart> distend ) trackPos = trackEnd;

  return trackPos;

}
 /*--------------------------------------------------
  
  Analyze...
  
  --------------------------------------------------*/
  void  SRCAna::analyze(const art::Event& event)
  {
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    isMC = !event.isRealData();
    isData= !isMC;
    Int_t nGeniePrimaries=0;
    Int_t nGeantParticles=0;
    TVector3 q_tran;
    
    
    if(isMC)
      { 
	fMCTruthMatching->Rebuild(event); 
	art::Handle< vector<simb::MCTruth> > mctruthListHandle;
	vector<art::Ptr<simb::MCTruth> > mclist;
	art::Ptr<simb::MCTruth> mctruth; 
	if (event.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
	  art::fill_ptr_vector(mclist, mctruthListHandle);      
	Nevt_truth=mclist.size();
	mctruth = mclist[0];
	if (mctruth-> NeutrinoSet()) 
	  nGeniePrimaries = mctruth->NParticles();
	if(mctruth-> NeutrinoSet() && mctruth->Origin() == simb::kBeamNeutrino) 
	  {
	    Truth_CCNC = mctruth -> GetNeutrino().CCNC();  //0 CC, 1 NC
	    Truth_channel = mctruth -> GetNeutrino().Mode();   //0 QE, 1 RES, 2 DIS, 3 COH
	    Truth_inttype = mctruth -> GetNeutrino().InteractionType();
	    Truth_nupdg = mctruth -> GetNeutrino().Nu().PdgCode();
	    Truth_Q2 = mctruth -> GetNeutrino().QSqr(); 
	    Truth_W = mctruth -> GetNeutrino().W();
	    Truth_x = mctruth -> GetNeutrino().X();
	    Truth_y = mctruth -> GetNeutrino().Y();	  

	    Truth_E_neutrino=mctruth->GetNeutrino().Nu().E(); 
	    Truth_px_neutrino = mctruth->GetNeutrino().Nu().Px(); 
	    Truth_py_neutrino = mctruth->GetNeutrino().Nu().Py(); 
	    Truth_pz_neutrino = mctruth->GetNeutrino().Nu().Pz();
	    Truth_vtxx_neutrino = mctruth->GetNeutrino().Nu().Vx(); 
	    Truth_vtxy_neutrino = mctruth->GetNeutrino().Nu().Vy();
	    Truth_vtxz_neutrino = mctruth->GetNeutrino().Nu().Vz(); 
	    
	    Truth_E_muon=mctruth->GetNeutrino().Lepton().E(); 
	    Truth_KE_muon=Truth_E_muon-muon_mass;
	    Truth_px_muon=mctruth->GetNeutrino().Lepton().Px(); 
	    Truth_py_muon=mctruth->GetNeutrino().Lepton().Py(); 
	    Truth_pz_muon=mctruth->GetNeutrino().Lepton().Pz(); 
	    Truth_p_muon=sqrt(Truth_px_muon*Truth_px_muon+Truth_py_muon*Truth_py_muon+Truth_pz_muon*Truth_pz_muon);
	    q_tran.SetXYZ(Truth_px_neutrino-Truth_px_muon,Truth_py_neutrino-Truth_py_muon,Truth_pz_neutrino-Truth_pz_muon);
	    
	    for (int igeniepart(0); igeniepart<nGeniePrimaries; igeniepart++)
	      {
		simb::MCParticle part = mctruth->GetParticle(igeniepart);
		if (part.PdgCode()==2212 && part.StatusCode()==1 && part.Mother()==0)
		  {
		    trueProtonsTrueMomentum->push_back(part.P());
		    trueProtonsTrueTheta->push_back(part.Momentum().Theta());
		    trueProtonsTruePhi->push_back(part.Momentum().Phi());
		    trueProtonsEndMomentum->push_back(part.EndMomentum().P());
		  }

	      }    
	  }
	art::Handle< vector<simb::MCParticle> > mcparticleListHandle;
	vector< art::Ptr<simb::MCParticle> > mcparticlelist;
	if (event.getByLabel(fG4ModuleLabel,mcparticleListHandle))
	  art::fill_ptr_vector(mcparticlelist, mcparticleListHandle);
	map<int, size_t> TrackIDtoIndex;
	vector<int> gpdg;   //vector to store the pdg id of geant4
	vector<int> gmother; //vector to store the mother id of geant4 
	Int_t nmuons = 0;
	Int_t npions = 0;
	Int_t npi0 = 0;
	Int_t nprotons = 0;
	Int_t nelectrons = 0;
	Int_t nneutrons = 0;
	Int_t proton_holder_1 = -999;
	Int_t proton_holder_2 = -999;
	Int_t proton_holder_more = -999;
	float e_1,px_1,py_1,pz_1,startx_1,starty_1,startz_1,endx_1,endy_1,endz_1,alpha_1,theta_beam_1,theta_q_1;
	float e_2,px_2,py_2,pz_2,startx_2,starty_2,startz_2,endx_2,endy_2,endz_2,alpha_2,theta_beam_2,theta_q_2;
	
	string pri("primary");

	nGeantParticles=mcparticleListHandle->size();
	nGEANTparticles=nGeantParticles;

	for(int g4pt=0; g4pt<nGEANTparticles; g4pt++)
	  {
	    if(mcparticleListHandle.isValid() && mcparticleListHandle->size()>0)
	      {
		simb::MCParticle const & pPart = mcparticleListHandle->at(g4pt);
		const art::Ptr<simb::MCTruth> mc_truth=fMCTruthMatching->TrackIDToMCTruth(pPart.TrackId());   
		gpdg.push_back(pPart.PdgCode());          //get the PDG ID of all the particles
		gmother.push_back(pPart.Mother());        //get the mother particles ID
		if (pPart.E()<fG4minE ) continue;       
		Bool_t isPrimary = pPart.Process() == pri;
		if( isPrimary && mctruth->NeutrinoSet() && pPart.StatusCode()==1 && pPart.Mother()==0 && mc_truth->Origin()== simb::kBeamNeutrino)  
		  { 
		    fhg4parpdg->push_back(pPart.PdgCode());
		    fhg4parstatus->push_back(pPart.StatusCode());
		    fhg4parpx->push_back(pPart.Px());
		    fhg4parpy->push_back(pPart.Py());
		    fhg4parpz->push_back(pPart.Pz());
		    fhg4partheta->push_back(pPart.Momentum().Theta());
		    fhg4parphi->push_back(pPart.Momentum().Phi());
		    fhg4parp->push_back(pPart.Momentum().Vect().Mag());
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==13) 
		      {
			nmuons++;
			Truth_start_vtxx_muon = pPart.Vx();
			Truth_start_vtxy_muon = pPart.Vy();
			Truth_start_vtxz_muon = pPart.Vz();
			Truth_end_vtxx_muon = pPart.EndPosition()[0];
			Truth_end_vtxy_muon = pPart.EndPosition()[1];
			Truth_end_vtxz_muon = pPart.EndPosition()[2];
			
			//	cout<<"g4 muon: "<<pPart.Px()<<" "<<pPart.Py()<<" "<<pPart.Pz()<<" "<<endl;
			//	cout<<pPart.EndPosition()[0]<<" "<<pPart.EndPosition()[1]<<" "<<pPart.EndPosition()[2]<<" "<<endl;
		      }
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==211) {npions++;}
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==111) {npi0++;}
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==11) {nelectrons++;}
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==2112) {nneutrons++;}
		    if(Truth_CCNC==0 && abs(pPart.PdgCode())==2212) 
		      {
			nprotons++; 
			//	cout<<"Caught proton : "<<g4pt<<endl;
			if(proton_holder_1<0)  
			  {
			    proton_holder_1 = g4pt;
			    e_1 = pPart.E();
			    px_1 = pPart.Px();
			    py_1 = pPart.Py();
			    pz_1 = pPart.Pz();
			    startx_1 = pPart.Vx();
			    starty_1 = pPart.Vy();
			    startz_1 = pPart.Vz();
			    endx_1 = pPart.EndPosition()[0];
			    endy_1 = pPart.EndPosition()[1];
			    endz_1 = pPart.EndPosition()[2];
			    TVector3 temp(px_1,py_1,pz_1);
			    alpha_1 = pPart.Momentum().Vect().Mag()/q_tran.Mag();
			    theta_beam_1 = pPart.Momentum().Theta();
			    theta_q_1= temp.Angle(q_tran);
			  }
			if(proton_holder_1>0&&proton_holder_2<0&&g4pt!=proton_holder_1)  
			  {
			    proton_holder_2 = g4pt;
			    e_2 = pPart.E();
			    px_2 = pPart.Px();
			    py_2 = pPart.Py();
			    pz_2 = pPart.Pz();
			    startx_2 = pPart.Vx();
			    starty_2 = pPart.Vy();
			    startz_2 = pPart.Vz();
			    endx_2 = pPart.EndPosition()[0];
			    endy_2 = pPart.EndPosition()[1];
			    endz_2 = pPart.EndPosition()[2];
			    TVector3 temp(px_2,py_2,pz_2);
			    alpha_2 = pPart.Momentum().Vect().Mag()/q_tran.Mag();
			    theta_beam_2 = pPart.Momentum().Theta();
			    theta_q_2 = temp.Angle(q_tran);
			  }
			if(proton_holder_1>0&&proton_holder_2>0&&g4pt!=proton_holder_1&&g4pt!=proton_holder_2) {proton_holder_more = g4pt;}
		      } 
		  }
		
	      }//end of is the g4 handle is valid and the size is greater than 0;
	  }//end of loop over all the geant 4 particles
       	if(proton_holder_1>0&&proton_holder_2>0&&proton_holder_more<0)
	  {
	    if(cos(theta_q_1)<=cos(theta_q_2))
	      {	      
	       Truth_E_proton_recoil = e_1;
	       Truth_KE_proton_recoil = Truth_E_proton_recoil-proton_mass;
	       Truth_px_proton_recoil = px_1;
	       Truth_py_proton_recoil = py_1;
	       Truth_pz_proton_recoil = pz_1;
	       Truth_start_vtxx_proton_recoil = startx_1;
	       Truth_start_vtxy_proton_recoil = starty_1;
	       Truth_start_vtxz_proton_recoil = startz_1;
	       Truth_end_vtxx_proton_recoil = endx_1;
	       Truth_end_vtxy_proton_recoil = endy_1;
	       Truth_end_vtxz_proton_recoil = endz_1;
	       Truth_theta_wrt_q_proton_recoil = theta_q_1;
	       Truth_theta_wrt_beam_proton_recoil = theta_beam_1;
	       Truth_alpha_proton_recoil = alpha_1;

	       Truth_E_proton_leading =e_2 ;
	       Truth_KE_proton_leading = Truth_E_proton_leading-proton_mass;
	       Truth_px_proton_leading = px_2;
	       Truth_py_proton_leading = py_2;
	       Truth_pz_proton_leading = pz_2;
	       Truth_start_vtxx_proton_leading = startx_2;
	       Truth_start_vtxy_proton_leading = starty_2;
	       Truth_start_vtxz_proton_leading = startz_2;
	       Truth_end_vtxx_proton_leading = endx_2;
	       Truth_end_vtxy_proton_leading = endy_2;
	       Truth_end_vtxz_proton_leading = endz_2;
	       Truth_theta_wrt_q_proton_leading = theta_q_2;
	       Truth_theta_wrt_beam_proton_leading = theta_beam_2;
	       Truth_alpha_proton_leading = alpha_2;	       
	      }
	    else
	      {
	      Truth_E_proton_recoil = e_2;
	       Truth_px_proton_recoil = px_2;
	       Truth_py_proton_recoil = py_2;
	       Truth_pz_proton_recoil = pz_2;
	       Truth_start_vtxx_proton_recoil = startx_2;
	       Truth_start_vtxy_proton_recoil = starty_2;
	       Truth_start_vtxz_proton_recoil = startz_2;
	       Truth_end_vtxx_proton_recoil = endx_2;
	       Truth_end_vtxy_proton_recoil = endy_2;
	       Truth_end_vtxz_proton_recoil = endz_2;
	       Truth_theta_wrt_q_proton_recoil = theta_q_2;
	       Truth_theta_wrt_beam_proton_recoil = theta_beam_2;
	       Truth_alpha_proton_recoil = alpha_2;

	       Truth_E_proton_leading =e_1 ;
	       Truth_px_proton_leading = px_1;
	       Truth_py_proton_leading = py_1;
	       Truth_pz_proton_leading = pz_1;
	       Truth_start_vtxx_proton_leading = startx_1;
	       Truth_start_vtxy_proton_leading = starty_1;
	       Truth_start_vtxz_proton_leading = startz_1;
	       Truth_end_vtxx_proton_leading = endx_1;
	       Truth_end_vtxy_proton_leading = endy_1;
	       Truth_end_vtxz_proton_leading = endz_1;
	       Truth_theta_wrt_q_proton_leading = theta_q_1;
	       Truth_theta_wrt_beam_proton_leading = theta_beam_1;
	       Truth_alpha_proton_leading = alpha_1;
	      }
	  }
	
	Truth_N_muon = nmuons;
	Truth_N_proton = nprotons;
	Truth_N_electron = nelectrons;
	Truth_N_pion = npions + npi0;
	Truth_N_neutron = nneutrons;
       
	if(inFV(Truth_vtxx_neutrino, Truth_vtxy_neutrino, Truth_vtxz_neutrino))
	  {
	    truth_vtx_inFV = 1;
	    cout<<"Truth vtx in FV!"<<endl;
	  }
	else 
	  {
	  truth_vtx_inFV = 0;
	  cout<<"Truth vtx not in FV!"<<endl;
	  }
	/*
	 Fill Truth for all
	 */
	Tr_ALL_Truth->Fill();
       	Int_t Hammer_Flag=Hammer(nmuons, nelectrons, npions, npi0, nprotons);
	if(truth_vtx_inFV == true && Hammer_Flag>0)
	  {
	    Tr_SRC_Truth->Fill();
	  }

      }//end of if MC
    /*--------------------------------------------------

      End of filling Truth tree (GENIE and GEANT4)

      Start with Reco tree

     --------------------------------------------------*/
  
    art::Handle<vector<recob::Track>>   trackVecHandle;
    art::Handle<vector<recob::OpFlash>> flashListHandle;
    art::Handle<vector<recob::MCSFitResult>>   MCSFitHandle;    
  
                                         
    event.getByLabel(fTrackModuleLabel, trackVecHandle);
    event.getByLabel(fTrackMCSFitLabel,     MCSFitHandle);
    vector<art::Ptr<recob::MCSFitResult> >  mcsfitlist;
    art::fill_ptr_vector(mcsfitlist, MCSFitHandle);

    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if(event.getByLabel(fHitsModuleLabel,hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
  
   
    vector<art::Ptr<recob::OpFlash> > flashlist;
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle)) 
      art::fill_ptr_vector(flashlist, flashListHandle);  
  
    vector<art::Ptr<recob::Vertex> >      neutrino_vertex;
    vector<art::Ptr<recob::PFParticle> >  neutrino_vertex_PFP;
  
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(event, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);  
    
    for (lar_pandora::PFParticlesToVertices::iterator it = particlesToVertices.begin(); it != particlesToVertices.end(); ++it)
      {
      art::Ptr<recob::PFParticle> pfParticle = it->first;
      lar_pandora::VertexVector   vertex_v   = it->second;
      if (vertex_v.empty()) continue;
      if (lar_pandora::LArPandoraHelper::IsNeutrino(pfParticle)) 
    	{
    	  if (vertex_v.size() == 1) 
    	    { // require 1 vtx associated to the neutrino PFP
    	      neutrino_vertex.emplace_back(vertex_v[0]);
    	      neutrino_vertex_PFP.emplace_back(pfParticle);
    	    }
    	}
      }

    if( neutrino_vertex.size() > 0 && trackVecHandle.isValid() && trackVecHandle -> size() > 0)
      {
	map< int,vector<int> > VertexTrackCollection;
	//	map< int,vector<int> > SelectedCollection;

	art::FindManyP<recob::PFParticle> trackToPFPartAssns(trackVecHandle,  event, fTrackModuleLabel);
	/*--Flash---*/
	bool flashtag(false);   	
	const recob::OpFlash* flashPtr(0);
	double flashmax = 0;
	for(const auto& opFlash : flashlist)   //loop over all the flashes
	  {
	    if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
	      {
		flashtag=true;
		if (opFlash->TotalPE() > flashmax)              
		  {
		    flash_Time = opFlash->Time(); 
		    flashPtr = opFlash.get();
		    flash_PE = opFlash->TotalPE();
		  }
	      }
	  } //end of loop over all the flash  
	flash_Tag = flashtag;
	/*--Vertex---*/
	int no_nu_vtx=neutrino_vertex.size();
       	//cout<<no_nu_vtx<<" vertex/vertices and "<<trackVecHandle->size()<<" tracks in this event. "<<endl;
	nvertex = no_nu_vtx;
	//	float temp_vtxID;
	float temp_vtxx[no_nu_vtx], temp_vtxy[no_nu_vtx], temp_vtxz[no_nu_vtx];
	//float temp_sum_dist=0;
	//int temp_vtx_inFV;
	int no_track_at_vtx = 0;
	int theone_lu = -1;
	int theone_rac = -1;
	int theone_libo = -1;
	double distmin = 9999;
	float VertexCosTheta = 0.0;
	float deltatheta=-9999;
	int trk_1=-1;
	int trk_2=-1;
	int trk_3=-1;
	for(int i_nuvtx = 0; i_nuvtx < no_nu_vtx; i_nuvtx++)//loop over neutrino vertex
	  {
	    // temp_vtxID=neutrino_vertex[i_nuvtx]->ID();
	    Double_t xyz[3]={};
	    neutrino_vertex[i_nuvtx]->XYZ(xyz);
	    temp_vtxx[i_nuvtx]=xyz[0];
	    temp_vtxy[i_nuvtx]=xyz[1];
	    temp_vtxz[i_nuvtx]=xyz[2];	    
	    if(inFV(temp_vtxx[i_nuvtx], temp_vtxy[i_nuvtx], temp_vtxz[i_nuvtx])) 
	      {
		TVector3 nuvertexPos(temp_vtxx[i_nuvtx], temp_vtxy[i_nuvtx], temp_vtxz[i_nuvtx]);
		int no_track = trackVecHandle->size();
		//	cout<<"inFV vtx ID: "<<temp_vtxID<<endl;
		//	    float temp_startx[no_track], temp_starty[no_track], temp_startz[no_track];
		    for(int i_trk = 0; i_trk < no_track; i_trk++)
		      {
			art::Ptr<recob::Track> track(trackVecHandle,i_trk);
			TVector3 trackPos = track->Vertex();
			TVector3 trackEnd = track->End();
			TVector3 trackdir = track->VertexDirection();
			//	float trackmom=track->VertexMomentum();
			double trackToVertexDist = (trackPos - nuvertexPos).Mag();
			if ((trackEnd - nuvertexPos).Mag() < trackToVertexDist)
			  {
			    trackPos          = track->End();
			    trackEnd          = track->Vertex();
			    trackToVertexDist = (trackPos - nuvertexPos).Mag();
			  }
			if(trackToVertexDist<fMinTrk2VtxDist)
			  {
			    no_track_at_vtx+=1;
			    VertexTrackCollection.insert(pair< int,vector<int> >(i_nuvtx,vector<int>()));
			    VertexTrackCollection.at(i_nuvtx).push_back(i_trk);
			    // 	    cout<<"what's collected: vertex"<<i_nuvtx<<", track: "<<i_trk<<endl;
			  }
		      }//end of looping over all tracks for a neutrino vertex
	      }	  
	  }//end of looping over neutrino vertex
	ntrack = no_track_at_vtx;
	/*---- Loop over all tracks for a neutrino vertex---*/	
	//cout<<ntrack<<" tracks in this event. "<<endl;
	trkf::TrackMomentumCalculator trkmomcal;
	if(no_track_at_vtx==3)
	  {
	    for (auto const& VtxTrack : VertexTrackCollection)
	      {
		int VertexID = VtxTrack.first;
		Reco_vtxID_neutrino=neutrino_vertex[VertexID]->ID();
		Double_t xyz[3]={};
		neutrino_vertex[VertexID]->XYZ(xyz);
		TVector3 nuvertexPos(temp_vtxx[VertexID], temp_vtxy[VertexID], temp_vtxz[VertexID]);
		
		/*  racquel's vertex  */
		float dist12= 0;
		float dist13= 0;
		float dist23 = 0;
		int itertrk=0;
		
		TVector3 trk1Start;
		TVector3 trk2Start;
		TVector3 trk3Start;
		for (auto const& TrackID : VtxTrack.second)
		  {
		    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
		    if(itertrk ==0)  trk1Start = GetStartTrack(track,nuvertexPos);
		    else if(itertrk ==1)  trk2Start = GetStartTrack(track,nuvertexPos);
		    else if(itertrk ==2)  trk3Start = GetStartTrack(track,nuvertexPos);
		    itertrk++;
		  }//loop track
		dist12 = (trk1Start - trk2Start).Mag();
		dist13 = (trk1Start - trk3Start).Mag();
		dist23 = (trk2Start - trk3Start).Mag();
		double dist123 = dist12 + dist13 + dist23;
		if (distmin>dist123) 
		  {
		    distmin = dist123;
		    theone_rac = VertexID;
		    //	  cout<<"Rac: "<<VertexID<<" "<<theone_rac<<" "<<distmin<<endl;
		  }
		
		/*  Libo's vertex  */
		
                float WeightedCosTheta = 0.0;
                float NormFactor = 0.0;
		
                for (auto const& TrackID : VtxTrack.second)
		  {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
                    double TrackRange= track->Length();
                    NormFactor += TrackRange;
                    WeightedCosTheta += TrackRange*cos(track->Theta());
		  }// track ID loop
                // Make average
                WeightedCosTheta /= NormFactor;
                if (fabs(WeightedCosTheta) > VertexCosTheta)
		  {
                    theone_libo = VertexID;     
		    //	    cout<<"Libo: "<<VertexID<<" "<<theone_libo<<" "<<fabs(WeightedCosTheta)<<endl;
		  }
		
		/*  Lu's vertex  */
		float delta_sum = 0.0;
		//		int count=0;
		//  NumofRecoTrack=0;
		//	TVector3 nuvertexPos(temp_vtxx, temp_vtxy, temp_vtxz);
                for (auto const& TrackID : VtxTrack.second)
		  {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
		    //
		    if(trk_1==-1&&trk_2==-1&&trk_3==-1)  trk_1 = TrackID;
		    if(trk_1!=-1&&trk_2==-1&&trk_3==-1&&TrackID!=trk_1)  trk_2 = TrackID;
		    if(trk_1!=-1&&trk_2!=-1&&trk_3==-1&&TrackID!=trk_1&&TrackID!=trk_2)  trk_3 = TrackID;
		    //   count++;
		    //
		    TVector3 trackPos = track->Vertex();
		    TVector3 trackEnd = track->End();
		    TVector3 trackdir = track->VertexDirection();
		    //      float trackmom=track->VertexMomentum();
		    double trackToVertexDist = (trackPos - nuvertexPos).Mag();
		    if ((trackEnd - nuvertexPos).Mag() < trackToVertexDist)
		      { 
			trackPos          = track->End();
			trackEnd          = track->Vertex();
			//	 trackToVertexDist = (trackPos - nuvertexPos).Mag();
		      }
		    
		    double delta_test= pow((acos((trackPos.Z()-temp_vtxz[VertexID])/(trackPos-nuvertexPos).Mag())-track->Theta()), 2);
		    delta_sum += delta_test;
		    //	    cout<<"Collected track ID: "<<TrackID<<" "<<trk_1<<" "<<trk_2<<" "<<trk_3<<endl;
		  }// track ID loop
		
                if (delta_sum > deltatheta)
		  {
                    theone_lu = VertexID;                  
		    //	    cout<<"Lu: "<<VertexID<<" "<<theone_lu<<" "<<delta_sum<<endl;
		  }
		
	      }//loop vertex
	    //  SelectedCollection.insert(pair< int,vector<int> >(i_nuvtx,vector<int>()));
	    //SelectedCollection.at(i_nuvtx).push_back(i_trk);
	  }// if 3 tracks
	
       	Reco_vtxx_neutrino=temp_vtxx[theone_lu];
	Reco_vtxy_neutrino=temp_vtxy[theone_lu];
	Reco_vtxz_neutrino=temp_vtxz[theone_lu];
	if(inFV(Reco_vtxx_neutrino, Reco_vtxy_neutrino, Reco_vtxz_neutrino))
	   Reco_vtx_inFV = 1;
	else
	   Reco_vtx_inFV = 0;
	/*  for vertex comparison  */
	TVector3 truth_vertex(Truth_vtxx_neutrino,Truth_vtxy_neutrino,Truth_vtxz_neutrino);
	TVector3 lu_vertex(temp_vtxx[theone_lu],temp_vtxy[theone_lu],temp_vtxz[theone_lu]);
	TVector3 libo_vertex(temp_vtxx[theone_libo],temp_vtxy[theone_libo],temp_vtxz[theone_libo]);
	TVector3 rac_vertex(temp_vtxx[theone_rac],temp_vtxy[theone_rac],temp_vtxz[theone_rac]);

	vtx_truth_reco_lu=(truth_vertex-lu_vertex).Mag();
	vtx_truth_reco_libo=(truth_vertex-libo_vertex).Mag();
	vtx_truth_reco_rac=(truth_vertex-rac_vertex).Mag();
	/*  track mometum  */

	Reco_TMD_muon_tag_muon=false;

	//	bool mcs_isBestFwd = false;
     
	if(no_track_at_vtx==3)
	  {
	    int trklist[3]={trk_1,trk_2,trk_3};
	    double TrackLength[3];
	    double trackstartx[3];
	    double trackstarty[3];
	    double trackstartz[3];
	    double trackendx[3];
	    double trackendy[3];
	    double trackendz[3];
	    TVector3 trackPos[3];
	    TVector3 trackEnd[3];
	    double TrackRange[3];
	    // float trackMom[3];
	    float trackTheta[3];
	    float trackPhi[3];
	    float tracktrundqdx[3];
	    float tracktrundqdx_0[3];
	    float tracktrundqdx_1[3];
	  
	    float tracknhits[3];
	    float tracknhits_0[3];
	    float tracknhits_1[3];
	    int track_matchedtrue_PDG[3];
	    int track_matchedtrue_origin[3];
	    float track_matchedtrue_p[3];	    
	    double trackMom_MCS[3];
	    double temp_length[3];
	    //int temp_id[3];
	    int sorted_id[3]={0,1,2};
	   
	    
	    for (int i=0;i<3;i++)
	      {
		
		art::Ptr<recob::Track> track(trackVecHandle,trklist[i]);
		temp_length[i]= track->Length();
		//	cout<<"before sorted: "<<i<<" "<<temp_length[i]<<endl;
		//	 temp_id [i]= trklist[i];
	      }
	    /* sorting 3 tracks */
	    //double temp_for_sort;
	    if(temp_length[0]<temp_length[1]){sorted_id[0]=1;sorted_id[1]=0;}
	    if(temp_length[1]<temp_length[2]){sorted_id[1]=2;sorted_id[2]=1;}
	    if(temp_length[0]<temp_length[2]){sorted_id[0]=2;sorted_id[1]=0;}
	    
	    
	    
	    vector<double>track_dEdx[3];
	    vector<double>track_dQdx[3];
	    vector<double>track_ResidualRange[3];
	    vector<double>track_dEdx_0[3];
	    vector<double>track_dQdx_0[3];
	    vector<double>track_ResidualRange_0[3];
	    vector<double>track_dEdx_1[3];
	    vector<double>track_dQdx_1[3];
	    vector<double>track_ResidualRange_1[3];
	    float PIDA[3];
	    float trkpidchi[3];
	    float trkpidchi_proton[3];
	    float trkpidchi_kaon[3];
	    float trkpidchi_pion[3];
	    float trkpidchi_muon[3];
	    float PIDA_0[3];
	    //	    float trkpidchi_0[3];
	    float trkpidchi_proton_0[3];
	    float PIDA_1[3];
	    //    float trkpidchi_1[3];
	    float trkpidchi_proton_1[3];

	    int   best_plane[3];
	    int   best_nhits[3];
	    float best_tmdqdx[3];
	    float best_PIDA[3];
	    //  float best_PIDchi2[3];
	    float best_PIDchi2proton[3];
	    vector<double>best_dEdx[3];
	    vector<double>best_dQdx[3];
	    vector<double>best_ResidualRange[3];
	    bool  best_TMD_tag[3];
	
	    for (int i=0;i<3;i++)
	      {
		Reco_TMD_muon_tag_track[i]=false;
		art::Ptr<recob::Track> track(trackVecHandle,trklist[i]);
		trackMom_MCS[i] = mcsfitlist[trklist[i]]->bestMomentum();
		//	mcs_isBestFwd = mcsfitlist[trklist[i]]->isBestFwd();
		//if (!mcs_isBestFwd){cout << "muon candidate MCS fits better backwards" << endl;}
		//cout<<"Ok here?"<<endl;
		art::FindManyP<recob::Hit> fmh(trackVecHandle, event, fTrackModuleLabel);
		vector<art::Ptr<recob::Hit>> hits = fmh.at(track.key());
		art::FindManyP<anab::Calorimetry> calos_from_track(trackVecHandle, event, fCalorimetryModuleLabel);
		vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(track.key());
		art::FindMany<anab::ParticleID> fmpid(trackVecHandle, event, fParticleIDModuleLabel);
		std::vector<const anab::ParticleID*> pids = fmpid.at(track.key());
		for (auto c : calos) 
		  {
		    if (!c) continue;
		    if (!c->PlaneID().isValid) continue;
		    int planenum = c->PlaneID().Plane;
		    if (planenum == 2)
		      {
			
			track_dEdx[i]= c->dEdx();
			track_dQdx[i]= c->dQdx();
			track_ResidualRange[i]=c->ResidualRange();
		      }
		    if (planenum == 0)
		      {
			
			track_dEdx_0[i]= c->dEdx();
			track_dQdx_0[i]= c->dQdx();
			track_ResidualRange_0[i]=c->ResidualRange();
		      }
		    if (planenum == 1)
		      {
			
			track_dEdx_1[i]= c->dEdx();
			track_dQdx_1[i]= c->dQdx();
			track_ResidualRange_1[i]=c->ResidualRange();
		      }
		    
		  }
		
       
		for (auto p : pids)
		  {
		    if (!p) continue;
		    if (!p->PlaneID().isValid) continue;
		    int planenum = p->PlaneID().Plane;
		    if(planenum ==0) 
		      {
			//	trkpidchi_0[i] = p->MinChi2();
			trkpidchi_proton_0[i] = p->Chi2Proton();
			PIDA_0[i] = p->PIDA();
		      }
		     if(planenum ==1) 
		      {
			//	trkpidchi_1[i] = p->MinChi2();
			trkpidchi_proton_1[i] = p->Chi2Proton();
			PIDA_1[i] = p->PIDA();
		      }
		      if(planenum ==2) 
		      {
			trkpidchi[i] = p->MinChi2();
			trkpidchi_proton[i] = p->Chi2Proton();
			trkpidchi_kaon[i] = p->Chi2Kaon();
			trkpidchi_pion[i] = p->Chi2Pion();
			trkpidchi_muon[i] = p->Chi2Muon();
			PIDA[i] = p->PIDA();
		      }

		    
		  }
		
		
		 /* truth information for muon tracks */
		if(isMC)
		  {
		    double tmpEfracm=0;
		    double tmpEcompletm=0;
		    const simb::MCParticle *mparticle;
		    truthMatcher(hitlist, hits, mparticle, tmpEfracm, tmpEcompletm );
		    const art::Ptr<simb::MCTruth> MCtruth = fMCTruthMatching->ParticleToMCTruth(mparticle);
		    track_matchedtrue_origin[i]=MCtruth->Origin();
		    if(mparticle)
		      {
			track_matchedtrue_PDG[i] = mparticle->PdgCode();
			track_matchedtrue_p[i] = sqrt(mparticle->Px()*mparticle->Px()+mparticle->Py()*mparticle->Py()+mparticle->Pz()*mparticle->Pz());
		      }
		  }
		
		tracknhits[i] = calos[2] -> dEdx().size();
		tracknhits_0[i] = calos[0] -> dEdx().size();
		tracknhits_1[i] = calos[1] -> dEdx().size();
		tracktrundqdx[i]=GetDqDxTruncatedMean_alt(calos,0);      
		tracktrundqdx_0[i]=GetDqDxTruncatedMean_0(calos);      
		tracktrundqdx_1[i]=GetDqDxTruncatedMean_1(calos);      
		TrackLength[i]= track->Length();
		trackPos[i] = track->Vertex();
		trackEnd[i] = track->End();
		TrackRange[i] = (trackPos[i]-trackEnd[i]).Mag();
		//	trackMom[i]=track->VertexMomentum();
		trackTheta[i]=track->Theta();
		trackPhi[i]=track->Phi();       
		trackstartx[i]=trackPos[i].x();
		trackstarty[i]=trackPos[i].y();
		trackstartz[i]=trackPos[i].z();
		trackendx[i]=trackEnd[i].x();
		trackendy[i]=trackEnd[i].y();
		trackendz[i]=trackEnd[i].z();
	
		
		if(tracknhits[i]>=tracknhits_0[i]&&tracknhits[i]>=tracknhits_1[i])
		  {
		    best_plane[i]=2;
		    best_nhits[i]=tracknhits[i];
		    best_tmdqdx[i]=tracktrundqdx[i];
		    best_PIDA[i]=PIDA[i];
		    //  best_PIDchi2[i]=trkpidchi[i];
		    best_PIDchi2proton[i]=trkpidchi_proton[i];
		    best_dEdx[i]=track_dEdx[i];
		    best_dQdx[i]=track_dQdx[i];
		    best_ResidualRange[i]=track_dQdx[i];
		  }
		else if(tracknhits_1[i]>tracknhits_0[i]&&tracknhits_1[i]>=tracknhits[i])
		  {
		    best_plane[i]=1;
		    best_nhits[i]=tracknhits_1[i];
		    best_tmdqdx[i]=tracktrundqdx_1[i];
		    best_PIDA[i]=PIDA_1[i];
		    //  best_PIDchi2[i]=trkpidchi_1[i];
		    best_PIDchi2proton[i]=trkpidchi_proton_1[i];
		    best_dEdx[i]=track_dEdx_1[i];
		    best_dQdx[i]=track_dQdx_1[i];
		    best_ResidualRange[i]=track_dQdx_1[i];
		  }
		else
		  {
		    best_plane[i]=0;
		    best_nhits[i]=tracknhits_0[i];
		    best_tmdqdx[i]=tracktrundqdx_0[i];
		    best_PIDA[i]=PIDA_0[i];
		    //  best_PIDchi2[i]=trkpidchi_0[i];
		    best_PIDchi2proton[i]=trkpidchi_proton_0[i];
		    best_dEdx[i]=track_dEdx_0[i];
		    best_dQdx[i]=track_dQdx_0[i];
		    best_ResidualRange[i]=track_dQdx_0[i];
		  }
		

		if(MIPConsistency(tracktrundqdx[i],TrackLength[i]))
		  {
		    Reco_TMD_muon_tag_track[i]=true;
		    //  cout<<"muon? "<<Reco_muon_tag_track[i]<<endl;
		  }
		
		if(MIPConsistency(best_tmdqdx[i],TrackLength[i]))
		  {
		    best_TMD_tag[i]=true;

		  }
		else{ best_TMD_tag[i]=false;}

	
		Reco_startx_track[i]=trackstartx[i];
		Reco_starty_track[i]=trackstarty[i];
		Reco_startz_track[i]=trackstartz[i];
		Reco_nhits_track[i]=tracknhits[i];
		Reco_TrunMean_dQdx_track[i]=tracktrundqdx[i];
		Reco_TrunMean_dQdx_scaled_track[i]=tracktrundqdx[i]*fdQdx_scale;		
		Reco_range_track[i]=TrackRange[i];
		Reco_length_track[i]=TrackLength[i];
		Reco_MC_PDG_track[i]=track_matchedtrue_PDG[i];
		Reco_MC_p_track[i]=track_matchedtrue_p[i];
		Reco_MC_origin_track[i]= track_matchedtrue_origin[i];
	
		Reco_endx_track[i]=trackendx[i];
		Reco_endy_track[i]=trackendy[i];
		Reco_endz_track[i]=trackendz[i];
	      } //end of loop over all the tracks associated to the vertex candidate

	     /* muon candidate information */
	     Reco_p_range_muon=trkmomcal.GetTrackMomentum(TrackLength[sorted_id[0]],13);
	     Reco_p_MCS_muon= trackMom_MCS[sorted_id[0]];
	     Reco_range_muon = TrackRange[sorted_id[0]];
	     Reco_length_muon = TrackLength[sorted_id[0]];
	     Reco_theta_muon = trackTheta[sorted_id[0]];
	     Reco_phi_muon = trackPhi[sorted_id[0]];   
	     Reco_MC_p_muon = track_matchedtrue_p[sorted_id[0]];
	     Reco_MC_PDG_muon =   track_matchedtrue_PDG[sorted_id[0]];
	     Reco_MC_origin_muon =  track_matchedtrue_origin[sorted_id[0]];
	     Reco_TMD_muon_tag_muon =  Reco_TMD_muon_tag_track[sorted_id[0]];
	     Reco_start_vtxx_muon= trackstartx[sorted_id[0]];
	     Reco_start_vtxy_muon= trackstarty[sorted_id[0]];
	     Reco_start_vtxz_muon= trackstartz[sorted_id[0]];
	     Reco_end_vtxx_muon= trackendx[sorted_id[0]];
	     Reco_end_vtxy_muon= trackendy[sorted_id[0]];
	     Reco_end_vtxz_muon= trackendz[sorted_id[0]];
	     Reco_TrunMean_dQdx_muon= tracktrundqdx[sorted_id[0]];
	     Reco_TrunMean_dQdx_scaled_muon= tracktrundqdx[sorted_id[0]]*fdQdx_scale;
	     Reco_nhits_muon=tracknhits[sorted_id[0]];
	     Reco_nhits_0_muon=tracknhits_0[sorted_id[0]];
	     Reco_nhits_1_muon=tracknhits_1[sorted_id[0]];
	     Reco_dQdx_muon= track_dQdx[sorted_id[0]];
	     Reco_dEdx_muon= track_dEdx[sorted_id[0]];
	     Reco_ResidualRange_muon= track_ResidualRange[sorted_id[0]];
	     Reco_dQdx_0_muon= track_dQdx_0[sorted_id[0]];
	     Reco_dEdx_0_muon= track_dEdx_0[sorted_id[0]];
	     Reco_ResidualRange_0_muon= track_ResidualRange_0[sorted_id[0]];
	     Reco_dQdx_1_muon= track_dQdx_1[sorted_id[0]];
	     Reco_dEdx_1_muon= track_dEdx_1[sorted_id[0]];
	     Reco_ResidualRange_1_muon= track_ResidualRange_1[sorted_id[0]];
	     Reco_PID_chi2_muon=trkpidchi[sorted_id[0]];
	     Reco_PID_chi2_proton_muon=trkpidchi_proton[sorted_id[0]];
	     Reco_PID_chi2_kaon_muon=trkpidchi_kaon[sorted_id[0]];
	     Reco_PID_chi2_pion_muon=trkpidchi_pion[sorted_id[0]];
	     Reco_PID_chi2_muon_muon=trkpidchi_muon[sorted_id[0]];
	     Reco_PIDA_muon=PIDA[sorted_id[0]];
	     /*!!!!!!!!*/     
	     Reco_best_plane_muon= best_plane[sorted_id[0]];
	     Reco_best_nhits_muon= best_nhits[sorted_id[0]];
	     Reco_best_TrunMean_dQdx_muon =  best_tmdqdx[sorted_id[0]];
	     Reco_best_TrunMean_dQdx_scaled_muon = best_tmdqdx[sorted_id[0]]*fdQdx_scale;
	     Reco_best_TMD_muon_tag_muon=best_TMD_tag[sorted_id[0]];
	     Reco_best_PIDA_muon= best_PIDA[sorted_id[0]];
	     Reco_best_PID_chi2_proton_muon= best_PIDchi2proton[sorted_id[0]];
	     Reco_best_dQdx_muon= best_dQdx[sorted_id[0]];
	     Reco_best_dEdx_muon=best_dEdx[sorted_id[0]];
	     Reco_best_ResidualRange_muon=best_ResidualRange[sorted_id[0]];
	     
	     
	     if(flash_Tag==true)
	       {
		 double flstrkdist=GetFlashTrackDist(flashPtr->ZCenter(), Reco_start_vtxz_muon, Reco_end_vtxz_muon);
		 Reco_flash_track_distance=flstrkdist;
		 if(flstrkdist < fFlashWidth) 
		   {
		     trackflash_Tag=true; 
		   }
	       }
	     if(inFV(Reco_start_vtxx_muon, Reco_start_vtxy_muon, Reco_start_vtxz_muon)&&inFV(Reco_end_vtxx_muon,Reco_end_vtxy_muon,Reco_end_vtxz_muon))
	       {
		 Reco_trk_inFV_muon=1;	   
		 Reco_p_muon= Reco_p_range_muon;
	       }
	     else
	       {
		 Reco_trk_inFV_muon=0;
		 Reco_p_muon= Reco_p_MCS_muon;
	       }
	     
	     /* Long proton block */
	     Reco_range_LP = TrackRange[sorted_id[1]];
	     Reco_length_LP = TrackLength[sorted_id[1]];
	     Reco_theta_LP = trackTheta[sorted_id[1]];
	     Reco_phi_LP = trackPhi[sorted_id[1]];   
	     Reco_MC_p_LP = track_matchedtrue_p[sorted_id[1]];
	     Reco_MC_PDG_LP =   track_matchedtrue_PDG[sorted_id[1]];
	     Reco_MC_origin_LP =  track_matchedtrue_origin[sorted_id[1]];
	     Reco_TMD_muon_tag_LP =  Reco_TMD_muon_tag_track[sorted_id[1]];
	     Reco_start_vtxx_LP= trackstartx[sorted_id[1]];
	     Reco_start_vtxy_LP= trackstarty[sorted_id[1]];
	     Reco_start_vtxz_LP= trackstartz[sorted_id[1]];
	     Reco_end_vtxx_LP= trackendx[sorted_id[1]];
	     Reco_end_vtxy_LP= trackendy[sorted_id[1]];
	     Reco_end_vtxz_LP= trackendz[sorted_id[1]];
	     Reco_TrunMean_dQdx_LP= tracktrundqdx[sorted_id[1]];
	     Reco_TrunMean_dQdx_scaled_LP= tracktrundqdx[sorted_id[1]]*fdQdx_scale;
	     Reco_nhits_LP=tracknhits[sorted_id[1]];
	      Reco_nhits_0_LP=tracknhits_0[sorted_id[1]];
	      Reco_nhits_1_LP=tracknhits_1[sorted_id[1]];
	     Reco_dQdx_LP= track_dQdx[sorted_id[1]];
	     Reco_dEdx_LP= track_dEdx[sorted_id[1]];
	     Reco_ResidualRange_LP= track_ResidualRange[sorted_id[1]];
	    Reco_dQdx_0_LP= track_dQdx_0[sorted_id[1]];
	     Reco_dEdx_0_LP= track_dEdx_0[sorted_id[1]];
	     Reco_ResidualRange_0_LP= track_ResidualRange_0[sorted_id[1]];
	    Reco_dQdx_1_LP= track_dQdx_1[sorted_id[1]];
	     Reco_dEdx_1_LP= track_dEdx_1[sorted_id[1]];
	     Reco_ResidualRange_1_LP= track_ResidualRange_1[sorted_id[1]];
	     Reco_PID_chi2_LP=trkpidchi[sorted_id[1]];
	     Reco_PID_chi2_proton_LP=trkpidchi_proton[sorted_id[1]];
	     Reco_PID_chi2_kaon_LP=trkpidchi_kaon[sorted_id[1]];
	     Reco_PID_chi2_pion_LP=trkpidchi_pion[sorted_id[1]];
	     Reco_PID_chi2_muon_LP=trkpidchi_muon[sorted_id[1]];
	     Reco_PIDA_LP=PIDA[sorted_id[1]];
	     Reco_best_plane_LP= best_plane[sorted_id[1]];
	     Reco_best_nhits_LP= best_nhits[sorted_id[1]];
	     Reco_best_TrunMean_dQdx_LP =  best_tmdqdx[sorted_id[1]];
	     Reco_best_TrunMean_dQdx_scaled_LP = best_tmdqdx[sorted_id[1]]*fdQdx_scale;
	     Reco_best_TMD_muon_tag_LP=best_TMD_tag[sorted_id[1]];
	     Reco_best_PIDA_LP= best_PIDA[sorted_id[1]];
	     Reco_best_PID_chi2_proton_LP= best_PIDchi2proton[sorted_id[1]];
	     Reco_best_dQdx_LP= best_dQdx[sorted_id[1]];
	     Reco_best_dEdx_LP=best_dEdx[sorted_id[1]];
	     Reco_best_ResidualRange_LP=best_ResidualRange[sorted_id[1]];
	     
	
	     if(inFV(Reco_start_vtxx_LP, Reco_start_vtxy_LP, Reco_start_vtxz_LP)&&inFV(Reco_end_vtxx_LP,Reco_end_vtxy_LP,Reco_end_vtxz_LP))
	       {
		 Reco_trk_inFV_LP=1;	   
		 Reco_p_LP=trkmomcal.GetTrackMomentum(TrackLength[sorted_id[1]],2212);
	       }
	     else
	       {
		 Reco_trk_inFV_LP=0;
	       }
	     
	     /* Short proton block */
	     Reco_range_SP = TrackRange[sorted_id[2]];
	     Reco_length_SP = TrackLength[sorted_id[2]];
	     Reco_theta_SP = trackTheta[sorted_id[2]];
	     Reco_phi_SP = trackPhi[sorted_id[2]];   
	     Reco_MC_p_SP = track_matchedtrue_p[sorted_id[2]];
	     Reco_MC_PDG_SP =   track_matchedtrue_PDG[sorted_id[2]];
	     Reco_MC_origin_SP =  track_matchedtrue_origin[sorted_id[2]];
	     Reco_TMD_muon_tag_SP =  Reco_TMD_muon_tag_track[sorted_id[2]];
	     Reco_start_vtxx_SP= trackstartx[sorted_id[2]];
	     Reco_start_vtxy_SP= trackstarty[sorted_id[2]];
	     Reco_start_vtxz_SP= trackstartz[sorted_id[2]];
	     Reco_end_vtxx_SP= trackendx[sorted_id[2]];
	     Reco_end_vtxy_SP= trackendy[sorted_id[2]];
	     Reco_end_vtxz_SP= trackendz[sorted_id[2]];
	     Reco_TrunMean_dQdx_SP= tracktrundqdx[sorted_id[2]];
	     Reco_TrunMean_dQdx_scaled_SP= tracktrundqdx[sorted_id[2]]*fdQdx_scale;
	     Reco_nhits_SP=tracknhits[sorted_id[2]];
	     Reco_dQdx_SP= track_dQdx[sorted_id[2]];
	     Reco_dEdx_SP= track_dEdx[sorted_id[2]];
	     Reco_ResidualRange_SP= track_ResidualRange[sorted_id[2]];
	     Reco_nhits_0_SP=tracknhits_0[sorted_id[2]];
	     Reco_dQdx_0_SP= track_dQdx_0[sorted_id[2]];
	     Reco_dEdx_0_SP= track_dEdx_0[sorted_id[2]];
	     Reco_ResidualRange_0_SP= track_ResidualRange_0[sorted_id[2]];
	     Reco_nhits_1_SP=tracknhits_1[sorted_id[2]];
	     Reco_dQdx_1_SP= track_dQdx_1[sorted_id[2]];
	     Reco_dEdx_1_SP= track_dEdx_1[sorted_id[2]];
	     Reco_ResidualRange_1_SP= track_ResidualRange_1[sorted_id[2]];
	     Reco_PID_chi2_SP=trkpidchi[sorted_id[2]];
	     Reco_PID_chi2_proton_SP=trkpidchi_proton[sorted_id[2]];
	     Reco_PID_chi2_kaon_SP=trkpidchi_kaon[sorted_id[2]];
	     Reco_PID_chi2_pion_SP=trkpidchi_pion[sorted_id[2]];
	     Reco_PID_chi2_muon_SP=trkpidchi_muon[sorted_id[2]];
	     Reco_PIDA_SP=PIDA[sorted_id[2]];
	     Reco_best_plane_SP= best_plane[sorted_id[2]];
	     Reco_best_nhits_SP= best_nhits[sorted_id[2]];
	     Reco_best_TrunMean_dQdx_SP =  best_tmdqdx[sorted_id[2]];
	     Reco_best_TrunMean_dQdx_scaled_SP = best_tmdqdx[sorted_id[2]]*fdQdx_scale;
	     Reco_best_TMD_muon_tag_SP=best_TMD_tag[sorted_id[2]];
	     Reco_best_PIDA_SP= best_PIDA[sorted_id[2]];
	     Reco_best_PID_chi2_proton_SP= best_PIDchi2proton[sorted_id[2]];
	     Reco_best_dQdx_SP= best_dQdx[sorted_id[2]];
	     Reco_best_dEdx_SP=best_dEdx[sorted_id[2]];
	     Reco_best_ResidualRange_SP=best_ResidualRange[sorted_id[2]];
	     if(inFV(Reco_start_vtxx_SP, Reco_start_vtxy_SP, Reco_start_vtxz_SP)&&inFV(Reco_end_vtxx_SP,Reco_end_vtxy_SP,Reco_end_vtxz_SP))
	       {
		 Reco_trk_inFV_SP=1;	   
		 Reco_p_SP=trkmomcal.GetTrackMomentum(TrackLength[sorted_id[2]],2212);
	       }
	     else
	       {
		 Reco_trk_inFV_SP=0;
	       }
	  }
	if(flash_Tag)
	  {
	    Tr_SRC_Reco->Fill();
	  }//selected sample reco
      }//valid neutrino vertex and valid track
    Tr_ALL_Reco->Fill();
    ClearLocalData();
    return;
  }// Analyze
  DEFINE_ART_MODULE(SRCAna)
} // namespace  SRCAna
#endif //  SRCAna_Module
