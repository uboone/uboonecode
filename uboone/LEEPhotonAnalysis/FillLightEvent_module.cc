////////////////////////////////////////////////////////////////////////
// Class:       FillLightEvent
// Plugin Type: analyzer (art v2_05_01)
// File:        FillLightEvent_module.cc
//
// Generated at Wed Feb 14 06:23:27 2018 by Robert Murrells using cetskelgen
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

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "uboone/EventWeight/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "../LLBasicTool/GeoAlgo/GeoAABox.h"
#include "../LLBasicTool/GeoAlgo/GeoAlgo.h"

#include "RecoMCMatching.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "FilterSignal.h"
#include "EnergyHelper.h"
#include "EnergyHelperNew.h"

#include "TTree.h"


class FillLightEvent;


class FillLightEvent : public art::EDAnalyzer {

public:

  explicit FillLightEvent(fhicl::ParameterSet const & p);
  FillLightEvent(FillLightEvent const &) = delete;
  FillLightEvent(FillLightEvent &&) = delete;
  FillLightEvent & operator = (FillLightEvent const &) = delete;
  FillLightEvent & operator = (FillLightEvent &&) = delete;

  
  void Reconfigure(fhicl::ParameterSet const & p);

  void SetupTrees();

  void beginSubRun(art::SubRun const & sr);
  void endJob();

  void ResetEvent();

  bool PassedSWTrigger(art::Event const & e, std::string const & swtrigger_product) const;
  void FillSWTriggerVectors(art::Event const & e);
  void FillRecoOpFlashVectors(art::Event const & e, int const producer_index, size_t const opflash_index_offset);
  void FillRecoOpFlashVectors(art::Event const & e);
  void FillRecoHitVectors(art::Event const & e);
  void FillRecoHitVectors(art::Event const & e, int const producer_index, size_t const hit_index_offset);
  std::pair<std::vector<double>,std::vector<double>> GetTrackCaloInfo(art::Event const & e,
								      std::string const & track_producer,
								      std::string const & track_calo_producer,
								      size_t const track_index,
								      double & energy);
  void FillRecoTrackVectors(art::Event const & e);
  void FillRecoTrackVectors(art::Event const & e, int const producer_index, size_t const track_index_offset, size_t const hit_index_offset);
  std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> GetTrackHelperEnergy(art::Event const & e,
											      std::string const & track_producer,
											      size_t const track_index,
											      std::vector<double> & track_energy_dedx);
  std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> GetShowerHelperEnergy(art::Event const & e,
											       std::string const & shower_producer,
											       size_t const shower_index,
											       bool const energy_helper_new = false);
  void FillRecoShowerVectors(art::Event const & e);
  void FillRecoShowerVectors(art::Event const & e, int const producer_index, size_t const shower_index_offset, size_t const hit_index_offset);
  void FillPandora(art::Event const & e);

  void FillWeights(art::Event const & e);

  void FillGenieParticleVectors(art::Event const & e);
  void FillGenieParticleVectors(art::Event const & e,
				size_t const mct_index);
  void FillMCParticleVectors(art::Event const & e);
  void FillMCTrackVectors(art::Event const & e);
  void FillMCShowerVectors(art::Event const & e);
  void FillTruth(art::Event const & e);

  void GetDeltaMCShowerMCTrackIndices(art::Event const & e,
				      size_t const delta_rad_mct_index);
  
  void analyze(art::Event const & e) override;

private:

  bool fheavy;
  bool fmc;
  bool fuse_eventweight;
  std::string fpot_producer;

  std::string fswtrigger_product;
  std::vector<std::string> fopflash_producers;
  std::vector<int> fopflash_producer_indices;
  size_t fopflashp_size;

  std::vector<std::string> fhit_producers;
  std::vector<int> fhit_producer_indices;
  size_t fhitp_size;
  std::vector<std::string> ftrack_producers;
  std::string ftrack_calo_producer;
  std::vector<int> ftrack_producer_indices;
  size_t ftrackp_size;
  std::vector<std::string> fshower_producers;
  std::vector<int> fshower_producer_indices;
  size_t fshowerp_size;

  std::vector<std::string> frmcmassociation_producers;
  size_t frmcm_size;
  RecoMCMatching * frmcm_first;

  geoalgo::AABox ftpc_volume;
  double foffset;
  geoalgo::AABox ffiducial_volume;
  geoalgo::GeoAlgo const falgo;

  std::vector<RecoMCMatching> frmcm;
  int fmc_type_shower;
  int fmc_type_track;
  int fmc_type_particle;

  FilterSignal ffs;
  lee::EnergyHelper fenergyHelper;
  lee::EnergyHelperNew fenergyHelperNew;

  bool fverbose;

  TTree * fpot_tree;
  int fnumber_of_events;
  double fpot;

  TTree * fevent_tree;

  //All
  int frun_number;
  int fsubrun_number;
  int fevent_number;

  //Software Trigger
  //std::vector<int> fswtrigger_producer_index;
  std::vector<std::string> fswtrigger_getListOfAlgorithms;
  std::vector<bool> fswtrigger_passedAlgo;
  std::vector<bool> fswtrigger_passedPrescaleAlgo;
  std::vector<bool> fswtrigger_vetoAlgo;
  int fswtrigger_passedAlgos;
  int fswtrigger_vetoAlgos;
  int fswtrigger_passedPrescaleAlgos;
  int fswtrigger_passed_swtrigger;

  //Reco Flash
  std::vector<double> freco_opflash_producer_index;
  std::vector<double> freco_opflash_Time;
  std::vector<double> freco_opflash_TimeWidth;
  std::vector<double> freco_opflash_AbsTime;
  std::vector<int> freco_opflash_Frame;
  //std::vector<std::vector<double>> freco_opflash_PEs;
  std::vector<double> freco_opflash_YCenter;
  std::vector<double> freco_opflash_YWidth;
  std::vector<double> freco_opflash_ZCenter;
  std::vector<double> freco_opflash_ZWidth;
  std::vector<bool> freco_opflash_InBeamFrame;
  std::vector<int> freco_opflash_OnBeamFrame;
  std::vector<std::vector<double>> freco_opflash_WireCenters;
  std::vector<std::vector<double>> freco_opflash_WireWidths;
  std::vector<double> freco_opflash_TotalPE;
  std::vector<double> freco_opflash_FastToTotal;

  //Reco Hit
  std::vector<int> freco_hit_producer_index;
  std::vector<int> freco_hit_StartTick;
  std::vector<int> freco_hit_EndTick;
  std::vector<float> freco_hit_PeakTime;
  std::vector<float> freco_hit_SigmaPeakTime;
  std::vector<float> freco_hit_RMS;
  std::vector<float> freco_hit_PeakAmplitude;
  std::vector<float> freco_hit_SigmaPeakAmplitude;
  std::vector<float> freco_hit_SummedADC;
  std::vector<float> freco_hit_Integral;
  std::vector<float> freco_hit_SigmaIntegral;
  std::vector<int> freco_hit_Multiplicity;
  std::vector<int> freco_hit_LocalIndex;
  std::vector<float> freco_hit_GoodnessOfFit;
  std::vector<int> freco_hit_DegreesOfFreedom;
  std::vector<int> freco_hit_View;
  std::vector<int> freco_hit_SignalType;
  std::vector<int> freco_hit_WireID_CryostatID;
  std::vector<int> freco_hit_WireID_TPCID;
  std::vector<int> freco_hit_WireID_PlaneID;
  std::vector<int> freco_hit_WireID_WireID;
  //Reco - MC matching
  std::vector<std::vector<int>> freco_hit_mc_type;
  std::vector<std::vector<int>> freco_hit_mc_index;
  /*
  std::vector<std::vector<float>> freco_hit_true_ideFraction;
  std::vector<std::vector<int>> freco_hit_true_isMaxIDE;
  std::vector<std::vector<float>> freco_hit_true_ideNFraction;
  std::vector<std::vector<int>> freco_hit_true_isMaxIDEN;
  */
  std::vector<std::vector<float>> freco_hit_true_numElectrons;
  std::vector<std::vector<float>> freco_hit_true_energy;

  //Reco Track
  std::vector<int> freco_track_producer_index;
  std::vector<int> freco_track_NumberTrajectoryPoints;
  std::vector<int> freco_track_NPoints;
  std::vector<int> freco_track_FirstPoint;
  std::vector<int> freco_track_LastPoint;
  std::vector<int> freco_track_FirstValidPoint;
  std::vector<int> freco_track_LastValidPoint;
  std::vector<int> freco_track_CountValidPoints;
  std::vector<std::vector<double>> freco_track_X;
  std::vector<std::vector<double>> freco_track_Y;
  std::vector<std::vector<double>> freco_track_Z;
  std::vector<std::vector<double>> freco_track_Px;
  std::vector<std::vector<double>> freco_track_Py;
  std::vector<std::vector<double>> freco_track_Pz;
  std::vector<bool> freco_track_HasMomentum;
  std::vector<std::vector<double>> freco_track_Length;
  std::vector<float> freco_track_Chi2;
  std::vector<float> freco_track_Chi2PerNdof;
  std::vector<int> freco_track_Ndof;
  std::vector<int> freco_track_ParticleId;
  std::vector<double> freco_track_Theta;
  std::vector<double> freco_track_Phi;
  std::vector<double> freco_track_ZenithAngle;
  std::vector<double> freco_track_AzimuthAngle;
  /*
  std::vector<std::vector<double>> freco_track_Theta;
  std::vector<std::vector<double>> freco_track_Phi;
  std::vector<std::vector<double>> freco_track_ZenithAngle;
  std::vector<std::vector<double>> freco_track_AzimuthAngle;
  */
  std::vector<double> freco_track_VertexDirection_X;
  std::vector<double> freco_track_VertexDirection_Y;
  std::vector<double> freco_track_VertexDirection_Z;
  std::vector<std::vector<int>> freco_track_to_reco_hit;
  std::vector<std::vector<double>> freco_track_EnergyHelper_resrange;
  std::vector<std::vector<double>> freco_track_EnergyHelper_dedx;
  std::vector<double> freco_track_EnergyHelper_energy;
  std::vector<double> freco_track_EnergyHelperNew_energy_legacy;
  std::vector<std::vector<double>> freco_track_EnergyHelperNew_energy;
  std::vector<std::vector<double>> freco_track_EnergyHelperNew_energy_from_dedx;
  std::vector<std::vector<double>> freco_track_EnergyHelperNew_dedx;

  //Reco - MC matching
  std::vector<int> freco_track_largest_mc_type;
  std::vector<int> freco_track_largest_mc_index;
  std::vector<double> freco_track_largest_ratio;
  std::vector<std::vector<int>> freco_track_mc_type;
  std::vector<std::vector<int>> freco_track_mc_index;
  std::vector<std::vector<double>> freco_track_charge_contribution;
  std::vector<double> freco_track_charge_total;

  //Reco Shower
  std::vector<int> freco_shower_producer_index;
  std::vector<double> freco_shower_Direction_x;
  std::vector<double> freco_shower_Direction_y;
  std::vector<double> freco_shower_Direction_z;
  std::vector<double> freco_shower_DirectionErr_x;
  std::vector<double> freco_shower_DirectionErr_y;
  std::vector<double> freco_shower_DirectionErr_z;
  std::vector<double> freco_shower_ShowerStart_x;
  std::vector<double> freco_shower_ShowerStart_y;
  std::vector<double> freco_shower_ShowerStart_z;
  std::vector<double> freco_shower_ShowerStartErr_x;
  std::vector<double> freco_shower_ShowerStartErr_y;
  std::vector<double> freco_shower_ShowerStartErr_z;
  std::vector<std::vector<double>> freco_shower_Energy;
  std::vector<std::vector<double>> freco_shower_EnergyErr;
  std::vector<std::vector<double>> freco_shower_MIPEnergy;
  std::vector<std::vector<double>> freco_shower_MIPEnergyErr;
  std::vector<int> freco_shower_best_plane;
  std::vector<double> freco_shower_Length;
  std::vector<double> freco_shower_OpenAngle;
  std::vector<std::vector<double>> freco_shower_dEdx;
  std::vector<std::vector<double>> freco_shower_dEdxErr;
  std::vector<bool> freco_shower_has_open_angle;
  std::vector<bool> freco_shower_has_length;
  std::vector<std::vector<int>> freco_shower_to_reco_hit;
  std::vector<double> freco_shower_EnergyHelper_energy_legacy;
  std::vector<std::vector<double>> freco_shower_EnergyHelper_energy;
  std::vector<std::vector<double>> freco_shower_EnergyHelper_dedx;
  std::vector<double> freco_shower_EnergyHelperNew_energy_legacy;
  std::vector<std::vector<double>> freco_shower_EnergyHelperNew_energy;
  std::vector<std::vector<double>> freco_shower_EnergyHelperNew_dedx;
  //Reco - MC matching
  std::vector<int> freco_shower_largest_mc_type;
  std::vector<int> freco_shower_largest_mc_index;
  std::vector<double> freco_shower_largest_ratio;
  std::vector<std::vector<int>> freco_shower_mc_type;
  std::vector<std::vector<int>> freco_shower_mc_index;
  std::vector<std::vector<double>> freco_shower_charge_contribution;
  std::vector<double> freco_shower_charge_total;

  //Pandora
  std::vector<int> fpfp_pdg;
  std::vector<double> fpfp_vertex_X;
  std::vector<double> fpfp_vertex_Y;
  std::vector<double> fpfp_vertex_Z;
  std::vector<int> fpfp_original_index;
  std::vector<std::vector<int>> fpfp_children;

  //Truth
  std::vector<int> fnu_pdg;
  std::vector<double> fnu_energy;
  std::vector<int> flep_pdg;
  std::vector<double> flep_energy;
  std::vector<int> fccnc;
  std::vector<int> fmode;
  std::vector<int> finteraction_type;

  std::vector<double> ftrue_nu_E;
  std::vector<double> ftrue_nuvertx;
  std::vector<double> ftrue_nuverty;
  std::vector<double> ftrue_nuvertz;

  std::vector<int> ftrue_nu_vtx_tpc_contained;
  std::vector<int> ftrue_nu_vtx_fid_contained;

  //GENIE MCParticle
  std::vector<std::vector<int>> fgenie_particle_TrackId;
  std::vector<std::vector<int>> fgenie_particle_StatusCode;
  std::vector<std::vector<int>> fgenie_particle_PdgCode;
  std::vector<std::vector<int>> fgenie_particle_Mother;
  std::vector<std::vector<double>> fgenie_particle_X;
  std::vector<std::vector<double>> fgenie_particle_Y;
  std::vector<std::vector<double>> fgenie_particle_Z;
  std::vector<std::vector<double>> fgenie_particle_T;
  std::vector<std::vector<double>> fgenie_particle_Px;
  std::vector<std::vector<double>> fgenie_particle_Py;
  std::vector<std::vector<double>> fgenie_particle_Pz;
  std::vector<std::vector<double>> fgenie_particle_E;

  //MCParticle
  std::vector<int> fmcparticle_TrackId;
  std::vector<int> fmcparticle_StatusCode;
  std::vector<int> fmcparticle_PdgCode;
  std::vector<int> fmcparticle_Mother;  
  /* 
  std::vector<double> fmcparticle_X;
  std::vector<double> fmcparticle_Y;
  std::vector<double> fmcparticle_Z;
  std::vector<double> fmcparticle_T;
  std::vector<double> fmcparticle_Px;
  std::vector<double> fmcparticle_Py;
  std::vector<double> fmcparticle_Pz;
  std::vector<double> fmcparticle_E;
  */

  //MCTrack
  std::vector<int> fmctrack_Origin;
  std::vector<int> fmctrack_PdgCode;
  std::vector<int> fmctrack_TrackID;
  std::vector<std::string> fmctrack_Process;
  std::vector<double> fmctrack_Start_X;
  std::vector<double> fmctrack_Start_Y;
  std::vector<double> fmctrack_Start_Z;
  std::vector<double> fmctrack_Start_T;
  std::vector<double> fmctrack_Start_Px;
  std::vector<double> fmctrack_Start_Py;
  std::vector<double> fmctrack_Start_Pz;
  std::vector<double> fmctrack_Start_E;
  std::vector<double> fmctrack_End_X;
  std::vector<double> fmctrack_End_Y;
  std::vector<double> fmctrack_End_Z;
  std::vector<double> fmctrack_End_T;
  std::vector<double> fmctrack_End_Px;
  std::vector<double> fmctrack_End_Py;
  std::vector<double> fmctrack_End_Pz;
  std::vector<double> fmctrack_End_E;
  std::vector<std::vector<double>> fmctrack_X;
  std::vector<std::vector<double>> fmctrack_Y;
  std::vector<std::vector<double>> fmctrack_Z;
  std::vector<std::vector<double>> fmctrack_T;
  std::vector<std::vector<double>> fmctrack_Px;
  std::vector<std::vector<double>> fmctrack_Py;
  std::vector<std::vector<double>> fmctrack_Pz;
  std::vector<std::vector<double>> fmctrack_E;
  std::vector<std::vector<std::vector<double>>> fmctrack_dQdx;
  std::vector<std::vector<double>> fmctrack_dEdx;
  std::vector<int> fmctrack_MotherPdgCode;
  std::vector<int> fmctrack_MotherTrackID;
  std::vector<std::string> fmctrack_MotherProcess;
  std::vector<int> fmctrack_AncestorPdgCode;
  std::vector<int> fmctrack_AncestorTrackID;
  std::vector<std::string> fmctrack_AncestorProcess;
  std::vector<std::vector<double>> fmctrack_contributed_charge;

  //MCShower
  std::vector<int> fmcshower_Origin;
  std::vector<int> fmcshower_PdgCode;
  std::vector<int> fmcshower_TrackID;
  std::vector<std::string> fmcshower_Process;
  std::vector<double> fmcshower_Start_X;
  std::vector<double> fmcshower_Start_Y;
  std::vector<double> fmcshower_Start_Z;
  std::vector<double> fmcshower_Start_T;
  std::vector<double> fmcshower_Start_Px;
  std::vector<double> fmcshower_Start_Py;
  std::vector<double> fmcshower_Start_Pz;
  std::vector<double> fmcshower_Start_E;
  std::vector<double> fmcshower_End_X;
  std::vector<double> fmcshower_End_Y;
  std::vector<double> fmcshower_End_Z;
  std::vector<double> fmcshower_End_T;
  std::vector<double> fmcshower_End_Px;
  std::vector<double> fmcshower_End_Py;
  std::vector<double> fmcshower_End_Pz;
  std::vector<double> fmcshower_End_E;
  std::vector<int> fmcshower_MotherPdgCode;
  std::vector<int> fmcshower_MotherTrackID;
  std::vector<std::string> fmcshower_MotherProcess;
  std::vector<int> fmcshower_AncestorPdgCode;
  std::vector<int> fmcshower_AncestorTrackID;
  std::vector<std::string> fmcshower_AncestorProcess;
  std::vector<double> fmcshower_DetProfile_X;
  std::vector<double> fmcshower_DetProfile_Y;
  std::vector<double> fmcshower_DetProfile_Z;
  std::vector<double> fmcshower_DetProfile_T;
  std::vector<double> fmcshower_DetProfile_Px;
  std::vector<double> fmcshower_DetProfile_Py;
  std::vector<double> fmcshower_DetProfile_Pz;
  std::vector<double> fmcshower_DetProfile_E;
  std::vector<std::vector<int>> fmcshower_DaughterTrackID;
  std::vector<std::vector<double>> fmcshower_Charge;
  std::vector<std::vector<double>> fmcshower_dQdx;
  std::vector<double> fmcshower_StartDir_X;
  std::vector<double> fmcshower_StartDir_Y;
  std::vector<double> fmcshower_StartDir_Z;
  std::vector<std::vector<double>> fmcshower_contributed_charge;

  //Delta radiative
  int fdelta_mct_index;
  std::vector<int> fis_delta_rad;
  std::vector<int> fdelta_index;
  std::vector<int> fdelta_photon_index;
  std::vector<int> fdelta_mcshower_index;
  std::vector<int> fdelta_proton_index;
  std::vector<int> fdelta_mctrack_index;

  std::map<std::string, std::vector<double *>> fweight_branch_map;
  
  double fweight_genie_ncelaxial_p1sigma;
  double fweight_genie_ncelaxial_m1sigma;
  double fweight_genie_nceleta_p1sigma;
  double fweight_genie_nceleta_m1sigma;
  double fweight_genie_qema_p1sigma;
  double fweight_genie_qema_m1sigma;
  double fweight_genie_qevec_p1sigma;
  double fweight_genie_qevec_m1sigma;
  double fweight_genie_ccresaxial_p1sigma;
  double fweight_genie_ccresaxial_m1sigma;
  double fweight_genie_ccresvector_p1sigma;
  double fweight_genie_ccresvector_m1sigma;
  double fweight_genie_resganged_p1sigma;
  double fweight_genie_resganged_m1sigma;
  double fweight_genie_ncresaxial_p1sigma;
  double fweight_genie_ncresaxial_m1sigma;
  double fweight_genie_ncresvector_p1sigma;
  double fweight_genie_ncresvector_m1sigma;
  double fweight_genie_cohma_p1sigma;
  double fweight_genie_cohma_m1sigma;
  double fweight_genie_cohr0_p1sigma;
  double fweight_genie_cohr0_m1sigma;
  double fweight_genie_nonresrvp1pi_p1sigma;
  double fweight_genie_nonresrvp1pi_m1sigma;
  double fweight_genie_nonresrvbarp1pi_p1sigma;
  double fweight_genie_nonresrvbarp1pi_m1sigma;
  double fweight_genie_nonresrvp2pi_p1sigma;
  double fweight_genie_nonresrvp2pi_m1sigma;
  double fweight_genie_nonresrvbarp2pi_p1sigma;
  double fweight_genie_nonresrvbarp2pi_m1sigma;
  double fweight_genie_resdecaygamma_p1sigma;
  double fweight_genie_resdecaygamma_m1sigma;
  double fweight_genie_resdecayeta_p1sigma;
  double fweight_genie_resdecayeta_m1sigma;
  double fweight_genie_resdecaytheta_p1sigma;
  double fweight_genie_resdecaytheta_m1sigma;
  double fweight_genie_nc_p1sigma;
  double fweight_genie_nc_m1sigma;
  double fweight_genie_disath_p1sigma;
  double fweight_genie_disath_m1sigma;
  double fweight_genie_disbth_p1sigma;
  double fweight_genie_disbth_m1sigma;
  double fweight_genie_discv1u_p1sigma;
  double fweight_genie_discv1u_m1sigma;
  double fweight_genie_discv2u_p1sigma;
  double fweight_genie_discv2u_m1sigma;
  double fweight_genie_disnucl_p1sigma;
  double fweight_genie_disnucl_m1sigma;
  double fweight_genie_agkyxf_p1sigma;
  double fweight_genie_agkyxf_m1sigma;
  double fweight_genie_agkypt_p1sigma;
  double fweight_genie_agkypt_m1sigma;
  double fweight_genie_formzone_p1sigma;
  double fweight_genie_formzone_m1sigma;
  double fweight_genie_fermigasmodelkf_p1sigma;
  double fweight_genie_fermigasmodelkf_m1sigma;
  double fweight_genie_fermigasmodelsf_p1sigma;
  double fweight_genie_fermigasmodelsf_m1sigma;
  double fweight_genie_intranukenmfp_p1sigma;
  double fweight_genie_intranukenmfp_m1sigma;
  double fweight_genie_intranukencex_p1sigma;
  double fweight_genie_intranukencex_m1sigma;
  double fweight_genie_intranukenel_p1sigma;
  double fweight_genie_intranukenel_m1sigma;
  double fweight_genie_intranukeninel_p1sigma;
  double fweight_genie_intranukeninel_m1sigma;
  double fweight_genie_intranukenabs_p1sigma;
  double fweight_genie_intranukenabs_m1sigma;
  double fweight_genie_intranukenpi_p1sigma;
  double fweight_genie_intranukenpi_m1sigma;
  double fweight_genie_intranukepimfp_p1sigma;
  double fweight_genie_intranukepimfp_m1sigma;
  double fweight_genie_intranukepicex_p1sigma;
  double fweight_genie_intranukepicex_m1sigma;
  double fweight_genie_intranukepiel_p1sigma;
  double fweight_genie_intranukepiel_m1sigma;
  double fweight_genie_intranukepiinel_p1sigma;
  double fweight_genie_intranukepiinel_m1sigma;
  double fweight_genie_intranukepiabs_p1sigma;
  double fweight_genie_intranukepiabs_m1sigma;
  double fweight_genie_intranukepipi_p1sigma;
  double fweight_genie_intranukepipi_m1sigma;

};


FillLightEvent::FillLightEvent(fhicl::ParameterSet const & p) :
  EDAnalyzer(p),
  fheavy(false),
  fuse_eventweight(true),
  ftpc_volume(0,
	      -lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
	      0,
	      2*lar::providerFrom<geo::Geometry>()->DetHalfWidth(),
	      lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
	      lar::providerFrom<geo::Geometry>()->DetLength()),
  foffset(10),
  ffiducial_volume(foffset,
		   foffset-lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
		   foffset,
		   -foffset+2*lar::providerFrom<geo::Geometry>()->DetHalfWidth(),
		   -foffset+lar::providerFrom<geo::Geometry>()->DetHalfHeight(),
		   -foffset+lar::providerFrom<geo::Geometry>()->DetLength()),
  fverbose(false),
  fpot_tree(nullptr),
  fnumber_of_events(0),
  fevent_tree(nullptr)
{

  Reconfigure(p);
  SetupTrees();
  
}


void FillLightEvent::Reconfigure(fhicl::ParameterSet const & p) {

  p.get_if_present<bool>("heavy", fheavy);
  p.get_if_present<bool>("use_eventweight",fuse_eventweight);
  fmc = p.get<bool>("mc");

  p.get_if_present<std::string>("pot_producer", fpot_producer);

  fswtrigger_product = p.get<std::string>("trigger_product");
  fopflash_producers = p.get<std::vector<std::string>>("opflash_producers");

  fhit_producers = p.get<std::vector<std::string>>("hit_producers");
  ftrack_producers = p.get<std::vector<std::string>>("track_producers");
  ftrack_calo_producer = p.get<std::string>("track_calo_producer");
  fshower_producers = p.get<std::vector<std::string>>("shower_producers");

  p.get_if_present<std::vector<std::string>>("rmcmassociation_producers", frmcmassociation_producers);

  fopflashp_size = fopflash_producers.size();
  fhitp_size = fhit_producers.size();
  ftrackp_size = ftrack_producers.size();
  fshowerp_size = fshower_producers.size();
  frmcm_size = frmcmassociation_producers.size();
  if(fmc && !frmcmassociation_producers.empty()) {
    if(frmcm_size != fhitp_size ||
       frmcm_size != ftrackp_size ||
       frmcm_size != fshowerp_size) {
      std::cout << "Must have same number of producers\n";
      exit(1);
    }
    frmcm.reserve(frmcm_size);
    for(size_t i = 0; i < frmcm_size; ++i) {
      RecoMCMatching rmcm;
      rmcm.Configure(fhit_producers.at(i),
		     ftrack_producers.at(i),
		     fshower_producers.at(i),
		     frmcmassociation_producers.at(i));
      frmcm.push_back(rmcm);
    }
    if(frmcm_size) frmcm_first = &frmcm.front();
    else frmcm_first = nullptr;
    fmc_type_shower = frmcm.front().fmc_type_shower;
    fmc_type_track = frmcm.front().fmc_type_track;
    fmc_type_particle = frmcm.front().fmc_type_particle;
  }

}


void FillLightEvent::SetupTrees() {
  
  art::ServiceHandle< art::TFileService > tfs;

  //Change producer tree to be per event and contain a map of producer - object indices

  TTree * meta_tree = tfs->make<TTree>("meta_tree", "");
  int fis_heavy;
  if(fheavy) fis_heavy = 1;
  else fis_heavy = 0;
  meta_tree->Branch("is_heavy", &fis_heavy, "is_heavy/I");
  int fis_mc;
  if(fmc) fis_mc = 1;
  else fis_mc = 0;
  meta_tree->Branch("is_mc", &fis_mc, "is_mc/I");
  meta_tree->Branch("swtrigger_product", &fswtrigger_product);
  meta_tree->Branch("opflash_producers", &fopflash_producers);
  meta_tree->Branch("hit_producers", &fhit_producers);
  meta_tree->Branch("track_producers", &ftrack_producers);
  meta_tree->Branch("shower_producers", &fshower_producers);
  meta_tree->Branch("rmcmassociation_producers", &frmcmassociation_producers);
  double DetHalfHeight = lar::providerFrom<geo::Geometry>()->DetHalfHeight();
  double DetHalfWidth = lar::providerFrom<geo::Geometry>()->DetHalfWidth();
  double DetLength = lar::providerFrom<geo::Geometry>()->DetLength();
  meta_tree->Branch("DetHalfHeight", &DetHalfHeight);
  meta_tree->Branch("DetHalfWidth", &DetHalfWidth);
  meta_tree->Branch("DetLength", &DetLength);
  int mc_type_shower = -1;
  int mc_type_track = -1;
  int mc_type_particle = -1;
  if(fmc && !frmcmassociation_producers.empty()) {
    mc_type_shower = frmcm_first->fmc_type_shower;
    mc_type_track = frmcm_first->fmc_type_track;
    mc_type_particle = frmcm_first->fmc_type_particle;
  }
  meta_tree->Branch("mc_type_shower", &mc_type_shower);
  meta_tree->Branch("mc_type_track", &mc_type_track);
  meta_tree->Branch("mc_type_particle", &mc_type_particle);
  meta_tree->Fill();

  if(fpot_producer != "") {
    fpot_tree = tfs->make<TTree>("pot_tree", "");
    fpot_tree->Branch("number_of_events", &fnumber_of_events, "number_of_events/I");
    fpot_tree->Branch("pot", &fpot, "pot/D");
  }

  fevent_tree = tfs->make<TTree>("event_tree", "");

  fevent_tree->Branch("opflash_producer_indices", &fopflash_producer_indices);
  fevent_tree->Branch("hit_producer_indices", &fhit_producer_indices);
  fevent_tree->Branch("track_producer_indices", &ftrack_producer_indices);
  fevent_tree->Branch("shower_producer_indices", &fshower_producer_indices);
  
  fevent_tree->Branch("run_number", &frun_number, "run_number/I");
  fevent_tree->Branch("subrun_number", &fsubrun_number, "subrun_number/I"); 
  fevent_tree->Branch("event_number", &fevent_number, "event_number/I");

  //fevent_tree->Branch("swtrigger_producer_index", &fswtrigger_producer_index);
  fevent_tree->Branch("swtrigger_getListOfAlgorithms", &fswtrigger_getListOfAlgorithms);
  fevent_tree->Branch("swtrigger_passedAlgo", &fswtrigger_passedAlgo);
  fevent_tree->Branch("swtrigger_passedPrescaleAlgo", &fswtrigger_passedPrescaleAlgo);
  fevent_tree->Branch("swtrigger_vetoAlgo", &fswtrigger_vetoAlgo);
  fevent_tree->Branch("swtrigger_passedAlgos", &fswtrigger_passedAlgos, "swtrigger_passedAlgos/I");
  fevent_tree->Branch("swtrigger_vetoAlgos", &fswtrigger_vetoAlgos, "swtrigger_vetoAlgos/I");
  fevent_tree->Branch("swtrigger_passedPrescaleAlgos", &fswtrigger_passedPrescaleAlgos, "swtrigger_passedPrescaleAlgos/I");
  fevent_tree->Branch("swtrigger_passed_swtrigger", &fswtrigger_passed_swtrigger, "swtrigger_passed_swtrigger/I");

  fevent_tree->Branch("reco_opflash_producer_index", &freco_opflash_producer_index);
  fevent_tree->Branch("reco_opflash_Time", &freco_opflash_Time);
  if(fheavy) {
    fevent_tree->Branch("reco_opflash_TimeWidth", &freco_opflash_TimeWidth);
    fevent_tree->Branch("reco_opflash_AbsTime", &freco_opflash_AbsTime);
    fevent_tree->Branch("reco_opflash_Frame", &freco_opflash_Frame);
    //fevent_tree->Branch("reco_opflash_PEs", &freco_opflash_PEs);
  }
  fevent_tree->Branch("reco_opflash_YCenter", &freco_opflash_YCenter);
  fevent_tree->Branch("reco_opflash_YWidth", &freco_opflash_YWidth);
  fevent_tree->Branch("reco_opflash_ZCenter", &freco_opflash_ZCenter);
  fevent_tree->Branch("reco_opflash_ZWidth", &freco_opflash_ZWidth);
  fevent_tree->Branch("reco_opflash_InBeamFrame", &freco_opflash_InBeamFrame);
  fevent_tree->Branch("reco_opflash_OnBeamFrame", &freco_opflash_OnBeamFrame);
  fevent_tree->Branch("reco_opflash_WireCenters", &freco_opflash_WireCenters);
  fevent_tree->Branch("reco_opflash_WireWidths", &freco_opflash_WireWidths);
  fevent_tree->Branch("reco_opflash_TotalPE", &freco_opflash_TotalPE);
  if(fheavy) fevent_tree->Branch("reco_opflash_FastToTotal", &freco_opflash_FastToTotal);

  if(fheavy) {
    fevent_tree->Branch("reco_hit_producer_index", &freco_hit_producer_index);
    fevent_tree->Branch("reco_hit_StartTick", &freco_hit_StartTick);
    fevent_tree->Branch("reco_hit_EndTick", &freco_hit_EndTick);
    fevent_tree->Branch("reco_hit_PeakTime", &freco_hit_PeakTime);
    if(fheavy) {
      fevent_tree->Branch("reco_hit_SigmaPeakTime", &freco_hit_SigmaPeakTime);
      fevent_tree->Branch("reco_hit_RMS", &freco_hit_RMS);
      fevent_tree->Branch("reco_hit_PeakAmplitude", &freco_hit_PeakAmplitude);
      fevent_tree->Branch("reco_hit_SigmaPeakAmplitude", &freco_hit_SigmaPeakAmplitude);
    }
    fevent_tree->Branch("reco_hit_SummedADC", &freco_hit_SummedADC);
    fevent_tree->Branch("reco_hit_Integral", &freco_hit_Integral);
    if(fheavy) {  
      fevent_tree->Branch("reco_hit_SigmaIntegral", &freco_hit_SigmaIntegral);
      fevent_tree->Branch("reco_hit_Multiplicity", &freco_hit_Multiplicity);
      fevent_tree->Branch("reco_hit_LocalIndex", &freco_hit_LocalIndex);
      fevent_tree->Branch("reco_hit_GoodnessOfFit", &freco_hit_GoodnessOfFit);
      fevent_tree->Branch("reco_hit_DegreesOfFreedom", &freco_hit_DegreesOfFreedom);
    }
    fevent_tree->Branch("reco_hit_View", &freco_hit_View);
    fevent_tree->Branch("reco_hit_SignalType", &freco_hit_SignalType);
    fevent_tree->Branch("reco_hit_WireID_CryostatID", &freco_hit_WireID_CryostatID);
    fevent_tree->Branch("reco_hit_WireID_TPCID", &freco_hit_WireID_TPCID);
    fevent_tree->Branch("reco_hit_WireID_PlaneID", &freco_hit_WireID_PlaneID);
    fevent_tree->Branch("reco_hit_WireID_WireID", &freco_hit_WireID_WireID);
    if(frmcm_first) {
      fevent_tree->Branch("reco_hit_mc_type", &freco_hit_mc_type);
      fevent_tree->Branch("reco_hit_mc_index", &freco_hit_mc_index);
      /*
	fevent_tree->Branch("reco_hit_true_ideFraction", &freco_hit_true_ideFraction);
	fevent_tree->Branch("reco_hit_true_isMaxIDE", &freco_hit_true_isMaxIDE);
	fevent_tree->Branch("reco_hit_true_ideNFraction", &freco_hit_true_ideNFraction);
	fevent_tree->Branch("reco_hit_true_isMaxIDEN", &freco_hit_true_isMaxIDEN);
      */
      fevent_tree->Branch("reco_hit_true_numElectrons", &freco_hit_true_numElectrons);
      if(fheavy) fevent_tree->Branch("reco_hit_true_energy", &freco_hit_true_energy);
    }
  }
  fevent_tree->Branch("reco_track_producer_index", &freco_track_producer_index);
  if(fheavy) {
    fevent_tree->Branch("reco_track_NumberTrajectoryPoints", &freco_track_NumberTrajectoryPoints);
    fevent_tree->Branch("reco_track_NPoints", &freco_track_NPoints);
    fevent_tree->Branch("reco_track_FirstPoint", &freco_track_FirstPoint);
    fevent_tree->Branch("reco_track_LastPoint", &freco_track_LastPoint);
    fevent_tree->Branch("reco_track_FirstValidPoint", &freco_track_FirstValidPoint);
    fevent_tree->Branch("reco_track_LastValidPoint", &freco_track_LastValidPoint);
    fevent_tree->Branch("reco_track_CountValidPoints", &freco_track_CountValidPoints);
  }
  fevent_tree->Branch("reco_track_X", &freco_track_X);
  fevent_tree->Branch("reco_track_Y", &freco_track_Y);
  fevent_tree->Branch("reco_track_Z", &freco_track_Z);
  if(fheavy) {
    fevent_tree->Branch("reco_track_Px", &freco_track_Px);
    fevent_tree->Branch("reco_track_Py", &freco_track_Py);
    fevent_tree->Branch("reco_track_Pz", &freco_track_Pz);
    fevent_tree->Branch("reco_track_HasMomentum", &freco_track_HasMomentum);
    fevent_tree->Branch("reco_track_Length", &freco_track_Length);
    fevent_tree->Branch("reco_track_Chi2", &freco_track_Chi2);
    fevent_tree->Branch("reco_track_Chi2PerNdof", &freco_track_Chi2PerNdof);
    fevent_tree->Branch("reco_track_Ndof", &freco_track_Ndof);
    fevent_tree->Branch("reco_track_ParticleId", &freco_track_ParticleId);
  }
  fevent_tree->Branch("reco_track_Theta", &freco_track_Theta);
  fevent_tree->Branch("reco_track_Phi", &freco_track_Phi);
  fevent_tree->Branch("reco_track_ZenithAngle", &freco_track_ZenithAngle);
  fevent_tree->Branch("reco_track_AzimuthAngle", &freco_track_AzimuthAngle);
  fevent_tree->Branch("reco_track_VertexDirection_X", &freco_track_VertexDirection_X);
  fevent_tree->Branch("reco_track_VertexDirection_Y", &freco_track_VertexDirection_Y);
  fevent_tree->Branch("reco_track_VertexDirection_Z", &freco_track_VertexDirection_Z);
  fevent_tree->Branch("reco_track_to_reco_hit", &freco_track_to_reco_hit);
  fevent_tree->Branch("reco_track_EnergyHelper_resrange", &freco_track_EnergyHelper_resrange);
  fevent_tree->Branch("reco_track_EnergyHelper_dedx", &freco_track_EnergyHelper_dedx);
  fevent_tree->Branch("reco_track_EnergyHelper_energy", &freco_track_EnergyHelper_energy);
  fevent_tree->Branch("reco_track_EnergyHelperNew_energy_legacy", &freco_track_EnergyHelperNew_energy_legacy);
  fevent_tree->Branch("reco_track_EnergyHelperNew_energy", &freco_track_EnergyHelperNew_energy);
  fevent_tree->Branch("reco_track_EnergyHelperNew_energy_from_dedx", &freco_track_EnergyHelperNew_energy_from_dedx);
  fevent_tree->Branch("reco_track_EnergyHelperNew_dedx", &freco_track_EnergyHelperNew_dedx);
  if(frmcm_first) {
    fevent_tree->Branch("reco_track_largest_mc_type", &freco_track_largest_mc_type);
    fevent_tree->Branch("reco_track_largest_mc_index", &freco_track_largest_mc_index);
    fevent_tree->Branch("reco_track_largest_ratio", &freco_track_largest_ratio);
    //if(fheavy) {
      fevent_tree->Branch("reco_track_mc_type", &freco_track_mc_type);
      fevent_tree->Branch("reco_track_mc_index", &freco_track_mc_index);
      fevent_tree->Branch("reco_track_charge_contribution", &freco_track_charge_contribution);
      //}
    fevent_tree->Branch("reco_track_charge_total", &freco_track_charge_total);
  }

  fevent_tree->Branch("reco_shower_producer_index", &freco_shower_producer_index);
  fevent_tree->Branch("reco_shower_Direction_x", &freco_shower_Direction_x);
  fevent_tree->Branch("reco_shower_Direction_y", &freco_shower_Direction_y);
  fevent_tree->Branch("reco_shower_Direction_z", &freco_shower_Direction_z);
  if(fheavy) {
    fevent_tree->Branch("reco_shower_DirectionErr_x", &freco_shower_DirectionErr_x);
    fevent_tree->Branch("reco_shower_DirectionErr_y", &freco_shower_DirectionErr_y);
    fevent_tree->Branch("reco_shower_DirectionErr_z", &freco_shower_DirectionErr_z);
  }
  fevent_tree->Branch("reco_shower_ShowerStart_x", &freco_shower_ShowerStart_x);
  fevent_tree->Branch("reco_shower_ShowerStart_y", &freco_shower_ShowerStart_y);
  fevent_tree->Branch("reco_shower_ShowerStart_z", &freco_shower_ShowerStart_z);
  if(fheavy) {
    fevent_tree->Branch("reco_shower_ShowerStartErr_x", &freco_shower_ShowerStartErr_x);
    fevent_tree->Branch("reco_shower_ShowerStartErr_y", &freco_shower_ShowerStartErr_y);
    fevent_tree->Branch("reco_shower_ShowerStartErr_z", &freco_shower_ShowerStartErr_z);
    fevent_tree->Branch("reco_shower_Energy", &freco_shower_Energy);
    fevent_tree->Branch("reco_shower_EnergyErr", &freco_shower_EnergyErr);
    fevent_tree->Branch("reco_shower_MIPEnergy", &freco_shower_MIPEnergy);
    fevent_tree->Branch("reco_shower_MIPEnergyErr", &freco_shower_MIPEnergyErr);
  }
  fevent_tree->Branch("reco_shower_best_plane", &freco_shower_best_plane);
  fevent_tree->Branch("reco_shower_Length", &freco_shower_Length);
  fevent_tree->Branch("reco_shower_OpenAngle", &freco_shower_OpenAngle);
  if(fheavy) {
    fevent_tree->Branch("reco_shower_dEdx", &freco_shower_dEdx);
    fevent_tree->Branch("reco_shower_dEdxErr", &freco_shower_dEdxErr);
    fevent_tree->Branch("reco_shower_has_open_angle", &freco_shower_has_open_angle);
    fevent_tree->Branch("reco_shower_has_length", &freco_shower_has_length);
  }
  fevent_tree->Branch("reco_shower_to_reco_hit", &freco_shower_to_reco_hit);
  fevent_tree->Branch("reco_shower_EnergyHelper_energy_legacy", &freco_shower_EnergyHelper_energy_legacy);
  fevent_tree->Branch("reco_shower_EnergyHelper_energy", &freco_shower_EnergyHelper_energy);
  fevent_tree->Branch("reco_shower_EnergyHelper_dedx", &freco_shower_EnergyHelper_dedx);
  fevent_tree->Branch("reco_shower_EnergyHelperNew_energy_legacy", &freco_shower_EnergyHelperNew_energy_legacy);
  fevent_tree->Branch("reco_shower_EnergyHelperNew_energy", &freco_shower_EnergyHelperNew_energy);
  fevent_tree->Branch("reco_shower_EnergyHelperNew_dedx", &freco_shower_EnergyHelperNew_dedx);
  if(frmcm_first) { 
    fevent_tree->Branch("reco_shower_largest_mc_type", &freco_shower_largest_mc_type);
    fevent_tree->Branch("reco_shower_largest_mc_index", &freco_shower_largest_mc_index);
    fevent_tree->Branch("reco_shower_largest_ratio", &freco_shower_largest_ratio);
    //if(fheavy) {
      fevent_tree->Branch("reco_shower_mc_type", &freco_shower_mc_type);
      fevent_tree->Branch("reco_shower_mc_index", &freco_shower_mc_index);
      fevent_tree->Branch("reco_shower_charge_contribution", &freco_shower_charge_contribution);
      //}
    fevent_tree->Branch("reco_shower_charge_total", &freco_shower_charge_total);
  }

  fevent_tree->Branch("pfp_pdg", &fpfp_pdg);
  fevent_tree->Branch("pfp_vertex_X", &fpfp_vertex_X);
  fevent_tree->Branch("pfp_vertex_Y", &fpfp_vertex_Y);
  fevent_tree->Branch("pfp_vertex_Z", &fpfp_vertex_Z);
  fevent_tree->Branch("pfp_original_index", &fpfp_original_index);
  fevent_tree->Branch("pfp_children", &fpfp_children);  

  if(fmc) {

    fevent_tree->Branch("nu_pdg", &fnu_pdg);
    fevent_tree->Branch("nu_energy", &fnu_energy);
    fevent_tree->Branch("lep_pdg", &flep_pdg);
    fevent_tree->Branch("lep_energy", &flep_energy);
    fevent_tree->Branch("ccnc", &fccnc);
    fevent_tree->Branch("mode", &fmode);
    fevent_tree->Branch("interaction_type", &finteraction_type);

    fevent_tree->Branch("true_nu_E", &ftrue_nu_E);

    fevent_tree->Branch("true_nuvertx", &ftrue_nuvertx);
    fevent_tree->Branch("true_nuverty", &ftrue_nuverty);
    fevent_tree->Branch("true_nuvertz", &ftrue_nuvertz);

    fevent_tree->Branch("true_nu_vtx_tpc_contained", &ftrue_nu_vtx_tpc_contained); 
    fevent_tree->Branch("true_nu_vtx_fid_contained", &ftrue_nu_vtx_fid_contained); 

    fevent_tree->Branch("genie_particle_TrackId", &fgenie_particle_TrackId);
    fevent_tree->Branch("genie_particle_StatusCode", &fgenie_particle_StatusCode);
    fevent_tree->Branch("genie_particle_PdgCode", &fgenie_particle_PdgCode);
    fevent_tree->Branch("genie_particle_Mother", &fgenie_particle_Mother);
    fevent_tree->Branch("genie_particle_X", &fgenie_particle_X);
    fevent_tree->Branch("genie_particle_Y", &fgenie_particle_Y);
    fevent_tree->Branch("genie_particle_Z", &fgenie_particle_Z);
    fevent_tree->Branch("genie_particle_Px", &fgenie_particle_Px);
    fevent_tree->Branch("genie_particle_Py", &fgenie_particle_Py);
    fevent_tree->Branch("genie_particle_Pz", &fgenie_particle_Pz);
    fevent_tree->Branch("genie_particle_E", &fgenie_particle_E);

    fevent_tree->Branch("mcparticle_TrackId", &fmcparticle_TrackId);
    fevent_tree->Branch("mcparticle_StatusCode", &fmcparticle_StatusCode);
    fevent_tree->Branch("mcparticle_PdgCode", &fmcparticle_PdgCode);
    fevent_tree->Branch("mcparticle_Mother", &fmcparticle_Mother);
    /*
    fevent_tree->Branch("mcparticle_X", &fmcparticle_X);
    fevent_tree->Branch("mcparticle_Y", &fmcparticle_Y);
    fevent_tree->Branch("mcparticle_Z", &fmcparticle_Z);
    fevent_tree->Branch("mcparticle_Px", &fmcparticle_Px);
    fevent_tree->Branch("mcparticle_Py", &fmcparticle_Py);
    fevent_tree->Branch("mcparticle_Pz", &fmcparticle_Pz);
    fevent_tree->Branch("mcparticle_E", &fmcparticle_E);
    */

    fevent_tree->Branch("mctrack_Origin", &fmctrack_Origin);
    fevent_tree->Branch("mctrack_PdgCode", &fmctrack_PdgCode);
    fevent_tree->Branch("mctrack_TrackID", &fmctrack_TrackID);
    fevent_tree->Branch("mctrack_Process", &fmctrack_Process);
    fevent_tree->Branch("mctrack_Start_X", &fmctrack_Start_X);
    fevent_tree->Branch("mctrack_Start_Y", &fmctrack_Start_Y);
    fevent_tree->Branch("mctrack_Start_Z", &fmctrack_Start_Z);
    fevent_tree->Branch("mctrack_Start_T", &fmctrack_Start_T);
    fevent_tree->Branch("mctrack_Start_Px", &fmctrack_Start_Px);
    fevent_tree->Branch("mctrack_Start_Py", &fmctrack_Start_Py);
    fevent_tree->Branch("mctrack_Start_Pz", &fmctrack_Start_Pz);
    fevent_tree->Branch("mctrack_Start_E", &fmctrack_Start_E);
    fevent_tree->Branch("mctrack_End_X", &fmctrack_End_X);
    fevent_tree->Branch("mctrack_End_Y", &fmctrack_End_Y);
    fevent_tree->Branch("mctrack_End_Z", &fmctrack_End_Z);
    fevent_tree->Branch("mctrack_End_T", &fmctrack_End_T);
    fevent_tree->Branch("mctrack_End_Px", &fmctrack_End_Px);
    fevent_tree->Branch("mctrack_End_Py", &fmctrack_End_Py);
    fevent_tree->Branch("mctrack_End_Pz", &fmctrack_End_Pz);
    fevent_tree->Branch("mctrack_End_E", &fmctrack_End_E);
    fevent_tree->Branch("mctrack_X", &fmctrack_X);
    fevent_tree->Branch("mctrack_Y", &fmctrack_Y);
    fevent_tree->Branch("mctrack_Z", &fmctrack_Z);
    fevent_tree->Branch("mctrack_T", &fmctrack_T);
    if(fheavy) {
      fevent_tree->Branch("mctrack_Px", &fmctrack_Px);
      fevent_tree->Branch("mctrack_Py", &fmctrack_Py);
      fevent_tree->Branch("mctrack_Pz", &fmctrack_Pz);
      fevent_tree->Branch("mctrack_E", &fmctrack_E);
      fevent_tree->Branch("mctrack_dQdx", &fmctrack_dQdx);
      fevent_tree->Branch("mctrack_dEdx", &fmctrack_dEdx);
    }
    fevent_tree->Branch("mctrack_MotherTrackID", &fmctrack_MotherTrackID);
    fevent_tree->Branch("mctrack_MotherPdgCode", &fmctrack_MotherPdgCode);
    fevent_tree->Branch("mctrack_MotherProcess", &fmctrack_MotherProcess);
    fevent_tree->Branch("mctrack_AncestorTrackID", &fmctrack_AncestorTrackID);
    fevent_tree->Branch("mctrack_AncestorPdgCode", &fmctrack_AncestorPdgCode);
    fevent_tree->Branch("mctrack_AncestorProcess", &fmctrack_AncestorProcess);
    if(frmcm_first) fevent_tree->Branch("mctrack_contributed_charge", &fmctrack_contributed_charge);

    fevent_tree->Branch("mcshower_Origin", &fmcshower_Origin);
    fevent_tree->Branch("mcshower_PdgCode", &fmcshower_PdgCode);
    fevent_tree->Branch("mcshower_TrackID", &fmcshower_TrackID);
    fevent_tree->Branch("mcshower_Process", &fmcshower_Process);
    fevent_tree->Branch("mcshower_Start_X", &fmcshower_Start_X);
    fevent_tree->Branch("mcshower_Start_Y", &fmcshower_Start_Y);
    fevent_tree->Branch("mcshower_Start_Z", &fmcshower_Start_Z);
    fevent_tree->Branch("mcshower_Start_T", &fmcshower_Start_T);
    fevent_tree->Branch("mcshower_Start_Px", &fmcshower_Start_Px);
    fevent_tree->Branch("mcshower_Start_Py", &fmcshower_Start_Py);
    fevent_tree->Branch("mcshower_Start_Pz", &fmcshower_Start_Pz);
    fevent_tree->Branch("mcshower_Start_E", &fmcshower_Start_E);
    fevent_tree->Branch("mcshower_End_X", &fmcshower_End_X);
    fevent_tree->Branch("mcshower_End_Y", &fmcshower_End_Y);
    fevent_tree->Branch("mcshower_End_Z", &fmcshower_End_Z);
    fevent_tree->Branch("mcshower_End_T", &fmcshower_End_T);
    fevent_tree->Branch("mcshower_End_Px", &fmcshower_End_Px);
    fevent_tree->Branch("mcshower_End_Py", &fmcshower_End_Py);
    fevent_tree->Branch("mcshower_End_Pz", &fmcshower_End_Pz);
    fevent_tree->Branch("mcshower_End_E", &fmcshower_End_E);
    fevent_tree->Branch("mcshower_MotherTrackID", &fmcshower_MotherTrackID);
    fevent_tree->Branch("mcshower_MotherPdgCode", &fmcshower_MotherPdgCode);
    fevent_tree->Branch("mcshower_MotherProcess", &fmcshower_MotherProcess);
    fevent_tree->Branch("mcshower_AncestorTrackID", &fmcshower_AncestorTrackID);
    fevent_tree->Branch("mcshower_AncestorPdgCode", &fmcshower_AncestorPdgCode);
    fevent_tree->Branch("mcshower_AncestorProcess", &fmcshower_AncestorProcess);
    fevent_tree->Branch("mcshower_DetProfile_X", &fmcshower_DetProfile_X);
    fevent_tree->Branch("mcshower_DetProfile_Y", &fmcshower_DetProfile_Y);
    fevent_tree->Branch("mcshower_DetProfile_Z", &fmcshower_DetProfile_Z);
    fevent_tree->Branch("mcshower_DetProfile_T", &fmcshower_DetProfile_T);
    fevent_tree->Branch("mcshower_DetProfile_Px", &fmcshower_DetProfile_Px);
    fevent_tree->Branch("mcshower_DetProfile_Py", &fmcshower_DetProfile_Py);
    fevent_tree->Branch("mcshower_DetProfile_Pz", &fmcshower_DetProfile_Pz);
    fevent_tree->Branch("mcshower_DetProfile_E", &fmcshower_DetProfile_E);
    fevent_tree->Branch("mcshower_DaughterTrackID", &fmcshower_DaughterTrackID);
    if(fheavy) {
      fevent_tree->Branch("mcshower_Charge", &fmcshower_Charge);
      fevent_tree->Branch("mcshower_dQdx", &fmcshower_dQdx);
    }
    fevent_tree->Branch("mcshower_StartDir_X", &fmcshower_StartDir_X);
    fevent_tree->Branch("mcshower_StartDir_Y", &fmcshower_StartDir_Y);
    fevent_tree->Branch("mcshower_StartDir_Z", &fmcshower_StartDir_Z);
    if(frmcm_first) fevent_tree->Branch("mcshower_contributed_charge", &fmcshower_contributed_charge);

    fevent_tree->Branch("delta_mct_index", &fdelta_mct_index);
    fevent_tree->Branch("is_delta_rad", &fis_delta_rad);
    fevent_tree->Branch("delta_index", &fdelta_index);
    fevent_tree->Branch("delta_photon_index", &fdelta_photon_index);
    fevent_tree->Branch("delta_mcshower_index", &fdelta_mcshower_index);
    fevent_tree->Branch("delta_proton_index", &fdelta_proton_index);
    fevent_tree->Branch("delta_mctrack_index", &fdelta_mctrack_index);

    if(fheavy) {

      fevent_tree->Branch("fweight_genie_ncelaxial_p1sigma", &fweight_genie_ncelaxial_p1sigma, "fweight_genie_ncelaxial_p1sigma/D");
      fevent_tree->Branch("fweight_genie_ncelaxial_m1sigma", &fweight_genie_ncelaxial_m1sigma, "fweight_genie_ncelaxial_m1sigma/D");
      fweight_branch_map.emplace("genie_NCELaxial_Genie", std::vector<double *>{&fweight_genie_ncelaxial_p1sigma, &fweight_genie_ncelaxial_m1sigma});

      fevent_tree->Branch("fweight_genie_nceleta_p1sigma", &fweight_genie_nceleta_p1sigma, "fweight_genie_nceleta_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nceleta_m1sigma", &fweight_genie_nceleta_m1sigma, "fweight_genie_nceleta_m1sigma/D");
      fweight_branch_map.emplace("genie_NCELeta_Genie", std::vector<double *>{&fweight_genie_nceleta_p1sigma, &fweight_genie_nceleta_m1sigma});

      fevent_tree->Branch("fweight_genie_qema_p1sigma", &fweight_genie_qema_p1sigma, "fweight_genie_qema_p1sigma/D");
      fevent_tree->Branch("fweight_genie_qema_m1sigma", &fweight_genie_qema_m1sigma, "fweight_genie_qema_m1sigma/D");
      fweight_branch_map.emplace("genie_QEMA_Genie", std::vector<double *>{&fweight_genie_qema_p1sigma, &fweight_genie_qema_m1sigma});

      fevent_tree->Branch("fweight_genie_qevec_p1sigma", &fweight_genie_qevec_p1sigma, "fweight_genie_qevec_p1sigma/D");
      fevent_tree->Branch("fweight_genie_qevec_m1sigma", &fweight_genie_qevec_m1sigma, "fweight_genie_qevec_m1sigma/D");
      fweight_branch_map.emplace("genie_QEVec_Genie", std::vector<double *>{&fweight_genie_qevec_p1sigma, &fweight_genie_qevec_m1sigma});

      fevent_tree->Branch("fweight_genie_ccresaxial_p1sigma", &fweight_genie_ccresaxial_p1sigma, "fweight_genie_ccresaxial_p1sigma/D");
      fevent_tree->Branch("fweight_genie_ccresaxial_m1sigma", &fweight_genie_ccresaxial_m1sigma, "fweight_genie_ccresaxial_m1sigma/D");
      fweight_branch_map.emplace("genie_CCResAxial_Genie", std::vector<double *>{&fweight_genie_ccresaxial_p1sigma, &fweight_genie_ccresaxial_m1sigma});

      fevent_tree->Branch("fweight_genie_ccresvector_p1sigma", &fweight_genie_ccresvector_p1sigma, "fweight_genie_ccresvector_p1sigma/D");
      fevent_tree->Branch("fweight_genie_ccresvector_m1sigma", &fweight_genie_ccresvector_m1sigma, "fweight_genie_ccresvector_m1sigma/D");
      fweight_branch_map.emplace("genie_CCResVector_Genie", std::vector<double *>{&fweight_genie_ccresvector_p1sigma, &fweight_genie_ccresvector_m1sigma});

      fevent_tree->Branch("fweight_genie_resganged_p1sigma", &fweight_genie_resganged_p1sigma, "fweight_genie_resganged_p1sigma/D");
      fevent_tree->Branch("fweight_genie_resganged_m1sigma", &fweight_genie_resganged_m1sigma, "fweight_genie_resganged_m1sigma/D");
      fweight_branch_map.emplace("genie_ResGanged_Genie", std::vector<double *>{&fweight_genie_resganged_p1sigma, &fweight_genie_resganged_m1sigma});

      fevent_tree->Branch("fweight_genie_ncresaxial_p1sigma", &fweight_genie_ncresaxial_p1sigma, "fweight_genie_ncresaxial_p1sigma/D");
      fevent_tree->Branch("fweight_genie_ncresaxial_m1sigma", &fweight_genie_ncresaxial_m1sigma, "fweight_genie_ncresaxial_m1sigma/D");
      fweight_branch_map.emplace("genie_NCResAxial_Genie", std::vector<double *>{&fweight_genie_ncresaxial_p1sigma, &fweight_genie_ncresaxial_m1sigma});

      fevent_tree->Branch("fweight_genie_ncresvector_p1sigma", &fweight_genie_ncresvector_p1sigma, "fweight_genie_ncresvector_p1sigma/D");
      fevent_tree->Branch("fweight_genie_ncresvector_m1sigma", &fweight_genie_ncresvector_m1sigma, "fweight_genie_ncresvector_m1sigma/D");
      fweight_branch_map.emplace("genie_NCResVector_Genie", std::vector<double *>{&fweight_genie_ncresvector_p1sigma, &fweight_genie_ncresvector_m1sigma});

      fevent_tree->Branch("fweight_genie_cohma_p1sigma", &fweight_genie_cohma_p1sigma, "fweight_genie_cohma_p1sigma/D");
      fevent_tree->Branch("fweight_genie_cohma_m1sigma", &fweight_genie_cohma_m1sigma, "fweight_genie_cohma_m1sigma/D");
      fweight_branch_map.emplace("genie_CohMA_Genie", std::vector<double *>{&fweight_genie_cohma_p1sigma, &fweight_genie_cohma_m1sigma});

      fevent_tree->Branch("fweight_genie_cohr0_p1sigma", &fweight_genie_cohr0_p1sigma, "fweight_genie_cohr0_p1sigma/D");
      fevent_tree->Branch("fweight_genie_cohr0_m1sigma", &fweight_genie_cohr0_m1sigma, "fweight_genie_cohr0_m1sigma/D");
      fweight_branch_map.emplace("genie_CohR0_Genie", std::vector<double *>{&fweight_genie_cohr0_p1sigma, &fweight_genie_cohr0_m1sigma});

      fevent_tree->Branch("fweight_genie_nonresrvp1pi_p1sigma", &fweight_genie_nonresrvp1pi_p1sigma, "fweight_genie_nonresrvp1pi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nonresrvp1pi_m1sigma", &fweight_genie_nonresrvp1pi_m1sigma, "fweight_genie_nonresrvp1pi_m1sigma/D");
      fweight_branch_map.emplace("genie_NonResRvp1pi_Genie", std::vector<double *>{&fweight_genie_nonresrvp1pi_p1sigma, &fweight_genie_nonresrvp1pi_m1sigma});

      fevent_tree->Branch("fweight_genie_nonresrvbarp1pi_p1sigma", &fweight_genie_nonresrvbarp1pi_p1sigma, "fweight_genie_nonresrvbarp1pi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nonresrvbarp1pi_m1sigma", &fweight_genie_nonresrvbarp1pi_m1sigma, "fweight_genie_nonresrvbarp1pi_m1sigma/D");
      fweight_branch_map.emplace("genie_NonResRvbarp1pi_Genie", std::vector<double *>{&fweight_genie_nonresrvbarp1pi_p1sigma, &fweight_genie_nonresrvbarp1pi_m1sigma});

      fevent_tree->Branch("fweight_genie_nonresrvp2pi_p1sigma", &fweight_genie_nonresrvp2pi_p1sigma, "fweight_genie_nonresrvp2pi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nonresrvp2pi_m1sigma", &fweight_genie_nonresrvp2pi_m1sigma, "fweight_genie_nonresrvp2pi_m1sigma/D");
      fweight_branch_map.emplace("genie_NonResRvp2pi_Genie", std::vector<double *>{&fweight_genie_nonresrvp2pi_p1sigma, &fweight_genie_nonresrvp2pi_m1sigma});

      fevent_tree->Branch("fweight_genie_nonresrvbarp2pi_p1sigma", &fweight_genie_nonresrvbarp2pi_p1sigma, "fweight_genie_nonresrvbarp2pi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nonresrvbarp2pi_m1sigma", &fweight_genie_nonresrvbarp2pi_m1sigma, "fweight_genie_nonresrvbarp2pi_m1sigma/D");
      fweight_branch_map.emplace("genie_NonResRvbarp2pi_Genie", std::vector<double *>{&fweight_genie_nonresrvbarp2pi_p1sigma, &fweight_genie_nonresrvbarp2pi_m1sigma});

      fevent_tree->Branch("fweight_genie_resdecaygamma_p1sigma", &fweight_genie_resdecaygamma_p1sigma, "fweight_genie_resdecaygamma_p1sigma/D");
      fevent_tree->Branch("fweight_genie_resdecaygamma_m1sigma", &fweight_genie_resdecaygamma_m1sigma, "fweight_genie_resdecaygamma_m1sigma/D");
      fweight_branch_map.emplace("genie_ResDecayGamma_Genie", std::vector<double *>{&fweight_genie_resdecaygamma_p1sigma, &fweight_genie_resdecaygamma_m1sigma});

      fevent_tree->Branch("fweight_genie_resdecayeta_p1sigma", &fweight_genie_resdecayeta_p1sigma, "fweight_genie_resdecayeta_p1sigma/D");
      fevent_tree->Branch("fweight_genie_resdecayeta_m1sigma", &fweight_genie_resdecayeta_m1sigma, "fweight_genie_resdecayeta_m1sigma/D");
      fweight_branch_map.emplace("genie_ResDecayEta_Genie", std::vector<double *>{&fweight_genie_resdecayeta_p1sigma, &fweight_genie_resdecayeta_m1sigma});

      fevent_tree->Branch("fweight_genie_resdecaytheta_p1sigma", &fweight_genie_resdecaytheta_p1sigma, "fweight_genie_resdecaytheta_p1sigma/D");
      fevent_tree->Branch("fweight_genie_resdecaytheta_m1sigma", &fweight_genie_resdecaytheta_m1sigma, "fweight_genie_resdecaytheta_m1sigma/D");
      fweight_branch_map.emplace("genie_ResDecayTheta_Genie", std::vector<double *>{&fweight_genie_resdecaytheta_p1sigma, &fweight_genie_resdecaytheta_m1sigma});

      fevent_tree->Branch("fweight_genie_nc_p1sigma", &fweight_genie_nc_p1sigma, "fweight_genie_nc_p1sigma/D");
      fevent_tree->Branch("fweight_genie_nc_m1sigma", &fweight_genie_nc_m1sigma, "fweight_genie_nc_m1sigma/D");
      fweight_branch_map.emplace("genie_NC_Genie", std::vector<double *>{&fweight_genie_nc_p1sigma, &fweight_genie_nc_m1sigma});

      fevent_tree->Branch("fweight_genie_disath_p1sigma", &fweight_genie_disath_p1sigma, "fweight_genie_disath_p1sigma/D");
      fevent_tree->Branch("fweight_genie_disath_m1sigma", &fweight_genie_disath_m1sigma, "fweight_genie_disath_m1sigma/D");
      fweight_branch_map.emplace("genie_DISAth_Genie", std::vector<double *>{&fweight_genie_disath_p1sigma, &fweight_genie_disath_m1sigma});

      fevent_tree->Branch("fweight_genie_disbth_p1sigma", &fweight_genie_disbth_p1sigma, "fweight_genie_disbth_p1sigma/D");
      fevent_tree->Branch("fweight_genie_disbth_m1sigma", &fweight_genie_disbth_m1sigma, "fweight_genie_disbth_m1sigma/D");
      fweight_branch_map.emplace("genie_DISBth_Genie", std::vector<double *>{&fweight_genie_disbth_p1sigma, &fweight_genie_disbth_m1sigma});

      fevent_tree->Branch("fweight_genie_discv1u_p1sigma", &fweight_genie_discv1u_p1sigma, "fweight_genie_discv1u_p1sigma/D");
      fevent_tree->Branch("fweight_genie_discv1u_m1sigma", &fweight_genie_discv1u_m1sigma, "fweight_genie_discv1u_m1sigma/D");
      fweight_branch_map.emplace("genie_DISCv1u_Genie", std::vector<double *>{&fweight_genie_discv1u_p1sigma, &fweight_genie_discv1u_m1sigma});

      fevent_tree->Branch("fweight_genie_discv2u_p1sigma", &fweight_genie_discv2u_p1sigma, "fweight_genie_discv2u_p1sigma/D");
      fevent_tree->Branch("fweight_genie_discv2u_m1sigma", &fweight_genie_discv2u_m1sigma, "fweight_genie_discv2u_m1sigma/D");
      fweight_branch_map.emplace("genie_DISCv2u_Genie", std::vector<double *>{&fweight_genie_discv2u_p1sigma, &fweight_genie_discv2u_m1sigma});

      fevent_tree->Branch("fweight_genie_disnucl_p1sigma", &fweight_genie_disnucl_p1sigma, "fweight_genie_disnucl_p1sigma/D");
      fevent_tree->Branch("fweight_genie_disnucl_m1sigma", &fweight_genie_disnucl_m1sigma, "fweight_genie_disnucl_m1sigma/D");
      fweight_branch_map.emplace("genie_DISnucl_Genie", std::vector<double *>{&fweight_genie_disnucl_p1sigma, &fweight_genie_disnucl_m1sigma});

      fevent_tree->Branch("fweight_genie_agkyxf_p1sigma", &fweight_genie_agkyxf_p1sigma, "fweight_genie_agkyxf_p1sigma/D");
      fevent_tree->Branch("fweight_genie_agkyxf_m1sigma", &fweight_genie_agkyxf_m1sigma, "fweight_genie_agkyxf_m1sigma/D");
      fweight_branch_map.emplace("genie_AGKYxF_Genie", std::vector<double *>{&fweight_genie_agkyxf_p1sigma, &fweight_genie_agkyxf_m1sigma});

      fevent_tree->Branch("fweight_genie_agkypt_p1sigma", &fweight_genie_agkypt_p1sigma, "fweight_genie_agkypt_p1sigma/D");
      fevent_tree->Branch("fweight_genie_agkypt_m1sigma", &fweight_genie_agkypt_m1sigma, "fweight_genie_agkypt_m1sigma/D");
      fweight_branch_map.emplace("genie_AGKYpT_Genie", std::vector<double *>{&fweight_genie_agkypt_p1sigma, &fweight_genie_agkypt_m1sigma});

      fevent_tree->Branch("fweight_genie_formzone_p1sigma", &fweight_genie_formzone_p1sigma, "fweight_genie_formzone_p1sigma/D");
      fevent_tree->Branch("fweight_genie_formzone_m1sigma", &fweight_genie_formzone_m1sigma, "fweight_genie_formzone_m1sigma/D");
      fweight_branch_map.emplace("genie_FormZone_Genie", std::vector<double *>{&fweight_genie_formzone_p1sigma, &fweight_genie_formzone_m1sigma});

      fevent_tree->Branch("fweight_genie_fermigasmodelkf_p1sigma", &fweight_genie_fermigasmodelkf_p1sigma, "fweight_genie_fermigasmodelkf_p1sigma/D");
      fevent_tree->Branch("fweight_genie_fermigasmodelkf_m1sigma", &fweight_genie_fermigasmodelkf_m1sigma, "fweight_genie_fermigasmodelkf_m1sigma/D");
      fweight_branch_map.emplace("genie_FermiGasModelKf_Genie", std::vector<double *>{&fweight_genie_fermigasmodelkf_p1sigma, &fweight_genie_fermigasmodelkf_m1sigma});

      fevent_tree->Branch("fweight_genie_fermigasmodelsf_p1sigma", &fweight_genie_fermigasmodelsf_p1sigma, "fweight_genie_fermigasmodelsf_p1sigma/D");
      fevent_tree->Branch("fweight_genie_fermigasmodelsf_m1sigma", &fweight_genie_fermigasmodelsf_m1sigma, "fweight_genie_fermigasmodelsf_m1sigma/D");
      fweight_branch_map.emplace("genie_FermiGasModelSf_Genie", std::vector<double *>{&fweight_genie_fermigasmodelsf_p1sigma, &fweight_genie_fermigasmodelsf_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukenmfp_p1sigma", &fweight_genie_intranukenmfp_p1sigma, "fweight_genie_intranukenmfp_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukenmfp_m1sigma", &fweight_genie_intranukenmfp_m1sigma, "fweight_genie_intranukenmfp_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNmfp_Genie", std::vector<double *>{&fweight_genie_intranukenmfp_p1sigma, &fweight_genie_intranukenmfp_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukencex_p1sigma", &fweight_genie_intranukencex_p1sigma, "fweight_genie_intranukencex_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukencex_m1sigma", &fweight_genie_intranukencex_m1sigma, "fweight_genie_intranukencex_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNcex_Genie", std::vector<double *>{&fweight_genie_intranukencex_p1sigma, &fweight_genie_intranukencex_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukenel_p1sigma", &fweight_genie_intranukenel_p1sigma, "fweight_genie_intranukenel_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukenel_m1sigma", &fweight_genie_intranukenel_m1sigma, "fweight_genie_intranukenel_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNel_Genie", std::vector<double *>{&fweight_genie_intranukenel_p1sigma, &fweight_genie_intranukenel_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukeninel_p1sigma", &fweight_genie_intranukeninel_p1sigma, "fweight_genie_intranukeninel_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukeninel_m1sigma", &fweight_genie_intranukeninel_m1sigma, "fweight_genie_intranukeninel_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNinel_Genie", std::vector<double *>{&fweight_genie_intranukeninel_p1sigma, &fweight_genie_intranukeninel_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukenabs_p1sigma", &fweight_genie_intranukenabs_p1sigma, "fweight_genie_intranukenabs_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukenabs_m1sigma", &fweight_genie_intranukenabs_m1sigma, "fweight_genie_intranukenabs_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNabs_Genie", std::vector<double *>{&fweight_genie_intranukenabs_p1sigma, &fweight_genie_intranukenabs_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukenpi_p1sigma", &fweight_genie_intranukenpi_p1sigma, "fweight_genie_intranukenpi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukenpi_m1sigma", &fweight_genie_intranukenpi_m1sigma, "fweight_genie_intranukenpi_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukeNpi_Genie", std::vector<double *>{&fweight_genie_intranukenpi_p1sigma, &fweight_genie_intranukenpi_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepimfp_p1sigma", &fweight_genie_intranukepimfp_p1sigma, "fweight_genie_intranukepimfp_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepimfp_m1sigma", &fweight_genie_intranukepimfp_m1sigma, "fweight_genie_intranukepimfp_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePImfp_Genie", std::vector<double *>{&fweight_genie_intranukepimfp_p1sigma, &fweight_genie_intranukepimfp_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepicex_p1sigma", &fweight_genie_intranukepicex_p1sigma, "fweight_genie_intranukepicex_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepicex_m1sigma", &fweight_genie_intranukepicex_m1sigma, "fweight_genie_intranukepicex_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePIcex_Genie", std::vector<double *>{&fweight_genie_intranukepicex_p1sigma, &fweight_genie_intranukepicex_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepiel_p1sigma", &fweight_genie_intranukepiel_p1sigma, "fweight_genie_intranukepiel_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepiel_m1sigma", &fweight_genie_intranukepiel_m1sigma, "fweight_genie_intranukepiel_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePIel_Genie", std::vector<double *>{&fweight_genie_intranukepiel_p1sigma, &fweight_genie_intranukepiel_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepiinel_p1sigma", &fweight_genie_intranukepiinel_p1sigma, "fweight_genie_intranukepiinel_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepiinel_m1sigma", &fweight_genie_intranukepiinel_m1sigma, "fweight_genie_intranukepiinel_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePIinel_Genie", std::vector<double *>{&fweight_genie_intranukepiinel_p1sigma, &fweight_genie_intranukepiinel_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepiabs_p1sigma", &fweight_genie_intranukepiabs_p1sigma, "fweight_genie_intranukepiabs_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepiabs_m1sigma", &fweight_genie_intranukepiabs_m1sigma, "fweight_genie_intranukepiabs_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePIabs_Genie", std::vector<double *>{&fweight_genie_intranukepiabs_p1sigma, &fweight_genie_intranukepiabs_m1sigma});

      fevent_tree->Branch("fweight_genie_intranukepipi_p1sigma", &fweight_genie_intranukepipi_p1sigma, "fweight_genie_intranukepipi_p1sigma/D");
      fevent_tree->Branch("fweight_genie_intranukepipi_m1sigma", &fweight_genie_intranukepipi_m1sigma, "fweight_genie_intranukepipi_m1sigma/D");
      fweight_branch_map.emplace("genie_IntraNukePIpi_Genie", std::vector<double *>{&fweight_genie_intranukepipi_p1sigma, &fweight_genie_intranukepipi_m1sigma});
    }

  }    

}


void FillLightEvent::beginSubRun(art::SubRun const & sr) {

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


void FillLightEvent::ResetEvent() {

  fopflash_producer_indices.clear();
  fhit_producer_indices.clear();
  ftrack_producer_indices.clear();
  fshower_producer_indices.clear();

  fopflash_producer_indices.reserve(fopflashp_size);
  fhit_producer_indices.reserve(fhitp_size);
  ftrack_producer_indices.reserve(ftrackp_size);
  fshower_producer_indices.reserve(fshowerp_size);

  frun_number = -1;
  fsubrun_number = -1;
  fevent_number = -1;

  //fswtrigger_producer_index.clear();
  fswtrigger_getListOfAlgorithms.clear();
  fswtrigger_passedAlgo.clear();
  fswtrigger_passedPrescaleAlgo.clear();
  fswtrigger_vetoAlgo.clear();
  fswtrigger_passedAlgos = -1;
  fswtrigger_vetoAlgos = -1;
  fswtrigger_passedPrescaleAlgos = -1;
  fswtrigger_passed_swtrigger = -1;

  freco_opflash_producer_index.clear();
  freco_opflash_Time.clear();
  freco_opflash_TimeWidth.clear();
  freco_opflash_AbsTime.clear();
  freco_opflash_Frame.clear();
  //freco_opflash_PEs.clear();
  freco_opflash_YCenter.clear();
  freco_opflash_YWidth.clear();
  freco_opflash_ZCenter.clear();
  freco_opflash_ZWidth.clear();
  freco_opflash_InBeamFrame.clear();
  freco_opflash_OnBeamFrame.clear();
  freco_opflash_WireCenters.clear();
  freco_opflash_WireWidths.clear();
  freco_opflash_TotalPE.clear();
  freco_opflash_FastToTotal.clear();

  freco_hit_producer_index.clear();
  freco_hit_StartTick.clear();
  freco_hit_EndTick.clear();
  freco_hit_PeakTime.clear();
  freco_hit_SigmaPeakTime.clear();
  freco_hit_RMS.clear();
  freco_hit_PeakAmplitude.clear();
  freco_hit_SigmaPeakAmplitude.clear();
  freco_hit_SummedADC.clear();
  freco_hit_Integral.clear();
  freco_hit_SigmaIntegral.clear();
  freco_hit_Multiplicity.clear();
  freco_hit_LocalIndex.clear();
  freco_hit_GoodnessOfFit.clear();
  freco_hit_DegreesOfFreedom.clear();
  freco_hit_View.clear();
  freco_hit_SignalType.clear();
  freco_hit_WireID_CryostatID.clear();
  freco_hit_WireID_TPCID.clear();
  freco_hit_WireID_PlaneID.clear();
  freco_hit_WireID_WireID.clear();
  if(frmcm_first) {
    freco_hit_mc_type.clear();
    freco_hit_mc_index.clear();
    /*
    freco_hit_true_ideFraction.clear();
    freco_hit_true_isMaxIDE.clear();
    freco_hit_true_ideNFraction.clear();
    freco_hit_true_isMaxIDEN.clear();
    */
    freco_hit_true_numElectrons.clear();
    freco_hit_true_energy.clear();
  }

  freco_track_producer_index.clear();
  freco_track_NumberTrajectoryPoints.clear();
  freco_track_NPoints.clear();
  freco_track_FirstPoint.clear();
  freco_track_LastPoint.clear();
  freco_track_FirstValidPoint.clear();
  freco_track_LastValidPoint.clear();
  freco_track_CountValidPoints.clear();
  freco_track_X.clear();
  freco_track_Y.clear();
  freco_track_Z.clear();
  freco_track_Px.clear();
  freco_track_Py.clear();
  freco_track_Pz.clear();
  freco_track_HasMomentum.clear();
  freco_track_Length.clear();
  freco_track_Chi2.clear();
  freco_track_Chi2PerNdof.clear();
  freco_track_Ndof.clear();
  freco_track_ParticleId.clear();
  freco_track_Theta.clear();
  freco_track_Phi.clear();
  freco_track_ZenithAngle.clear();
  freco_track_AzimuthAngle.clear();
  freco_track_VertexDirection_X.clear();
  freco_track_VertexDirection_Y.clear();
  freco_track_VertexDirection_Z.clear();
  freco_track_to_reco_hit.clear();
  freco_track_EnergyHelper_resrange.clear();
  freco_track_EnergyHelper_dedx.clear();
  freco_track_EnergyHelper_energy.clear();
  freco_track_EnergyHelperNew_energy_legacy.clear();
  freco_track_EnergyHelperNew_energy.clear();
  freco_track_EnergyHelperNew_energy_from_dedx.clear();
  freco_track_EnergyHelperNew_dedx.clear();
  if(frmcm_first) {
    freco_track_largest_mc_type.clear();
    freco_track_largest_mc_index.clear();
    freco_track_largest_ratio.clear();
    freco_track_mc_type.clear();
    freco_track_mc_index.clear();
    freco_track_charge_contribution.clear();
    freco_track_charge_total.clear();
  }    

  freco_shower_producer_index.clear();
  freco_shower_Direction_x.clear();
  freco_shower_Direction_y.clear();
  freco_shower_Direction_z.clear();
  freco_shower_DirectionErr_x.clear();
  freco_shower_DirectionErr_y.clear();
  freco_shower_DirectionErr_z.clear();
  freco_shower_ShowerStart_x.clear();
  freco_shower_ShowerStart_y.clear();
  freco_shower_ShowerStart_z.clear();
  freco_shower_ShowerStartErr_x.clear();
  freco_shower_ShowerStartErr_y.clear();
  freco_shower_ShowerStartErr_z.clear();
  freco_shower_Energy.clear();
  freco_shower_EnergyErr.clear();
  freco_shower_MIPEnergy.clear();
  freco_shower_MIPEnergyErr.clear();
  freco_shower_best_plane.clear();
  freco_shower_Length.clear();
  freco_shower_OpenAngle.clear();
  freco_shower_dEdx.clear();
  freco_shower_dEdxErr.clear();
  freco_shower_has_open_angle.clear();
  freco_shower_has_length.clear();
  freco_shower_to_reco_hit.clear();
  freco_shower_EnergyHelper_energy_legacy.clear();
  freco_shower_EnergyHelper_energy.clear();
  freco_shower_EnergyHelper_dedx.clear();
  freco_shower_EnergyHelperNew_energy_legacy.clear();
  freco_shower_EnergyHelperNew_energy.clear();
  freco_shower_EnergyHelperNew_dedx.clear();
  if(frmcm_first) {
    freco_shower_largest_mc_type.clear();
    freco_shower_largest_mc_index.clear();
    freco_shower_largest_ratio.clear();
    freco_shower_mc_type.clear();
    freco_shower_mc_index.clear();
    freco_shower_charge_contribution.clear();
    freco_shower_charge_total.clear();
  }

  fpfp_pdg.clear();
  fpfp_vertex_X.clear();
  fpfp_vertex_Y.clear();
  fpfp_vertex_Z.clear();
  fpfp_original_index.clear();
  fpfp_children.clear();

  if(fmc) {

    fnu_pdg.clear();
    fnu_energy.clear();
    flep_pdg.clear();
    flep_energy.clear();
    fccnc.clear();
    fmode.clear();
    finteraction_type.clear();

    ftrue_nuvertx.clear();
    ftrue_nuverty.clear();
    ftrue_nuvertz.clear();

    ftrue_nu_E.clear();

    ftrue_nu_vtx_tpc_contained.clear();
    ftrue_nu_vtx_fid_contained.clear();

    fgenie_particle_TrackId.clear();
    fgenie_particle_StatusCode.clear();
    fgenie_particle_PdgCode.clear();
    fgenie_particle_Mother.clear();
    fgenie_particle_X.clear();
    fgenie_particle_Y.clear();
    fgenie_particle_Z.clear();
    fgenie_particle_T.clear();
    fgenie_particle_Px.clear();
    fgenie_particle_Py.clear();
    fgenie_particle_Pz.clear();
    fgenie_particle_E.clear();

    fmcparticle_TrackId.clear();
    fmcparticle_StatusCode.clear();
    fmcparticle_PdgCode.clear();
    fmcparticle_Mother.clear();      

    /*  
    fmcparticle_X.clear();
    fmcparticle_Y.clear();
    fmcparticle_Z.clear();
    fmcparticle_T.clear();
    fmcparticle_Px.clear();
    fmcparticle_Py.clear();
    fmcparticle_Pz.clear();
    fmcparticle_E.clear();
    */

    fmctrack_Origin.clear();
    fmctrack_PdgCode.clear();
    fmctrack_TrackID.clear();
    fmctrack_Process.clear();
    fmctrack_Start_X.clear();
    fmctrack_Start_Y.clear();
    fmctrack_Start_Z.clear();
    fmctrack_Start_T.clear();
    fmctrack_Start_Px.clear();
    fmctrack_Start_Py.clear();
    fmctrack_Start_Pz.clear();
    fmctrack_Start_E.clear();
    fmctrack_End_X.clear();
    fmctrack_End_Y.clear();
    fmctrack_End_Z.clear();
    fmctrack_End_T.clear();
    fmctrack_End_Px.clear();
    fmctrack_End_Py.clear();
    fmctrack_End_Pz.clear();
    fmctrack_End_E.clear();
    fmctrack_X.clear();
    fmctrack_Y.clear();
    fmctrack_Z.clear();
    fmctrack_T.clear();
    fmctrack_Px.clear();
    fmctrack_Py.clear();
    fmctrack_Pz.clear();
    fmctrack_E.clear();
    fmctrack_dQdx.clear();
    fmctrack_dEdx.clear();
    fmctrack_MotherPdgCode.clear();
    fmctrack_MotherTrackID.clear();
    fmctrack_MotherProcess.clear();
    fmctrack_AncestorPdgCode.clear();
    fmctrack_AncestorTrackID.clear();
    fmctrack_AncestorProcess.clear();
    fmctrack_contributed_charge.clear();

    fmcshower_Origin.clear();
    fmcshower_PdgCode.clear();
    fmcshower_TrackID.clear();
    fmcshower_Process.clear();
    fmcshower_Start_X.clear();
    fmcshower_Start_Y.clear();
    fmcshower_Start_Z.clear();
    fmcshower_Start_T.clear();
    fmcshower_Start_Px.clear();
    fmcshower_Start_Py.clear();
    fmcshower_Start_Pz.clear();
    fmcshower_Start_E.clear();
    fmcshower_End_X.clear();
    fmcshower_End_Y.clear();
    fmcshower_End_Z.clear();
    fmcshower_End_T.clear();
    fmcshower_End_Px.clear();
    fmcshower_End_Py.clear();
    fmcshower_End_Pz.clear();
    fmcshower_End_E.clear();
    fmcshower_MotherPdgCode.clear();
    fmcshower_MotherTrackID.clear();
    fmcshower_MotherProcess.clear();
    fmcshower_AncestorPdgCode.clear();
    fmcshower_AncestorTrackID.clear();
    fmcshower_AncestorProcess.clear();
    fmcshower_DetProfile_X.clear();
    fmcshower_DetProfile_Y.clear();
    fmcshower_DetProfile_Z.clear();
    fmcshower_DetProfile_T.clear();
    fmcshower_DetProfile_Px.clear();
    fmcshower_DetProfile_Py.clear();
    fmcshower_DetProfile_Pz.clear();
    fmcshower_DetProfile_E.clear();
    fmcshower_DaughterTrackID.clear();
    fmcshower_Charge.clear();
    fmcshower_dQdx.clear();
    fmcshower_StartDir_X.clear();
    fmcshower_StartDir_Y.clear();
    fmcshower_StartDir_Z.clear();
    fmcshower_contributed_charge.clear();

    fdelta_mct_index = -1;
    fis_delta_rad.clear();
    fdelta_index.clear();
    fdelta_photon_index.clear();
    fdelta_mcshower_index.clear();
    fdelta_proton_index.clear();
    fdelta_mctrack_index.clear();

    fweight_genie_ncelaxial_p1sigma = -1;
    fweight_genie_ncelaxial_m1sigma = -1;
    fweight_genie_nceleta_p1sigma = -1;
    fweight_genie_nceleta_m1sigma = -1;
    fweight_genie_qema_p1sigma = -1;
    fweight_genie_qema_m1sigma = -1;
    fweight_genie_qevec_p1sigma = -1;
    fweight_genie_qevec_m1sigma = -1;
    fweight_genie_ccresaxial_p1sigma = -1;
    fweight_genie_ccresaxial_m1sigma = -1;
    fweight_genie_ccresvector_p1sigma = -1;
    fweight_genie_ccresvector_m1sigma = -1;
    fweight_genie_resganged_p1sigma = -1;
    fweight_genie_resganged_m1sigma = -1;
    fweight_genie_ncresaxial_p1sigma = -1;
    fweight_genie_ncresaxial_m1sigma = -1;
    fweight_genie_ncresvector_p1sigma = -1;
    fweight_genie_ncresvector_m1sigma = -1;
    fweight_genie_cohma_p1sigma = -1;
    fweight_genie_cohma_m1sigma = -1;

    fweight_genie_cohr0_p1sigma = -1;
    fweight_genie_cohr0_m1sigma = -1;
    fweight_genie_nonresrvp1pi_p1sigma = -1;
    fweight_genie_nonresrvp1pi_m1sigma = -1;
    fweight_genie_nonresrvbarp1pi_p1sigma = -1;
    fweight_genie_nonresrvbarp1pi_m1sigma = -1;
    fweight_genie_nonresrvp2pi_p1sigma = -1;
    fweight_genie_nonresrvp2pi_m1sigma = -1;
    fweight_genie_nonresrvbarp2pi_p1sigma = -1;
    fweight_genie_nonresrvbarp2pi_m1sigma = -1;
    fweight_genie_resdecaygamma_p1sigma = -1;
    fweight_genie_resdecaygamma_m1sigma = -1;
    fweight_genie_resdecayeta_p1sigma = -1;
    fweight_genie_resdecayeta_m1sigma = -1;
    fweight_genie_resdecaytheta_p1sigma = -1;
    fweight_genie_resdecaytheta_m1sigma = -1;
    fweight_genie_nc_p1sigma = -1;
    fweight_genie_nc_m1sigma = -1;
    fweight_genie_disath_p1sigma = -1;
    fweight_genie_disath_m1sigma = -1;
    fweight_genie_disbth_p1sigma = -1;
    fweight_genie_disbth_m1sigma = -1;
    fweight_genie_discv1u_p1sigma = -1;
    fweight_genie_discv1u_m1sigma = -1;
    fweight_genie_discv2u_p1sigma = -1;
    fweight_genie_discv2u_m1sigma = -1;
    fweight_genie_disnucl_p1sigma = -1;
    fweight_genie_disnucl_m1sigma = -1;
    fweight_genie_agkyxf_p1sigma = -1;
    fweight_genie_agkyxf_m1sigma = -1;
    fweight_genie_agkypt_p1sigma = -1;
    fweight_genie_agkypt_m1sigma = -1;
    fweight_genie_formzone_p1sigma = -1;
    fweight_genie_formzone_m1sigma = -1;
    fweight_genie_fermigasmodelkf_p1sigma = -1;
    fweight_genie_fermigasmodelkf_m1sigma = -1;
    fweight_genie_fermigasmodelsf_p1sigma = -1;
    fweight_genie_fermigasmodelsf_m1sigma = -1;
    fweight_genie_intranukenmfp_p1sigma = -1;
    fweight_genie_intranukenmfp_m1sigma = -1;
    fweight_genie_intranukencex_p1sigma = -1;
    fweight_genie_intranukencex_m1sigma = -1;
    fweight_genie_intranukenel_p1sigma = -1;
    fweight_genie_intranukenel_m1sigma = -1;
    fweight_genie_intranukeninel_p1sigma = -1;
    fweight_genie_intranukeninel_m1sigma = -1;
    fweight_genie_intranukenabs_p1sigma = -1;
    fweight_genie_intranukenabs_m1sigma = -1;
    fweight_genie_intranukenpi_p1sigma = -1;
    fweight_genie_intranukenpi_m1sigma = -1;
    fweight_genie_intranukepimfp_p1sigma = -1;
    fweight_genie_intranukepimfp_m1sigma = -1;
    fweight_genie_intranukepicex_p1sigma = -1;
    fweight_genie_intranukepicex_m1sigma = -1;
    fweight_genie_intranukepiel_p1sigma = -1;
    fweight_genie_intranukepiel_m1sigma = -1;
    fweight_genie_intranukepiinel_p1sigma = -1;
    fweight_genie_intranukepiinel_m1sigma = -1;
    fweight_genie_intranukepiabs_p1sigma = -1;
    fweight_genie_intranukepiabs_m1sigma = -1;
    fweight_genie_intranukepipi_p1sigma = -1;
    fweight_genie_intranukepipi_m1sigma = -1;

  }

}


bool FillLightEvent::PassedSWTrigger(art::Event const & e, std::string const & swtrigger_product) const {
  
  art::ValidHandle<::raw::ubdaqSoftwareTriggerData> const & swt = e.getValidHandle<::raw::ubdaqSoftwareTriggerData>(swtrigger_product);
  //art::ValidHandle<::raw::ubdaqSoftwareTriggerData> const & swt =    e.getValidHandle<::raw::ubdaqSoftwareTriggerData>("daq");
  //art::ValidHandle<::raw::ubdaqSoftwareTriggerData> const & swt =    e.getValidHandle<::raw::ubdaqSoftwareTriggerData>("swtrigger");

  std::vector<std::string> const & algo_v = swt->getListOfAlgorithms();

  std::cout<<"ALGO_V: of size: "<<algo_v.size()<<std::endl;
 // for(std::string const &a: algo_v){
	//std::cout<<a<<std::endl;	
 // }

  std::string const int_str = "BNB_FEMBeamTriggerAlgo";
  std::string const ext_str = "EXT_BNBwin_FEMBeamTriggerAlgo";

  auto const int_it = std::find(algo_v.begin(), algo_v.end(), int_str);
  auto const ext_it = std::find(algo_v.begin(), algo_v.end(), ext_str);

  if(int_it != algo_v.end() && ext_it != algo_v.end()) {
    std::cout << "function: " << __PRETTY_FUNCTION__ << " line: " << __LINE__ << std::endl
	      << "Found both swtriggers\n";
    return false;
  }
  else if(int_it == algo_v.end() && ext_it == algo_v.end()) {
    std::cout << "function: " << __PRETTY_FUNCTION__ << " line: " << __LINE__ << std::endl
	      << "Found neither swtrigger\n";
    return false;
  }
  else if(int_it != algo_v.end()) {
    return swt->passedAlgo(int_str);
  }      
  else if(ext_it != algo_v.end()) {
    return swt->passedAlgo(ext_str);
  }

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  return false;

}


void FillLightEvent::FillSWTriggerVectors(art::Event const & e) {

  art::ValidHandle<::raw::ubdaqSoftwareTriggerData> const & swt = e.getValidHandle<::raw::ubdaqSoftwareTriggerData>(fswtrigger_product);
  std::vector<std::string> const & algo_v = swt->getListOfAlgorithms();
  size_t const algo_v_size = algo_v.size();

  //fswtrigger_producer_index.push_back(producer_index);

  fswtrigger_getListOfAlgorithms = algo_v;

  fswtrigger_passedAlgo.reserve(algo_v_size);
  fswtrigger_passedPrescaleAlgo.reserve(algo_v_size);
  fswtrigger_vetoAlgo.reserve(algo_v_size);
  
  for(std::string const & algo : algo_v) {

    fswtrigger_passedAlgo.push_back(swt->passedAlgo(algo));
    fswtrigger_passedPrescaleAlgo.push_back(swt->passedPrescaleAlgo(algo));
    fswtrigger_vetoAlgo.push_back(swt->vetoAlgo(algo));

  }

  fswtrigger_passedAlgos = swt->passedAlgos(algo_v);
  fswtrigger_vetoAlgos = swt->vetoAlgos(algo_v);
  fswtrigger_passedPrescaleAlgos = swt->passedPrescaleAlgos(algo_v);

  fswtrigger_passed_swtrigger = PassedSWTrigger(e, fswtrigger_product);

}


void FillLightEvent::FillRecoOpFlashVectors(art::Event const & e, int const producer_index, size_t const opflash_index_offset) {

  art::ValidHandle<std::vector<recob::OpFlash>> const & ev_opf = e.getValidHandle<std::vector<recob::OpFlash>>(fopflash_producers.at(producer_index));

  for(size_t i = 0; i < ev_opf->size(); ++i) {

    recob::OpFlash const & opf = ev_opf->at(i);

    //size_t const index = i + opflash_index_offset;

    freco_opflash_producer_index.push_back(producer_index);
    freco_opflash_Time.push_back(opf.Time());
    freco_opflash_TimeWidth.push_back(opf.TimeWidth());
    freco_opflash_AbsTime.push_back(opf.AbsTime());
    freco_opflash_Frame.push_back(opf.Frame());
    //freco_opflash_PE.push_back(opf.PEs());
    freco_opflash_YCenter.push_back(opf.YCenter());
    freco_opflash_YWidth.push_back(opf.YWidth());
    freco_opflash_ZCenter.push_back(opf.ZCenter());
    freco_opflash_ZWidth.push_back(opf.ZWidth());
    freco_opflash_InBeamFrame.push_back(opf.InBeamFrame());
    freco_opflash_OnBeamFrame.push_back(opf.OnBeamTime());
    freco_opflash_WireCenters.push_back(opf.WireCenters());
    freco_opflash_WireWidths.push_back(opf.WireWidths());
    freco_opflash_TotalPE.push_back(opf.TotalPE());
    freco_opflash_FastToTotal.push_back(opf.FastToTotal());      

  }

}


void FillLightEvent::FillRecoOpFlashVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < fopflashp_size; ++i) {
    art::ValidHandle<std::vector<recob::OpFlash>> const & ev_opf = e.getValidHandle<std::vector<recob::OpFlash>>(fopflash_producers.at(i));
    size += ev_opf->size();
    fopflash_producer_indices.push_back(size);
  }

  freco_opflash_producer_index.reserve(size);
  freco_opflash_Time.reserve(size);
  freco_opflash_TimeWidth.reserve(size);
  freco_opflash_AbsTime.reserve(size);
  freco_opflash_Frame.reserve(size);
  //freco_opflash_PEs.reserve(size);
  freco_opflash_YCenter.reserve(size);
  freco_opflash_YWidth.reserve(size);
  freco_opflash_ZCenter.reserve(size);
  freco_opflash_ZWidth.reserve(size);
  freco_opflash_InBeamFrame.reserve(size);
  freco_opflash_OnBeamFrame.reserve(size);
  freco_opflash_WireCenters.reserve(size);
  freco_opflash_WireWidths.reserve(size);
  freco_opflash_TotalPE.reserve(size);
  freco_opflash_FastToTotal.reserve(size);  

  size_t opflash_index_offset = 0;
  for(size_t i = 0; i < fopflashp_size; ++i) {
    FillRecoOpFlashVectors(e, i, opflash_index_offset);
    art::ValidHandle<std::vector<recob::OpFlash>> const & ev_opf = e.getValidHandle<std::vector<recob::OpFlash>>(fopflash_producers.at(i));
    opflash_index_offset += ev_opf->size();
  }    

}


void FillLightEvent::FillRecoHitVectors(art::Event const & e, int const producer_index, size_t const hit_index_offset) {

  art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(producer_index));  

  std::unordered_map<int, size_t> const * tp_map = nullptr;
  std::unordered_map<int, size_t> const * sp_map = nullptr;
  std::unordered_map<int, size_t> const * mcp_map = nullptr;
  art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> * particles_per_hit = nullptr;
  std::cout<<"JustBefore "<<std::endl;
  if(frmcm_first) {
  std::cout<<"JustAFTER "<<std::endl;
    tp_map = &frmcm_first->GetMCTrackMap();
    sp_map = &frmcm_first->GetMCShowerMap();
    mcp_map = &frmcm_first->GetMCParticleMap();
    particles_per_hit = new art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>(ev_h, e, frmcmassociation_producers.at(producer_index));
  }
  std::vector<simb::MCParticle const *> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

  for(size_t i = 0; i < ev_h->size(); ++i) {

    recob::Hit const & h = ev_h->at(i);

    size_t const index = i + hit_index_offset;

    freco_hit_producer_index.push_back(producer_index);
    freco_hit_StartTick.push_back(h.StartTick());
    freco_hit_EndTick.push_back(h.EndTick());
    freco_hit_PeakTime.push_back(h.PeakTime());
    freco_hit_SigmaPeakTime.push_back(h.SigmaPeakTime());
    freco_hit_RMS.push_back(h.RMS());
    freco_hit_PeakAmplitude.push_back(h.PeakAmplitude());
    freco_hit_SigmaPeakAmplitude.push_back(h.SigmaPeakAmplitude());
    freco_hit_SummedADC.push_back(h.SummedADC());
    freco_hit_Integral.push_back(h.Integral());
    freco_hit_SigmaIntegral.push_back(h.SigmaIntegral());
    freco_hit_Multiplicity.push_back(h.Multiplicity());
    freco_hit_LocalIndex.push_back(h.LocalIndex());
    freco_hit_GoodnessOfFit.push_back(h.GoodnessOfFit());
    freco_hit_DegreesOfFreedom.push_back(h.DegreesOfFreedom());
    freco_hit_View.push_back(h.View());
    freco_hit_SignalType.push_back(h.SignalType());
    freco_hit_WireID_CryostatID.push_back(h.WireID().Cryostat);
    freco_hit_WireID_TPCID.push_back(h.WireID().TPC);
    freco_hit_WireID_PlaneID.push_back(h.WireID().Plane);
    freco_hit_WireID_WireID.push_back(h.WireID().Wire);

    if(frmcm_first) {

      particle_vec.clear(); 
      match_vec.clear();
      particles_per_hit->get(i, particle_vec, match_vec);
      size_t const particle_vec_size = particle_vec.size();

      std::vector<int> & reco_hit_mc_type = freco_hit_mc_type.at(index);
      std::vector<int> & reco_hit_mc_index = freco_hit_mc_index.at(index);
      /*
      std::vector<float> & reco_hit_true_ideFraction = freco_hit_true_ideFraction.at(index);
      std::vector<int> & reco_hit_true_isMaxIDE = freco_hit_true_isMaxIDE.at(index);
      std::vector<float> & reco_hit_true_ideNFraction = freco_hit_true_ideNFraction.at(index);
      std::vector<int> & reco_hit_true_isMaxIDEN = freco_hit_true_isMaxIDEN.at(index);
      */ 
      std::vector<float> & reco_hit_true_numElectrons = freco_hit_true_numElectrons.at(index);
      std::vector<float> & reco_hit_true_energy = freco_hit_true_energy.at(index);
      
      freco_hit_mc_type.reserve(particle_vec_size);
      freco_hit_mc_index.reserve(particle_vec_size);
      /*  
      freco_hit_true_ideFraction.reserve(particle_vec_size);
      freco_hit_true_isMaxIDE.reserve(particle_vec_size);
      freco_hit_true_ideNFraction.reserve(particle_vec_size);
      freco_hit_true_isMaxIDEN.reserve(particle_vec_size);
      */
      freco_hit_true_numElectrons.reserve(particle_vec_size);
      freco_hit_true_energy.reserve(particle_vec_size);      

      for(size_t i_p = 0; i_p < particle_vec_size; ++i_p) {
	anab::BackTrackerHitMatchingData const * match = match_vec[i_p];
	/*
	reco_hit_true_ideFraction.push_back(match->ideFraction);
	reco_hit_true_isMaxIDE.push_back(match->isMaxIDE);
	reco_hit_true_ideNFraction.push_back(match->ideNFraction);
	reco_hit_true_isMaxIDEN.push_back(match->isMaxIDEN);
	*/
	reco_hit_true_numElectrons.push_back(match->numElectrons);
	reco_hit_true_energy.push_back(match->energy);      	
	int const trackid = particle_vec[i_p]->TrackId();
	auto const tp_it = tp_map->find(trackid);
	if(tp_it != tp_map->end()) {
	  reco_hit_mc_type.push_back(fmc_type_track);
	  reco_hit_mc_index.push_back(tp_it->second);	  
	  continue;
	}
	auto const sp_it = sp_map->find(trackid);
	if(sp_it != sp_map->end()) {
	  reco_hit_mc_type.push_back(fmc_type_shower);
	  reco_hit_mc_index.push_back(sp_it->second);
	  continue;
	}
	auto const mcp_it = mcp_map->find(trackid);
	if(mcp_it != mcp_map->end()) {
	  reco_hit_mc_type.push_back(fmc_type_particle);
	  reco_hit_mc_index.push_back(mcp_it->second);
	  continue;
	}

      }

    }

  }

  if(frmcm_first) delete particles_per_hit;

}


void FillLightEvent::FillRecoHitVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < fhitp_size; ++i) {
    art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(i));
    size += ev_h->size();
    fhit_producer_indices.push_back(size);
  }

  freco_hit_producer_index.reserve(size);
  freco_hit_StartTick.reserve(size);
  freco_hit_EndTick.reserve(size);
  freco_hit_PeakTime.reserve(size);
  freco_hit_SigmaPeakTime.reserve(size);
  freco_hit_RMS.reserve(size);
  freco_hit_PeakAmplitude.reserve(size);
  freco_hit_SigmaPeakAmplitude.reserve(size);
  freco_hit_SummedADC.reserve(size);
  freco_hit_Integral.reserve(size);
  freco_hit_SigmaIntegral.reserve(size);
  freco_hit_Multiplicity.reserve(size);
  freco_hit_LocalIndex.reserve(size);
  freco_hit_GoodnessOfFit.reserve(size);
  freco_hit_DegreesOfFreedom.reserve(size);
  freco_hit_View.reserve(size);
  freco_hit_SignalType.reserve(size);
  freco_hit_WireID_CryostatID.reserve(size);
  freco_hit_WireID_TPCID.reserve(size);
  freco_hit_WireID_PlaneID.reserve(size);
  freco_hit_WireID_WireID.reserve(size);
  if(frmcm_first) {
    freco_hit_mc_type.resize(size, std::vector<int>());
    freco_hit_mc_index.resize(size, std::vector<int>());
    /*
    freco_hit_true_ideFraction.resize(size, std::vector<float>());
    freco_hit_true_isMaxIDE.resize(size, std::vector<int>());
    freco_hit_true_ideNFraction.resize(size, std::vector<float>());
    freco_hit_true_isMaxIDEN.resize(size, std::vector<int>());
    */
    freco_hit_true_numElectrons.resize(size, std::vector<float>());
    freco_hit_true_energy.resize(size, std::vector<float>());
  }

  size_t hit_index_offset = 0;
  for(size_t i = 0; i < fhitp_size; ++i) {
    FillRecoHitVectors(e, i, hit_index_offset);
    art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(i));
    hit_index_offset += ev_h->size();
  }

}


std::pair<std::vector<double>,std::vector<double>> FillLightEvent::GetTrackCaloInfo(art::Event const & e,
										    std::string const & track_producer,
										    std::string const & track_calo_producer,
										    size_t const track_index,
										    double & energy) {

  art::ValidHandle<std::vector<recob::PFParticle>> const & ev_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(track_producer);    
  art::FindManyP<recob::Track> PFPToTrack(ev_pfp, e, track_producer);  

  art::Ptr<recob::Track> const * track_ptr = nullptr;

  for(size_t i = 0; i < ev_pfp->size(); ++i) {
    std::vector<art::Ptr<recob::Track>> const & tracks = PFPToTrack.at(i);
    for(art::Ptr<recob::Track> const & track : tracks) {
      if(track.key() == track_index) {
	track_ptr = &track;
	break;
      }
      if(track_ptr) break;
    }
  }

  if(!track_ptr) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: no track pointer found\n";
    std::pair<std::vector<double>,std::vector<double>> fail({-999},{-999});
    return fail;
  }

  energy = fenergyHelper.trackEnergy(*track_ptr, e, track_producer, track_calo_producer);

  return fenergyHelper.trackdEdx(*track_ptr, e, track_producer, track_calo_producer);
}


std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> FillLightEvent::GetTrackHelperEnergy(art::Event const & e,
													    std::string const & track_producer,
													    size_t const track_index,
													    std::vector<double> & track_energy_dedx) {
  
  art::ValidHandle<std::vector<recob::PFParticle>> const & ev_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(track_producer);
  art::FindManyP<recob::Track> PFPToTrack(ev_pfp, e, track_producer);  

  art::Ptr<recob::Track> const * track_ptr = nullptr;

  for(size_t i = 0; i < ev_pfp->size(); ++i) {
    std::vector<art::Ptr<recob::Track>> const & tracks = PFPToTrack.at(i);
    for(art::Ptr<recob::Track> const & track : tracks) {
      if(track.key() == track_index) {
	track_ptr = &track;
	break;
      }
      if(track_ptr) break;
    }
  }

  if(!track_ptr) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: no track pointer found\n";
    return std::make_pair(std::make_pair(std::vector<double>(), SIZE_MAX), std::vector<double>());
  }

  art::ValidHandle<std::vector<recob::Track>> const & ev_s = e.getValidHandle<std::vector<recob::Track>>(track_producer);
  art::FindManyP<recob::PFParticle> TrackToPFP(ev_s, e, track_producer);
  std::vector<art::Ptr<recob::PFParticle>> track_pfp = TrackToPFP.at(track_index);
  if(track_pfp.size() != 1) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: track_pfp.size() != 1: " << track_pfp.size() << "\n";
  }
  std::vector<double> dqdx_v(3, -1);
  std::vector<double> dqdx_hits_v(3, -1);
  std::vector<double> reco_track_dedx_vector(3, -1);
  std::pair<std::vector<double>, size_t> energy_pair;

  std::vector<int> nHits;
  std::vector<double> pfenergy;
  fenergyHelperNew.energyFromHits(*track_pfp.front(), nHits, pfenergy, e, track_producer);
  track_energy_dedx = fenergyHelperNew.trackEnergy_dedx(*track_ptr, e, track_producer);
  energy_pair = std::make_pair(pfenergy, std::distance(nHits.begin(), std::max_element(nHits.begin(), nHits.end())));
  fenergyHelperNew.dQdx(track_pfp.front().key(), e, dqdx_v, dqdx_hits_v, 4, 1, track_producer);
  fenergyHelperNew.dEdxFromdQdx(reco_track_dedx_vector, dqdx_v);
  
  return std::make_pair(energy_pair, reco_track_dedx_vector);

}


void FillLightEvent::FillRecoTrackVectors(art::Event const & e, int const producer_index, size_t const track_index_offset, size_t const hit_index_offset) {

  std::string const & track_producer = ftrack_producers.at(producer_index);
  art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(track_producer);  
  art::FindManyP<recob::Hit> hits_per_track(ev_t, e, track_producer);
  std::vector<RecoMCMatch> const * track_matches = nullptr;
  if(frmcm_first) {
    RecoMCMatching const & rmcm = frmcm.at(producer_index);
    track_matches = &rmcm.GetTrackMatches();
  }

  for(size_t i = 0; i < ev_t->size(); ++i) {

    recob::Track const & t = ev_t->at(i);
    size_t const traj_size = t.NumberTrajectoryPoints();

    size_t const index = i + track_index_offset;

    freco_track_producer_index.push_back(producer_index);
    freco_track_NumberTrajectoryPoints.push_back(t.NumberTrajectoryPoints());
    freco_track_NPoints.push_back(t.NPoints());
    freco_track_FirstPoint.push_back(t.FirstPoint());
    freco_track_LastPoint.push_back(t.LastPoint());
    freco_track_FirstValidPoint.push_back(t.FirstValidPoint());
    freco_track_LastValidPoint.push_back(t.LastValidPoint());
    freco_track_CountValidPoints.push_back(t.CountValidPoints());
    freco_track_HasMomentum.push_back(t.HasMomentum());
    freco_track_Chi2.push_back(t.Chi2());
    freco_track_Chi2PerNdof.push_back(t.Chi2PerNdof());
    freco_track_Ndof.push_back(t.Ndof());
    freco_track_ParticleId.push_back(t.ParticleId());

    std::vector<double> & reco_track_X_traj = freco_track_X.at(index);
    std::vector<double> & reco_track_Y_traj = freco_track_Y.at(index);
    std::vector<double> & reco_track_Z_traj = freco_track_Z.at(index);
    std::vector<double> & reco_track_Px_traj = freco_track_Px.at(index);
    std::vector<double> & reco_track_Py_traj = freco_track_Py.at(index);
    std::vector<double> & reco_track_Pz_traj = freco_track_Pz.at(index);
    std::vector<double> & reco_track_Length_traj = freco_track_Length.at(index);

    reco_track_X_traj.reserve(traj_size);
    reco_track_Y_traj.reserve(traj_size);
    reco_track_Z_traj.reserve(traj_size);
    reco_track_Px_traj.reserve(traj_size);
    reco_track_Py_traj.reserve(traj_size);
    reco_track_Pz_traj.reserve(traj_size);
    reco_track_Length_traj.reserve(traj_size);
    
    for(size_t j = 0; j < traj_size; ++j) {

      recob::Track::TrajectoryPoint_t const & trajp = t.TrajectoryPoint(j);
      recob::Track::Point_t const & pos = trajp.position;
      recob::Track::Vector_t const & mom = trajp.momentum;

      reco_track_X_traj.push_back(pos.X());
      reco_track_Y_traj.push_back(pos.Y());
      reco_track_Z_traj.push_back(pos.Z());
      reco_track_Px_traj.push_back(mom.X());
      reco_track_Py_traj.push_back(mom.Y());
      reco_track_Pz_traj.push_back(mom.Z());
      reco_track_Length_traj.push_back(t.Length(j));

    }

    freco_track_Theta.push_back(t.Theta());
    freco_track_Phi.push_back(t.Phi());
    freco_track_ZenithAngle.push_back(t.ZenithAngle());
    freco_track_AzimuthAngle.push_back(t.AzimuthAngle());
    freco_track_VertexDirection_X.push_back(t.VertexDirection().X());
    freco_track_VertexDirection_Y.push_back(t.VertexDirection().Y());
    freco_track_VertexDirection_Z.push_back(t.VertexDirection().Z());

    std::vector<int> & reco_track_to_reco_hit = freco_track_to_reco_hit.at(index);
    reco_track_to_reco_hit.reserve(hits_per_track.size());
    for(art::Ptr<recob::Hit> const & hit_ptr : hits_per_track.at(i)) {
      reco_track_to_reco_hit.push_back(hit_ptr.key() + hit_index_offset);
    }

    double track_energy = -1;
    std::pair<std::vector<double>, std::vector<double>> const calo_pair = GetTrackCaloInfo(e, track_producer, ftrack_calo_producer, i, track_energy);
    freco_track_EnergyHelper_resrange.push_back(calo_pair.first);
    freco_track_EnergyHelper_dedx.push_back(calo_pair.second);
    freco_track_EnergyHelper_energy.push_back(track_energy);

    std::vector<double> track_energy_from_dedx;
    std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> const energynew_pair = GetTrackHelperEnergy(e, track_producer, i, track_energy_from_dedx);
    freco_track_EnergyHelperNew_energy_legacy.push_back(energynew_pair.first.first.at(energynew_pair.first.second));
    freco_track_EnergyHelperNew_energy.push_back(energynew_pair.first.first);
    freco_track_EnergyHelperNew_energy_from_dedx.push_back(track_energy_from_dedx);
    freco_track_EnergyHelperNew_dedx.push_back(energynew_pair.second);    

    if(track_matches) {
      
      RecoMCMatch const & track_match = track_matches->at(i);
      
      freco_track_largest_mc_type.push_back(track_match.mc_type);
      freco_track_largest_mc_index.push_back(track_match.mc_index);
      freco_track_largest_ratio.push_back(track_match.ratio);
      
      std::vector<int> & reco_track_mc_type = freco_track_mc_type.at(index);
      std::vector<int> & reco_track_mc_index = freco_track_mc_index.at(index);
      std::vector<double> & reco_track_charge_contribution = freco_track_charge_contribution.at(index);    
      
      size_t const matching_map_size = track_match.mctrack_map.size() + track_match.mcshower_map.size() + track_match.mcparticle_map.size();
      reco_track_mc_type.reserve(matching_map_size);
      reco_track_mc_index.reserve(matching_map_size);
      reco_track_charge_contribution.reserve(matching_map_size);
      
      double charge_total = 0;
      
      for(auto const & p : track_match.mctrack_map) {
	reco_track_mc_type.push_back(fmc_type_track);
	reco_track_mc_index.push_back(p.first);
	reco_track_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }
      
      for(auto const & p : track_match.mcshower_map) {
	reco_track_mc_type.push_back(fmc_type_shower);
	reco_track_mc_index.push_back(p.first);
	reco_track_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }
      
      for(auto const & p : track_match.mcparticle_map) {
	reco_track_mc_type.push_back(fmc_type_particle);
	reco_track_mc_index.push_back(p.first);
	reco_track_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }    
      
      freco_track_charge_total.push_back(charge_total);

    }

  }
  
}


void FillLightEvent::FillRecoTrackVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < ftrackp_size; ++i) {
    art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producers.at(i));
    size += ev_t->size();
    ftrack_producer_indices.push_back(size);
  }

  freco_track_producer_index.reserve(size);
  freco_track_NumberTrajectoryPoints.reserve(size);
  freco_track_NPoints.reserve(size);
  freco_track_FirstPoint.reserve(size);
  freco_track_LastPoint.reserve(size);
  freco_track_FirstValidPoint.reserve(size);
  freco_track_LastValidPoint.reserve(size);
  freco_track_CountValidPoints.reserve(size);
  freco_track_X.resize(size, std::vector<double>());
  freco_track_Y.resize(size, std::vector<double>());
  freco_track_Z.resize(size, std::vector<double>());
  freco_track_Px.resize(size, std::vector<double>());
  freco_track_Py.resize(size, std::vector<double>());
  freco_track_Pz.resize(size, std::vector<double>());
  freco_track_HasMomentum.reserve(size);
  freco_track_Length.resize(size, std::vector<double>());
  freco_track_Chi2.reserve(size);
  freco_track_Chi2PerNdof.reserve(size);
  freco_track_Ndof.reserve(size);
  freco_track_ParticleId.reserve(size);
  freco_track_Theta.reserve(size);
  freco_track_Phi.reserve(size);
  freco_track_ZenithAngle.reserve(size);
  freco_track_AzimuthAngle.reserve(size);
  freco_track_VertexDirection_X.reserve(size);
  freco_track_VertexDirection_Y.reserve(size);
  freco_track_VertexDirection_Z.reserve(size);
  freco_track_to_reco_hit.resize(size, std::vector<int>());
  freco_track_EnergyHelper_resrange.reserve(size);
  freco_track_EnergyHelper_dedx.reserve(size);
  freco_track_EnergyHelper_energy.reserve(size);
  freco_track_EnergyHelperNew_energy_legacy.reserve(size);
  freco_track_EnergyHelperNew_energy.reserve(size);
  freco_track_EnergyHelperNew_energy_from_dedx.reserve(size);
  freco_track_EnergyHelperNew_dedx.reserve(size);
  if(frmcm_first) {
    freco_track_largest_mc_type.reserve(size);
    freco_track_largest_mc_index.reserve(size);
    freco_track_largest_ratio.reserve(size);
    freco_track_mc_type.resize(size, std::vector<int>());
    freco_track_mc_index.resize(size, std::vector<int>());
    freco_track_charge_contribution.resize(size, std::vector<double>());
    freco_track_charge_total.reserve(size);
  }

  size_t track_index_offset = 0;
  size_t hit_index_offset = 0;
  for(size_t i = 0; i < ftrackp_size; ++i) {
    FillRecoTrackVectors(e, i, track_index_offset, hit_index_offset);
    art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producers.at(i));
    track_index_offset += ev_t->size();
    art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(i));
    hit_index_offset += ev_h->size();
  }

}


std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> FillLightEvent::GetShowerHelperEnergy(art::Event const & e,
													     std::string const & shower_producer,
													     size_t const shower_index,
													     bool const energy_helper_new) {
  
  art::ValidHandle<std::vector<recob::PFParticle>> const & ev_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(shower_producer);
  art::FindManyP<recob::Shower> PFPToShower(ev_pfp, e, shower_producer);  
  art::Ptr<recob::Shower> const * shower_ptr = nullptr;

  for(size_t i = 0; i < ev_pfp->size(); ++i) {
    std::vector<art::Ptr<recob::Shower>> const & showers = PFPToShower.at(i);
    for(art::Ptr<recob::Shower> const & shower : showers) {
      if(shower.key() == shower_index) {
	shower_ptr = &shower;
	break;
      }
      if(shower_ptr) break;
    }
  }

  if(!shower_ptr) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: no shower pointer found\n";
    return std::make_pair(std::make_pair(std::vector<double>(), SIZE_MAX), std::vector<double>());
  }

  art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(shower_producer);
  art::FindManyP<recob::PFParticle> ShowerToPFP(ev_s, e, shower_producer);
  std::vector<art::Ptr<recob::PFParticle>> shower_pfp = ShowerToPFP.at(shower_index);
  if(shower_pfp.size() != 1) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: shower_pfp.size() != 1: " << shower_pfp.size() << "\n";
  }
  std::vector<double> dqdx_v(3, -1);
  std::vector<double> dqdx_hits_v(3, -1);
  std::vector<double> reco_shower_dedx_vector(3, -1);
  std::pair<std::vector<double>, size_t> energy_pair;
  
  if(energy_helper_new) {
    std::vector<int> nHits;
    std::vector<double> pfenergy;
    fenergyHelperNew.energyFromHits(*shower_pfp.front(), nHits, pfenergy, e, shower_producer);
    energy_pair = std::make_pair(pfenergy, std::distance(nHits.begin(), std::max_element(nHits.begin(), nHits.end())));
    fenergyHelperNew.dQdx(shower_pfp.front().key(), e, dqdx_v, dqdx_hits_v, 4, 1, shower_producer);
    fenergyHelperNew.dEdxFromdQdx(reco_shower_dedx_vector, dqdx_v);
  }
  else {
    energy_pair = fenergyHelper.showerEnergyV(*shower_ptr, e, shower_producer);
    fenergyHelper.dQdx(shower_pfp.front().key(), e, dqdx_v, dqdx_hits_v, 4, 2, shower_producer);
    fenergyHelper.dEdxFromdQdx(reco_shower_dedx_vector, dqdx_v);
  }

  return std::make_pair(energy_pair, reco_shower_dedx_vector);

}


void FillLightEvent::FillRecoShowerVectors(art::Event const & e, int const producer_index, size_t const shower_index_offset, size_t const hit_index_offset) {

  std::string const & shower_producer = fshower_producers.at(producer_index);
  art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(shower_producer);  
  art::FindManyP<recob::Hit> hits_per_shower(ev_s, e, shower_producer);
  std::vector<RecoMCMatch> const * shower_matches = nullptr;
  if(frmcm_first) {
    RecoMCMatching const & rmcm = frmcm.at(producer_index);
    shower_matches = &rmcm.GetShowerMatches();
  }
  
  for(size_t i = 0; i < ev_s->size(); ++i) {

    recob::Shower const & s = ev_s->at(i);

    size_t const index = i + shower_index_offset;

    freco_shower_producer_index.push_back(producer_index);
    freco_shower_Direction_x.push_back(s.Direction().X());
    freco_shower_Direction_y.push_back(s.Direction().Y());
    freco_shower_Direction_z.push_back(s.Direction().Z());
    freco_shower_DirectionErr_x.push_back(s.DirectionErr().X());
    freco_shower_DirectionErr_y.push_back(s.DirectionErr().Y());
    freco_shower_DirectionErr_z.push_back(s.DirectionErr().Z());
    freco_shower_ShowerStart_x.push_back(s.ShowerStart().X());
    freco_shower_ShowerStart_y.push_back(s.ShowerStart().Y());
    freco_shower_ShowerStart_z.push_back(s.ShowerStart().Z());
    freco_shower_ShowerStartErr_x.push_back(s.ShowerStartErr().X());
    freco_shower_ShowerStartErr_y.push_back(s.ShowerStartErr().Y());
    freco_shower_ShowerStartErr_z.push_back(s.ShowerStartErr().Z());
    freco_shower_Energy.push_back(s.Energy());
    freco_shower_EnergyErr.push_back(s.EnergyErr());
    freco_shower_MIPEnergy.push_back(s.MIPEnergy());
    freco_shower_MIPEnergyErr.push_back(s.MIPEnergyErr());
    freco_shower_best_plane.push_back(s.best_plane());
    freco_shower_Length.push_back(s.Length());
    freco_shower_OpenAngle.push_back(s.OpenAngle());
    freco_shower_dEdx.push_back(s.dEdx());
    freco_shower_dEdxErr.push_back(s.dEdxErr());
    freco_shower_has_open_angle.push_back(s.OpenAngle());
    freco_shower_has_length.push_back(s.has_length());

    std::vector<int> & reco_shower_to_reco_hit = freco_shower_to_reco_hit.at(index);
    reco_shower_to_reco_hit.reserve(hits_per_shower.size());
    for(art::Ptr<recob::Hit> const & hit_ptr : hits_per_shower.at(i)) {
      reco_shower_to_reco_hit.push_back(hit_ptr.key() + hit_index_offset);
    }

    std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> const energy_pair = GetShowerHelperEnergy(e, shower_producer, i);
    freco_shower_EnergyHelper_energy_legacy.push_back(energy_pair.first.first.at(energy_pair.first.second));
    freco_shower_EnergyHelper_energy.push_back(energy_pair.first.first);
    freco_shower_EnergyHelper_dedx.push_back(energy_pair.second);

    std::pair<std::pair<std::vector<double>, size_t>, std::vector<double>> const energynew_pair = GetShowerHelperEnergy(e, shower_producer, i, true);
    freco_shower_EnergyHelperNew_energy_legacy.push_back(energynew_pair.first.first.at(energynew_pair.first.second));
    freco_shower_EnergyHelperNew_energy.push_back(energynew_pair.first.first);
    freco_shower_EnergyHelperNew_dedx.push_back(energynew_pair.second);

    if(shower_matches) {
      
      RecoMCMatch const & shower_match = shower_matches->at(i);

      freco_shower_largest_mc_type.push_back(shower_match.mc_type);
      freco_shower_largest_mc_index.push_back(shower_match.mc_index);
      freco_shower_largest_ratio.push_back(shower_match.ratio);
      
      std::vector<int> & reco_shower_mc_type = freco_shower_mc_type.at(index);
      std::vector<int> & reco_shower_mc_index = freco_shower_mc_index.at(index);
      std::vector<double> & reco_shower_charge_contribution = freco_shower_charge_contribution.at(index);    
      
      size_t const matching_map_size = shower_match.mctrack_map.size() + shower_match.mcshower_map.size() + shower_match.mcparticle_map.size();
      reco_shower_mc_type.reserve(matching_map_size);
      reco_shower_mc_index.reserve(matching_map_size);
      reco_shower_charge_contribution.reserve(matching_map_size);
      
      double charge_total = 0;
      
      for(auto const & p : shower_match.mctrack_map) {
	reco_shower_mc_type.push_back(fmc_type_track);
	reco_shower_mc_index.push_back(p.first);
	reco_shower_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }
      
      for(auto const & p : shower_match.mcshower_map) {
	reco_shower_mc_type.push_back(fmc_type_shower);
	reco_shower_mc_index.push_back(p.first);
	reco_shower_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }
      
      for(auto const & p : shower_match.mcparticle_map) {
	reco_shower_mc_type.push_back(fmc_type_particle);
	reco_shower_mc_index.push_back(p.first);
	reco_shower_charge_contribution.push_back(p.second);
	charge_total += p.second;
      }    
      
      freco_shower_charge_total.push_back(charge_total);

    }

  }  
  
}


void FillLightEvent::FillRecoShowerVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < fshowerp_size; ++i) {
    art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producers.at(i));
    size += ev_s->size();
    fshower_producer_indices.push_back(size);
  }

  freco_shower_producer_index.reserve(size);
  freco_shower_Direction_x.reserve(size);
  freco_shower_Direction_y.reserve(size);
  freco_shower_Direction_z.reserve(size);
  freco_shower_DirectionErr_x.reserve(size);
  freco_shower_DirectionErr_y.reserve(size);
  freco_shower_DirectionErr_z.reserve(size);
  freco_shower_ShowerStart_x.reserve(size);
  freco_shower_ShowerStart_y.reserve(size);
  freco_shower_ShowerStart_z.reserve(size);
  freco_shower_ShowerStartErr_x.reserve(size);
  freco_shower_ShowerStartErr_y.reserve(size);
  freco_shower_ShowerStartErr_z.reserve(size);
  freco_shower_Energy.reserve(size);
  freco_shower_EnergyErr.reserve(size);
  freco_shower_MIPEnergy.reserve(size);
  freco_shower_MIPEnergyErr.reserve(size);
  freco_shower_best_plane.reserve(size);
  freco_shower_Length.reserve(size);
  freco_shower_OpenAngle.reserve(size);
  freco_shower_dEdx.reserve(size);
  freco_shower_dEdxErr.reserve(size);
  freco_shower_has_open_angle.reserve(size);
  freco_shower_has_length.reserve(size);
  freco_shower_to_reco_hit.resize(size, std::vector<int>());
  freco_shower_EnergyHelper_energy_legacy.reserve(size);
  freco_shower_EnergyHelper_energy.reserve(size);
  freco_shower_EnergyHelper_dedx.reserve(size);
  freco_shower_EnergyHelperNew_energy_legacy.reserve(size);
  freco_shower_EnergyHelperNew_energy.reserve(size);
  freco_shower_EnergyHelperNew_dedx.reserve(size);
  if(frmcm_first) {
    freco_shower_largest_mc_type.reserve(size);
    freco_shower_largest_mc_index.reserve(size);
    freco_shower_largest_ratio.reserve(size);
    freco_shower_mc_type.resize(size, std::vector<int>());
    freco_shower_mc_index.resize(size, std::vector<int>());
    freco_shower_charge_contribution.resize(size, std::vector<double>());
    freco_shower_charge_total.reserve(size);
  }

  size_t shower_index_offset = 0;
  size_t hit_index_offset = 0;
  for(size_t i = 0; i < fshowerp_size; ++i) {
    FillRecoShowerVectors(e, i, shower_index_offset, hit_index_offset);
    art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producers.at(i));
    shower_index_offset += ev_s->size();
    art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(i));
    hit_index_offset += ev_h->size();
  }

}


void FillLightEvent::FillPandora(art::Event const & e) {

  std::string const fpfp_producer = "pandoraNu";
  art::ValidHandle<std::vector<recob::PFParticle>> const & ev_pfp = e.getValidHandle<std::vector<recob::PFParticle>>(fpfp_producer);
  size_t const size = ev_pfp->size();

  fpfp_pdg.reserve(size);
  fpfp_vertex_X.reserve(size);
  fpfp_vertex_Y.reserve(size);
  fpfp_vertex_Z.reserve(size);
  fpfp_original_index.reserve(size);
  fpfp_children.reserve(size);

  size_t track_offset = 0;
  for(size_t i = 0; i < ftrackp_size; ++i) {
    if(ftrack_producers.at(i) == "pandoraNu") break;
    art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producers.at(i));
    track_offset += ev_t->size();
  }

  size_t shower_offset = 0;
  for(size_t i = 0; i < fshowerp_size; ++i) {
    if(fshower_producers.at(i) == "pandoraNu") break;
    art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producers.at(i));
    shower_offset += ev_s->size();
  }

  art::FindManyP<recob::Vertex> PFPToVertex(ev_pfp, e, fpfp_producer);
  art::FindManyP<recob::Shower> PFPToShower(ev_pfp, e, fpfp_producer);
  art::FindManyP<recob::Track> PFPToTrack(ev_pfp, e, fpfp_producer);

  for(size_t i = 0; i < size; ++i) {

    recob::PFParticle const & pfp = ev_pfp->at(i);
    int const pdg = pfp.PdgCode();
    int original_index = -1;

    std::vector<art::Ptr<recob::Vertex>> asso_vertices = PFPToVertex.at(i);
    if(asso_vertices.size() > 1) {
      std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one vertex associated with pfp\n";
    }
    std::vector<double> xyz(3, -10000);

    if(abs(pdg) == 12 || abs(pdg) == 14) {
      asso_vertices.front()->XYZ(&(xyz[0]));
    }
    else if(abs(pdg) == 13) {
      std::vector<art::Ptr<recob::Track>> asso_tracks = PFPToTrack.at(i);
      if(asso_tracks.size() > 1) {
	std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one track associated with pfp\n";
      }
      else if(asso_tracks.size() == 0) {
	if(fverbose) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: no track associated with pfp\n";
      }
      else {
	original_index = asso_tracks.front().key() + track_offset;
      }
    }

    else if(abs(pdg) == 11) {
      std::vector<art::Ptr<recob::Shower>> asso_showers = PFPToShower.at(i);
      if(asso_showers.size() > 1) {
	std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: more than one shower associated with pfp\n";
      }
      else if(asso_showers.size() == 0) {
	if(fverbose) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nWarning: no shower associated with pfp\n";
      }
      else {
	original_index = asso_showers.front().key() + shower_offset;
      }
    }

    fpfp_pdg.push_back(pdg);
    fpfp_vertex_X.push_back(xyz[0]);
    fpfp_vertex_Y.push_back(xyz[1]);
    fpfp_vertex_Z.push_back(xyz[2]);
    fpfp_original_index.push_back(original_index);
    std::vector<int> daughters;
    daughters.reserve(pfp.Daughters().size());
    for(unsigned int const dtid : pfp.Daughters()) daughters.push_back(dtid);
    fpfp_children.push_back(daughters);
    
  }
  
}


void FillLightEvent::FillWeights(art::Event const & e) {
  if(fuse_eventweight){

  art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw =
    e.getValidHandle<std::vector<evwgh::MCEventWeight>>("eventweight");

  std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;

  if(ev_evw->size() > 1) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: eventweight has more than one entry\n";
  }

  for(auto const & p : weight_map) {
    auto const wbm_it = fweight_branch_map.find(p.first);
    if(wbm_it == fweight_branch_map.end()) {
      std::cout << "Could not find weight: " << p.first << "\n";
      continue;
    }
    *wbm_it->second.at(0) = p.second.at(0);
    *wbm_it->second.at(1) = p.second.at(1);
  }
}
}


void FillLightEvent::FillGenieParticleVectors(art::Event const & e,
					      size_t const mct_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  simb::MCTruth const & mct = ev_mct->at(mct_index);

  std::vector<int> & genie_particle_TrackId = fgenie_particle_TrackId.at(mct_index);
  std::vector<int> & genie_particle_StatusCode = fgenie_particle_StatusCode.at(mct_index);
  std::vector<int> & genie_particle_PdgCode = fgenie_particle_PdgCode.at(mct_index);
  std::vector<int> & genie_particle_Mother = fgenie_particle_Mother.at(mct_index);
  std::vector<double> & genie_particle_X = fgenie_particle_X.at(mct_index);
  std::vector<double> & genie_particle_Y = fgenie_particle_Y.at(mct_index);
  std::vector<double> & genie_particle_Z = fgenie_particle_Z.at(mct_index);
  std::vector<double> & genie_particle_T = fgenie_particle_T.at(mct_index);
  std::vector<double> & genie_particle_Px = fgenie_particle_Px.at(mct_index);
  std::vector<double> & genie_particle_Py = fgenie_particle_Py.at(mct_index);
  std::vector<double> & genie_particle_Pz = fgenie_particle_Pz.at(mct_index);
  std::vector<double> & genie_particle_E = fgenie_particle_E.at(mct_index);

  genie_particle_TrackId.reserve(mct.NParticles());
  genie_particle_StatusCode.reserve(mct.NParticles());
  genie_particle_PdgCode.reserve(mct.NParticles());
  genie_particle_Mother.reserve(mct.NParticles());
  genie_particle_X.reserve(mct.NParticles());
  genie_particle_Y.reserve(mct.NParticles());
  genie_particle_Z.reserve(mct.NParticles());
  genie_particle_T.reserve(mct.NParticles());
  genie_particle_Px.reserve(mct.NParticles());
  genie_particle_Py.reserve(mct.NParticles());
  genie_particle_Pz.reserve(mct.NParticles());
  genie_particle_E.reserve(mct.NParticles());

  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    genie_particle_TrackId.push_back(mcp.TrackId());
    genie_particle_StatusCode.push_back(mcp.StatusCode());
    genie_particle_PdgCode.push_back(mcp.PdgCode());
    genie_particle_Mother.push_back(mcp.Mother());
    simb::MCTrajectory const & mctraj = mcp.Trajectory();
    genie_particle_X.push_back(mctraj.Position(0).X());
    genie_particle_Y.push_back(mctraj.Position(0).Y());
    genie_particle_Z.push_back(mctraj.Position(0).Z());
    genie_particle_T.push_back(mctraj.Position(0).T());
    genie_particle_Px.push_back(mctraj.Momentum(0).Px());
    genie_particle_Py.push_back(mctraj.Momentum(0).Py());
    genie_particle_Pz.push_back(mctraj.Momentum(0).Pz());
    genie_particle_E.push_back(mctraj.Momentum(0).E());    

    //std::cout << mcp.TrackId() << " " << mcp.Mother() << " " << mcp.PdgCode() << " " << mctraj.Position(0).X() << " " << mctraj.Position(0).Y() << " " << mctraj.Position(0).Z() << "\n";
    
  }
  
}


void FillLightEvent::FillGenieParticleVectors(art::Event const & e) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  size_t const size = ev_mct->size();

  fgenie_particle_TrackId.resize(size, {});
  fgenie_particle_StatusCode.resize(size, {});
  fgenie_particle_PdgCode.resize(size, {});
  fgenie_particle_Mother.resize(size, {});
  fgenie_particle_X.resize(size, {});
  fgenie_particle_Y.resize(size, {});
  fgenie_particle_Z.resize(size, {});
  fgenie_particle_T.resize(size, {});
  fgenie_particle_Px.resize(size, {});
  fgenie_particle_Py.resize(size, {});
  fgenie_particle_Pz.resize(size, {});
  fgenie_particle_E.resize(size, {});
  
  for(size_t mct_index = 0; mct_index < size; ++mct_index) {
    FillGenieParticleVectors(e, mct_index);
  }

}


void FillLightEvent::FillMCParticleVectors(art::Event const & e) {

  art::ValidHandle<std::vector<simb::MCParticle>> const & ev_mcp = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  size_t const size = ev_mcp->size();

  fmcparticle_TrackId.reserve(size);
  fmcparticle_StatusCode.reserve(size);
  fmcparticle_PdgCode.reserve(size);
  fmcparticle_Mother.reserve(size);  
  /*
  fmcparticle_X.reserve(size);
  fmcparticle_Y.reserve(size);
  fmcparticle_Z.reserve(size);
  fmcparticle_T.reserve(size);
  fmcparticle_Px.reserve(size);
  fmcparticle_Py.reserve(size);
  fmcparticle_Pz.reserve(size);
  fmcparticle_E.reserve(size);
  */

  for(simb::MCParticle const & mcp : *ev_mcp) {
    fmcparticle_TrackId.push_back(mcp.TrackId());
    fmcparticle_StatusCode.push_back(mcp.StatusCode());
    fmcparticle_PdgCode.push_back(mcp.PdgCode());
    fmcparticle_Mother.push_back(mcp.Mother());
    /*
    simb::MCTrajectory const & mctraj = mcp.Trajectory();
    fmcparticle_X.push_back(mctraj.Position(0).X());
    fmcparticle_Y.push_back(mctraj.Position(0).Y());
    fmcparticle_Z.push_back(mctraj.Position(0).Z());
    fmcparticle_T.push_back(mctraj.Position(0).T());
    fmcparticle_Px.push_back(mctraj.Momentum(0).Px());
    fmcparticle_Py.push_back(mctraj.Momentum(0).Py());
    fmcparticle_Pz.push_back(mctraj.Momentum(0).Pz());
    fmcparticle_E.push_back(mctraj.Momentum(0).Pz());
    */
  }

}
  
  
void FillLightEvent::FillMCTrackVectors(art::Event const & e) {

  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctr = e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");  
  size_t const size = ev_mctr->size();

  fmctrack_Origin.reserve(size);
  fmctrack_PdgCode.reserve(size);
  fmctrack_TrackID.reserve(size);
  fmctrack_Process.reserve(size);
  fmctrack_Start_X.reserve(size);
  fmctrack_Start_Y.reserve(size);
  fmctrack_Start_Z.reserve(size);
  fmctrack_Start_T.reserve(size);
  fmctrack_Start_Px.reserve(size);
  fmctrack_Start_Py.reserve(size);
  fmctrack_Start_Pz.reserve(size);
  fmctrack_Start_E.reserve(size);
  fmctrack_End_X.reserve(size);
  fmctrack_End_Y.reserve(size);
  fmctrack_End_Z.reserve(size);
  fmctrack_End_T.reserve(size);
  fmctrack_End_Px.reserve(size);
  fmctrack_End_Py.reserve(size);
  fmctrack_End_Pz.reserve(size);
  fmctrack_End_E.reserve(size);
  fmctrack_X.resize(size, std::vector<double>());
  fmctrack_Y.resize(size, std::vector<double>());
  fmctrack_Z.resize(size, std::vector<double>());
  fmctrack_T.resize(size, std::vector<double>());
  fmctrack_Px.resize(size, std::vector<double>());
  fmctrack_Py.resize(size, std::vector<double>());
  fmctrack_Pz.resize(size, std::vector<double>());
  fmctrack_E.resize(size, std::vector<double>());
  fmctrack_dQdx.reserve(size);
  fmctrack_dEdx.reserve(size);
  fmctrack_MotherPdgCode.reserve(size);
  fmctrack_MotherTrackID.reserve(size);
  fmctrack_MotherProcess.reserve(size);
  fmctrack_AncestorPdgCode.reserve(size);
  fmctrack_AncestorTrackID.reserve(size);
  fmctrack_AncestorProcess.reserve(size);
  if(frmcm_first) fmctrack_contributed_charge.resize(size, std::vector<double>());  

  for(size_t i = 0; i < size; ++i) {

    sim::MCTrack const & mctr = ev_mctr->at(i);
    size_t const mctr_size = mctr.size();

    fmctrack_Origin.push_back(mctr.Origin());
    fmctrack_PdgCode.push_back(mctr.PdgCode());
    fmctrack_TrackID.push_back(mctr.TrackID());
    fmctrack_Process.push_back(mctr.Process());
    sim::MCStep const & start = mctr.Start();
    fmctrack_Start_X.push_back(start.X());
    fmctrack_Start_Y.push_back(start.Y());
    fmctrack_Start_Z.push_back(start.Z());
    fmctrack_Start_T.push_back(start.T());
    fmctrack_Start_Px.push_back(start.Px());
    fmctrack_Start_Py.push_back(start.Py());
    fmctrack_Start_Pz.push_back(start.Pz());
    fmctrack_Start_E.push_back(start.E());
    sim::MCStep const & end = mctr.End();
    fmctrack_End_X.push_back(end.X());
    fmctrack_End_Y.push_back(end.Y());
    fmctrack_End_Z.push_back(end.Z());
    fmctrack_End_T.push_back(end.T());
    fmctrack_End_Px.push_back(end.Px());
    fmctrack_End_Py.push_back(end.Py());
    fmctrack_End_Pz.push_back(end.Pz());
    fmctrack_End_E.push_back(end.E());
    
    std::vector<double> & mctrack_X_mcstp = fmctrack_X.at(i);
    std::vector<double> & mctrack_Y_mcstp = fmctrack_Y.at(i);
    std::vector<double> & mctrack_Z_mcstp = fmctrack_Z.at(i);
    std::vector<double> & mctrack_T_mcstp = fmctrack_T.at(i);
    std::vector<double> & mctrack_Px_mcstp = fmctrack_Px.at(i);
    std::vector<double> & mctrack_Py_mcstp = fmctrack_Py.at(i);
    std::vector<double> & mctrack_Pz_mcstp = fmctrack_Pz.at(i);
    std::vector<double> & mctrack_E_mcstp = fmctrack_E.at(i);

    mctrack_X_mcstp.reserve(mctr_size);
    mctrack_Y_mcstp.reserve(mctr_size);
    mctrack_Z_mcstp.reserve(mctr_size);
    mctrack_T_mcstp.reserve(mctr_size);
    mctrack_Px_mcstp.reserve(mctr_size);
    mctrack_Py_mcstp.reserve(mctr_size);
    mctrack_Pz_mcstp.reserve(mctr_size);
    mctrack_E_mcstp.reserve(mctr_size);

    for(sim::MCStep const & mcstp : mctr) {
      mctrack_X_mcstp.push_back(mcstp.X());
      mctrack_Y_mcstp.push_back(mcstp.Y());
      mctrack_Z_mcstp.push_back(mcstp.Z());
      mctrack_T_mcstp.push_back(mcstp.T());
      mctrack_Px_mcstp.push_back(mcstp.Px());
      mctrack_Py_mcstp.push_back(mcstp.Py());
      mctrack_Pz_mcstp.push_back(mcstp.Pz());
      mctrack_E_mcstp.push_back(mcstp.E());
    }

    fmctrack_dQdx.push_back(mctr.dQdx());
    fmctrack_dEdx.push_back(mctr.dEdx());
    fmctrack_MotherPdgCode.push_back(mctr.MotherPdgCode());
    fmctrack_MotherTrackID.push_back(mctr.MotherTrackID());
    fmctrack_MotherProcess.push_back(mctr.MotherProcess());
    fmctrack_AncestorPdgCode.push_back(mctr.AncestorPdgCode());
    fmctrack_AncestorTrackID.push_back(mctr.AncestorTrackID());
    fmctrack_AncestorProcess.push_back(mctr.AncestorProcess());

    if(frmcm_first) {

      std::vector<double> & mctrack_contributed_charge = fmctrack_contributed_charge.at(i);
      mctrack_contributed_charge.reserve(frmcm_size);
      
      for(size_t j = 0; j < frmcm_size; ++j) {
	RecoMCMatching const & rmcm = frmcm.at(j);
	mctrack_contributed_charge.push_back(rmcm.GetMCTrackCharge().at(i));
      }

    }
    
  }

}


void FillLightEvent::FillMCShowerVectors(art::Event const & e) {

  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcs = e.getValidHandle<std::vector<sim::MCShower>>("mcreco");  
  size_t const size = ev_mcs->size();

  fmcshower_Origin.reserve(size);
  fmcshower_PdgCode.reserve(size);
  fmcshower_TrackID.reserve(size);
  fmcshower_Process.reserve(size);
  fmcshower_Start_X.reserve(size);
  fmcshower_Start_Y.reserve(size);
  fmcshower_Start_Z.reserve(size);
  fmcshower_Start_T.reserve(size);
  fmcshower_Start_Px.reserve(size);
  fmcshower_Start_Py.reserve(size);
  fmcshower_Start_Pz.reserve(size);
  fmcshower_Start_E.reserve(size);
  fmcshower_End_X.reserve(size);
  fmcshower_End_Y.reserve(size);
  fmcshower_End_Z.reserve(size);
  fmcshower_End_T.reserve(size);
  fmcshower_End_Px.reserve(size);
  fmcshower_End_Py.reserve(size);
  fmcshower_End_Pz.reserve(size);
  fmcshower_End_E.reserve(size);
  fmcshower_MotherPdgCode.reserve(size);
  fmcshower_MotherTrackID.reserve(size);
  fmcshower_MotherProcess.reserve(size);
  fmcshower_AncestorPdgCode.reserve(size);
  fmcshower_AncestorTrackID.reserve(size);
  fmcshower_AncestorProcess.reserve(size);
  fmcshower_DetProfile_X.reserve(size);
  fmcshower_DetProfile_Y.reserve(size);
  fmcshower_DetProfile_Z.reserve(size);
  fmcshower_DetProfile_T.reserve(size);
  fmcshower_DetProfile_Px.reserve(size);
  fmcshower_DetProfile_Py.reserve(size);
  fmcshower_DetProfile_Pz.reserve(size);
  fmcshower_DetProfile_E.reserve(size);
  fmcshower_DaughterTrackID.reserve(size);
  fmcshower_Charge.reserve(size);
  fmcshower_dQdx.reserve(size);
  fmcshower_StartDir_X.reserve(size);
  fmcshower_StartDir_Y.reserve(size);
  fmcshower_StartDir_Z.reserve(size);
  if(frmcm_first) fmcshower_contributed_charge.resize(size, std::vector<double>());

  for(size_t i = 0; i < ev_mcs->size(); ++i) {

    sim::MCShower const & mcs = ev_mcs->at(i);

    fmcshower_Origin.push_back(mcs.Origin());
    fmcshower_PdgCode.push_back(mcs.PdgCode());
    fmcshower_TrackID.push_back(mcs.TrackID());
    fmcshower_Process.push_back(mcs.Process());
    sim::MCStep const & start = mcs.Start();
    fmcshower_Start_X.push_back(start.X());
    fmcshower_Start_Y.push_back(start.Y());
    fmcshower_Start_Z.push_back(start.Z());
    fmcshower_Start_T.push_back(start.T());
    fmcshower_Start_Px.push_back(start.Px());
    fmcshower_Start_Py.push_back(start.Py());
    fmcshower_Start_Pz.push_back(start.Pz());
    fmcshower_Start_E.push_back(start.E());
    sim::MCStep const & end = mcs.End();
    fmcshower_End_X.push_back(end.X());
    fmcshower_End_Y.push_back(end.Y());
    fmcshower_End_Z.push_back(end.Z());
    fmcshower_End_T.push_back(end.T());
    fmcshower_End_Px.push_back(end.Px());
    fmcshower_End_Py.push_back(end.Py());
    fmcshower_End_Pz.push_back(end.Pz());
    fmcshower_End_E.push_back(end.E());
    fmcshower_MotherPdgCode.push_back(mcs.MotherPdgCode());
    fmcshower_MotherTrackID.push_back(mcs.MotherTrackID());
    fmcshower_MotherProcess.push_back(mcs.MotherProcess());
    fmcshower_AncestorPdgCode.push_back(mcs.AncestorPdgCode());
    fmcshower_AncestorTrackID.push_back(mcs.AncestorTrackID());
    fmcshower_AncestorProcess.push_back(mcs.AncestorProcess());
    sim::MCStep const & detprofile = mcs.DetProfile();
    fmcshower_DetProfile_X.push_back(detprofile.X());
    fmcshower_DetProfile_Y.push_back(detprofile.Y());
    fmcshower_DetProfile_Z.push_back(detprofile.Z());
    fmcshower_DetProfile_T.push_back(detprofile.T());
    fmcshower_DetProfile_Px.push_back(detprofile.Px());
    fmcshower_DetProfile_Py.push_back(detprofile.Py());
    fmcshower_DetProfile_Pz.push_back(detprofile.Pz());
    fmcshower_DetProfile_E.push_back(detprofile.E());
    std::vector<int> mcs_DaughterTrackID;
    mcs_DaughterTrackID.reserve(mcs.DaughterTrackID().size());
    for(unsigned const int dtid : mcs.DaughterTrackID()) mcs_DaughterTrackID.push_back(dtid);
    fmcshower_DaughterTrackID.push_back(mcs_DaughterTrackID);
    fmcshower_Charge.push_back(mcs.Charge());
    fmcshower_dQdx.push_back(mcs.Charge());
    fmcshower_StartDir_X.push_back(mcs.StartDir().X());
    fmcshower_StartDir_Y.push_back(mcs.StartDir().Y());
    fmcshower_StartDir_Z.push_back(mcs.StartDir().Z());

    if(frmcm_first) {

      std::vector<double> & mcshower_contributed_charge = fmcshower_contributed_charge.at(i);
      mcshower_contributed_charge.reserve(fshowerp_size);
      
      for(size_t j = 0; j < frmcm_size; ++j) {
	RecoMCMatching const & rmcm = frmcm.at(j);
	mcshower_contributed_charge.push_back(rmcm.GetMCShowerCharge().at(i));
      }
      
    }

  }

}


void FillLightEvent::FillTruth(art::Event const & e) {
  
  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");  
  size_t const size = ev_mct->size();

  fnu_pdg.reserve(size);
  fnu_energy.reserve(size);
  flep_pdg.reserve(size);
  flep_energy.reserve(size);
  fccnc.reserve(size);
  fmode.reserve(size);
  finteraction_type.reserve(size);

  ftrue_nuvertx.reserve(size);
  ftrue_nuverty.reserve(size);
  ftrue_nuvertz.reserve(size);

  ftrue_nu_E.reserve(size);

  ftrue_nu_vtx_tpc_contained.reserve(size);
  ftrue_nu_vtx_fid_contained.reserve(size);
  
  fis_delta_rad.reserve(size);
  fdelta_index.reserve(size);
  fdelta_photon_index.reserve(size);
  fdelta_mcshower_index.reserve(size);
  fdelta_proton_index.reserve(size);
  fdelta_mctrack_index.reserve(size);
  
  size_t temp;
  for(size_t i = 0; i < ev_mct->size(); ++i) {

    if(ffs.RunSingle(e, i, temp)) {
      fdelta_mct_index = i;
      fis_delta_rad.push_back(1);
      GetDeltaMCShowerMCTrackIndices(e, i);
    }
    else {
      fis_delta_rad.push_back(0);
      fdelta_index.push_back(-1);
      fdelta_photon_index.push_back(-1);
      fdelta_mcshower_index.push_back(-1);
      fdelta_proton_index.push_back(-1);
      fdelta_mctrack_index.push_back(-1);      
    }

    simb::MCTruth const & mct = ev_mct->at(i);
    simb::MCNeutrino const & mcn = mct.GetNeutrino();

    fnu_pdg.push_back(mcn.Nu().PdgCode());
    fnu_energy.push_back(mcn.Nu().Trajectory().E(0));
    flep_pdg.push_back(mcn.Lepton().PdgCode());
    flep_energy.push_back(mcn.Lepton().Trajectory().E(0));
    fccnc.push_back(mcn.CCNC());
    fmode.push_back(mcn.Mode());
    finteraction_type.push_back(mcn.InteractionType());

    TLorentzVector const & true_nu_pos = mcn.Nu().Trajectory().Position(0);

    ftrue_nuvertx.push_back(true_nu_pos.X());
    ftrue_nuverty.push_back(true_nu_pos.Y());
    ftrue_nuvertz.push_back(true_nu_pos.Z());

    ftrue_nu_E.push_back(mcn.Nu().Trajectory().E(0));

    if(ftpc_volume.Contain(geoalgo::Point_t(true_nu_pos))) ftrue_nu_vtx_tpc_contained.push_back(1);
    else ftrue_nu_vtx_tpc_contained.push_back(0);
    if(ffiducial_volume.Contain(geoalgo::Point_t(true_nu_pos))) ftrue_nu_vtx_fid_contained.push_back(1);
    else ftrue_nu_vtx_fid_contained.push_back(0);
        
  }
  
}


void FillLightEvent::GetDeltaMCShowerMCTrackIndices(art::Event const & e,
						    size_t const delta_rad_mct_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");  
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack =
    e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower =
    e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

  simb::MCTruth const & mct = ev_mctruth->at(delta_rad_mct_index);
  
  int delta_index = -1;
  int delta_photon_index = -1;
  int delta_mcshower_index = -1;
  int delta_proton_index = -1;
  int delta_mctrack_index = -1;
    
  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(abs(mcp.PdgCode()) == 2214 || abs(mcp.PdgCode()) == 2114) {
      delta_index = i;
      break;
    }
  }

  if(delta_index == -1) {
    std::cout << "No delta\n";
    return;
  }

  std::vector<int> children_pos;
  std::vector<int> children_ext_pos;
  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(mcp.Mother() == mct.GetParticle(delta_index).TrackId()) {
      if(mcp.StatusCode() == 1) {
	children_ext_pos.push_back(i);  
      }
      else {
	children_pos.push_back(i);
      }
    }
  }

  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(mcp.StatusCode() == 1) {
      for(int const j : children_pos) {
	simb::MCParticle const & mcp_c = mct.GetParticle(j);
	if(mcp.Mother() == mcp_c.TrackId()) {
	  children_ext_pos.push_back(i);
	}
      }
    }
  }

  int photon_counter = 0;
  int proton_counter = 0;
  int iphoton_index = -1;
  int iproton_index = -1;
  for(int const i : children_ext_pos) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(mcp.PdgCode() == 22) {
      ++photon_counter;
      iphoton_index = i;
    }
    if(mcp.PdgCode() == 2212) {
      ++proton_counter;
      iproton_index = i;
    }
  }

  if(iphoton_index == -1) {
    std::cout << "could not find photon\n";
    return;
  }

  delta_photon_index = iphoton_index;
  delta_proton_index = iproton_index;

  double const diffd = 1e-10;
  double const diffp = 1e-10;
  double const diffE = 1e-2;

  std::vector<int> exiting_visible_particles;
  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(mcp.StatusCode() == 1 && 
       abs(mcp.PdgCode()) != 12 &&
       abs(mcp.PdgCode()) != 14 &&
       abs(mcp.PdgCode()) != 2112 &&
       abs(mcp.PdgCode()) != 111)
      exiting_visible_particles.push_back(i);
  }

  std::map<int, int> mcp_mcs;    
  for(size_t i = 0; i < ev_mcshower->size(); ++i) {
    sim::MCShower const & mcs = ev_mcshower->at(i);
    if(mcs.TrackID() != mcs.AncestorTrackID() || mcs.Origin() != 1) continue;
    double dist = 2000;
    double E = 1e12;
    int mcp_index = -1;
    geoalgo::Point_t mcsp(mcs.Start().Position());     
    for(int const i : exiting_visible_particles) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      if(mcp.PdgCode() != mcs.PdgCode()) continue;
      double tdist = mcsp.Dist(mcp.Trajectory().Position(0));
      double tE = fabs(mcs.Start().E() - mcp.Trajectory().E(0) * 1e3);
      if(tdist <= dist && tE < E && tE < diffE) {
	dist = tdist;
	E = tE;
	mcp_index = i;
      }
    }
    if(dist < diffd) {
      if(mcp_mcs.find(mcp_index) != mcp_mcs.end()) {
	auto mm_it = mcp_mcs.find(mcp_index);
	std::cout << "mcp_index already added, current mcsack trackid: "
		  << mcs.TrackID() << " previous mcsack trackid: "
		  << ev_mcshower->at(mm_it->second).TrackID() << "\nmcp_index: "
		  << mcp_index << " mcp trackid: " << mct.GetParticle(mcp_index).TrackId() << std::endl;
	exit(1);
      }
      mcp_mcs.emplace(mcp_index, i);
    }
  }

  auto const mcs_it = mcp_mcs.find(delta_photon_index);
  if(mcs_it != mcp_mcs.end()) {
    delta_mcshower_index = mcs_it->second;
  }

  std::map<int, int> mcp_mctr;
  for(size_t i = 0; i < ev_mctrack->size(); ++i) {
    sim::MCTrack const & mctr = ev_mctrack->at(i);
    if(mctr.TrackID() != mctr.AncestorTrackID() || mctr.Origin() != 1) continue;
    double mom = 2000;
    int mcp_index = -1;
    geoalgo::Point_t mctrp(mctr.Start().Position());
    for(int const i : exiting_visible_particles) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      if(mcp.PdgCode() != mctr.PdgCode()) continue;
      double tdist = mctrp.Dist(mcp.Trajectory().Position(0));
      double tmom = geoalgo::Vector_t(mctr.Start().Momentum()*1e-3).Dist(mcp.Trajectory().Momentum(0));
      if(tmom < mom && tdist < diffd) {
	mom = tmom;
	mcp_index = i;
      }
    }

    if(mom < diffp) {
      if(mcp_mctr.find(mcp_index) != mcp_mctr.end()) {
	auto mm_it = mcp_mctr.find(mcp_index);
	std::cout << "mcp_index already added, current mctrack trackid: "
		  << mctr.TrackID() << " previous mctrack trackid: "
		  << ev_mctrack->at(mm_it->second).TrackID() << "\nmcp_index: "
		  << mcp_index << " mcp trackid: " << mct.GetParticle(mcp_index).TrackId() << std::endl;
	std::cout << "Ancestor MCTrack list:\n";
	for(size_t i = 0; i < ev_mctrack->size(); ++i) {
	  sim::MCTrack const & mctr = ev_mctrack->at(i);
	  if(mctr.TrackID() != mctr.AncestorTrackID() || mctr.Origin() != 1) continue;
	  std::cout << i << " tid: " << mctr.TrackID() << " pdg: " << mctr.PdgCode() << " origin: " << mctr.Origin() << " Pos: " << geoalgo::Point_t(mctr.Start().Position()) << " Mom: " << geoalgo::Point_t(mctr.Start().Momentum()*1e-3) << "\n";
	}
	std::cout << "Gen Origin: " << mct.Origin() << "\nMCParticles:\n";
	for(int const i : exiting_visible_particles) {
	  simb::MCParticle const & mcp = mct.GetParticle(i);
	  std::cout << i << " " << mcp.TrackId() << " " << mcp.PdgCode() << " Pos: " << geoalgo::Point_t(mcp.Trajectory().Position(0)) << " Mom: " << geoalgo::Point_t(mcp.Trajectory().Momentum(0)) << "\n";
	}
	exit(1);
      }
      mcp_mctr.emplace(mcp_index, i);
    }
  }

  if(fverbose && exiting_visible_particles.size() != mcp_mctr.size() + mcp_mcs.size()) {
    std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
	      << "WARNING: not all particles matched\n";
    std::cout << "Ancestor MCShowers\n";
    for(size_t i = 0; i < ev_mcshower->size(); ++i) {
      sim::MCShower const & mcs = ev_mcshower->at(i);
      if(mcs.TrackID() != mcs.AncestorTrackID() || mcs.Origin() != 1) continue;
      std::cout << i << " tid: " << mcs.TrackID() << " pdg: " << mcs.PdgCode() << " origin: " << mcs.Origin() << " Pos: " << geoalgo::Point_t(mcs.Start().Position()) << " Mom: " << geoalgo::Point_t(mcs.Start().Momentum()*1e-3) << "\n";
    }
    std::cout << "Ancestor MCTracks\n";
    for(size_t i = 0; i < ev_mctrack->size(); ++i) {
      sim::MCTrack const & mctr = ev_mctrack->at(i);
      if(mctr.TrackID() != mctr.AncestorTrackID() || mctr.Origin() != 1) continue;
      std::cout << i << " tid: " << mctr.TrackID() << " pdg: " << mctr.PdgCode() << " origin: " << mctr.Origin() << " Pos: " << geoalgo::Point_t(mctr.Start().Position()) << " Mom: " << geoalgo::Point_t(mctr.Start().Momentum()*1e-3) << "\n";
    }
    std::cout << "Gen Origin: " << mct.Origin() << "\nExiting MCParticles\n";
    for(int const i : exiting_visible_particles) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout << i << " " << mcp.TrackId() << " " << mcp.PdgCode() << " Pos: " << geoalgo::Point_t(mcp.Trajectory().Position(0)) << " Mom: " << geoalgo::Point_t(mcp.Trajectory().Momentum(0)) << " Energy: " << mcp.Trajectory().E(0) << "\n";
    }
    std::cout << "Shower matches\n";
    for(auto const & p : mcp_mcs) {
      std::cout << p.first << " " << p.second << "\n";
    }
    std::cout << "Track matches\n";
    for(auto const & p : mcp_mctr) {
      std::cout << p.first << " " << p.second << "\n";
    }

  }

  auto const mctr_it = mcp_mctr.find(delta_proton_index);
  if(mctr_it != mcp_mctr.end()) {
    delta_mctrack_index = mctr_it->second;    
  }

  fdelta_index.push_back(delta_index);
  fdelta_photon_index.push_back(delta_photon_index);
  fdelta_mcshower_index.push_back(delta_mcshower_index);
  fdelta_proton_index.push_back(delta_proton_index);
  fdelta_mctrack_index.push_back(delta_mctrack_index);

}


void FillLightEvent::analyze(art::Event const & e) {

  ++fnumber_of_events;

  ResetEvent();

  frun_number = e.id().run();
  fsubrun_number = e.id().subRun();
  fevent_number = e.id().event();

  if(fmc) for(RecoMCMatching & rmcm : frmcm) rmcm.MatchWAssociations(e);

  if(fverbose) std::cout << "swtrigger\n";
  FillSWTriggerVectors(e);
  if(fverbose) std::cout << "opflash\n";
  FillRecoOpFlashVectors(e);
  if(fverbose) std::cout << "hit\n";
  FillRecoHitVectors(e);
  if(fverbose) std::cout << "track\n";
  FillRecoTrackVectors(e);
  if(fverbose) std::cout << "shower\n";
  FillRecoShowerVectors(e);
  if(fverbose) std::cout << "pandora\n";
  FillPandora(e);
  if(fverbose) std::cout << "reco done\n";

  if(fmc) {

    art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");  

    FillWeights(e);

    if(ev_mct->front().GetNeutrino().Nu().Mother() == -1) FillTruth(e);

    FillGenieParticleVectors(e);
    FillMCParticleVectors(e);
    FillMCTrackVectors(e);
    FillMCShowerVectors(e);

  }  
  
  fevent_tree->Fill();

}


void FillLightEvent::endJob() {

  if(fpot_producer != "") {
    fpot_tree->Fill();
  }

}



DEFINE_ART_MODULE(FillLightEvent)
