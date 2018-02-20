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

  void ResetEvent();

  void FillRecoHitVectors(art::Event const & e);
  void FillRecoHitVectors(art::Event const & e, int const producer_index);
  void FillRecoTrackVectors(art::Event const & e);
  void FillRecoTrackVectors(art::Event const & e, int const producer_index);

  void FillWeights(art::Event const & e);

  void FillGenieParticleVectors(simb::MCTruth const & mct);
  void FillMCTrackVectors(art::Event const & e);
  void FillMCShowerVectors(art::Event const & e);
  void FillTruth(art::Event const & e,
		 size_t & delta_rad_mct_index);

  void GetDeltaMCShowerMCTrackIndices(art::Event const & e,
				      size_t const delta_rad_mct_index,
				      size_t & delta_photon_index,
				      size_t & delta_mcshower_index,
				      size_t & delta_proton_index,
				      size_t & delta_mctrack_index);
  
  void analyze(art::Event const & e) override;

private:

  bool fmc;
  std::string fpot_producer;

  std::vector<std::string> fhit_producers;
  std::vector<std::string> frmcm_hit_producers;
  std::vector<std::string> ftrack_producers;
  std::vector<std::string> fshower_producers;

  std::string fswtrigger_product;
  std::string fopflash_producer;

  std::vector<std::string> frmcmassociation_producers;
  size_t frmcm_size;

  geoalgo::AABox ftpc_volume;
  double foffset;
  geoalgo::AABox ffiducial_volume;
  geoalgo::GeoAlgo const falgo;

  std::vector<RecoMCMatching> frmcm;
  int fmc_type_shower;
  int fmc_type_track;
  int fmc_type_particle;

  FilterSignal ffs;
  lee::EnergyHelper energyHelper;

  bool fverbose;

  TTree * fpot_tree;
  int fnumber_of_events;
  double fpot;

  TTree * fevent_tree;

  //All
  int frun_number;
  int fsubrun_number;
  int fevent_number;

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
  std::vector<unsigned int> freco_hit_View;
  std::vector<unsigned int> freco_hit_SignalType;
  std::vector<unsigned int> freco_hit_WireID_CryostatID;
  std::vector<unsigned int> freco_hit_WireID_TPCID;
  std::vector<unsigned int> freco_hit_WireID_PlaneID;
  std::vector<unsigned int> freco_hit_WireID_WireID;

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
  std::vector<std::vector<double>> freco_track_Theta;
  std::vector<std::vector<double>> freco_track_Phi;
  std::vector<std::vector<double>> freco_track_ZenithAngle;
  std::vector<std::vector<double>> freco_track_AzimuthAngle;
  //Reco - MC matching
  std::vector<std::vector<int>> freco_track_mc_type;
  std::vector<std::vector<unsigned int>> freco_track_mc_index;
  std::vector<std::vector<double>> freco_track_charge_contribution;
  std::vector<double> freco_track_charge_total;

  //Flash
  int fpassed_swtrigger;
  double ftotalpe_sum;  
  double ftotalpe_ibg_sum;
  double ftotalpe_bbg_sum;

  //Truth
  int fnu_pdg;
  double fnu_energy;
  int flep_pdg;
  double flep_energy;
  int fccnc;
  int fmode;
  int finteraction_type;

  double ftrue_nu_E;
  double ftrue_nuvertx;
  double ftrue_nuverty;
  double ftrue_nuvertz;

  int ftrue_nu_vtx_tpc_contained;
  int ftrue_nu_vtx_fid_contained;

  //GENIE MCParticle
  std::vector<int> fgenie_particle_TrackId;
  std::vector<int> fgenie_particle_StatusCode;
  std::vector<int> fgenie_particle_PdgCode;
  std::vector<int> fgenie_particle_Mother;
  std::vector<double> fgenie_particle_X;
  std::vector<double> fgenie_particle_Y;
  std::vector<double> fgenie_particle_Z;
  std::vector<double> fgenie_particle_T;
  std::vector<double> fgenie_particle_Px;
  std::vector<double> fgenie_particle_Py;
  std::vector<double> fgenie_particle_Pz;
  std::vector<double> fgenie_particle_E;

  //MCTrack
  std::vector<int> fmctrack_Origin;
  std::vector<int> fmctrack_PdgCode;
  std::vector<int> fmctrack_TrackID;
  std::vector<std::string> fmctrack_Process;
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
  std::vector<std::vector<unsigned int>> fmcshower_DaughterTrackID;
  std::vector<std::vector<double>> fmcshower_Charge;
  std::vector<std::vector<double>> fmcshower_dQdx;
  std::vector<double> fmcshower_StartDir_X;
  std::vector<double> fmcshower_StartDir_Y;
  std::vector<double> fmcshower_StartDir_Z;
  std::vector<std::vector<double>> fmcshower_contributed_charge;

  //Delta radiative
  int fis_delta_rad;
  int fdelta_true_pdg;
  double fdelta_true_energy;
  int fdelta_photon_index;
  int fdelta_mcshower_index;
  int fdelta_proton_index;
  int fdelta_mctrack_index;
  double fdelta_photon_energy;
  double fdelta_proton_energy;
  int fdelta_mcshower_true_pdg;
  double fdelta_mcshower_true_energy;
  double fdelta_mcshower_detprofile_energy;
  int fdelta_mctrack_true_pdg;
  double fdelta_mctrack_true_energy;
  double fdelta_mctrack_true_length;

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
  fevent_tree(nullptr) 
{

  Reconfigure(p);
  SetupTrees();
  
}


void FillLightEvent::Reconfigure(fhicl::ParameterSet const & p) {

  fmc = p.get<bool>("mc");

  p.get_if_present<std::string>("pot_producer", fpot_producer);

  ftrack_producers = p.get<std::vector<std::string>>("track_producers");
  fshower_producers = p.get<std::vector<std::string>>("shower_producers");

  fopflash_producer = p.get<std::string>("opflash_producer");
  fswtrigger_product = p.get<std::string>("trigger_product");

  p.get_if_present<std::vector<std::string>>("hit_producers", fhit_producers);
  p.get_if_present<std::vector<std::string>>("rmcm_hit_producers", frmcm_hit_producers);
  p.get_if_present<std::vector<std::string>>("rmcmassociation_producers", frmcmassociation_producers);

  //Split rmcm into track and shower sections for some reason, need to seperate hit producer and hit producer for rmcm association

  frmcm_size = frmcmassociation_producers.size();
  if(!frmcmassociation_producers.empty()) {
    if(frmcm_size != frmcm_hit_producers.size() ||
       frmcm_size != ftrack_producers.size() ||
       frmcm_size != fshower_producers.size()) {
      std::cout << "Must have same number of producers\n";
      exit(1);
    }
    frmcm.reserve(frmcm_size);
    for(size_t i = 0; i < frmcm_size; ++i) {
      RecoMCMatching rmcm;
      rmcm.Configure(frmcm_hit_producers.at(i),
		     ftrack_producers.at(i),
		     fshower_producers.at(i),
		     frmcmassociation_producers.at(i));
      frmcm.push_back(rmcm);
    }
    fmc_type_shower = frmcm.front().fmc_type_shower;
    fmc_type_track = frmcm.front().fmc_type_track;
    fmc_type_particle = frmcm.front().fmc_type_particle;
  }

}


void FillLightEvent::SetupTrees() {
  
  art::ServiceHandle< art::TFileService > tfs;

  TTree * producer_tree = tfs->make<TTree>("producer_tree", "");
  producer_tree->Branch("track_producers", &ftrack_producers);
  producer_tree->Branch("shower_producers", &fshower_producers);
  producer_tree->Fill();

  if(fpot_producer != "") {
    fpot_tree = tfs->make<TTree>("pot_tree", "");
    fpot_tree->Branch("number_of_events", &fnumber_of_events, "number_of_events/I");
    fpot_tree->Branch("pot", &fpot, "pot/D");
  }

  fevent_tree = tfs->make<TTree>("event_tree", "");
  
  fevent_tree->Branch("run_number", &frun_number, "run_number/I");
  fevent_tree->Branch("subrun_number", &fsubrun_number, "subrun_number/I"); 
  fevent_tree->Branch("event_number", &fevent_number, "event_number/I");

  fevent_tree->Branch("reco_hit_producer_index", &freco_hit_producer_index);
  fevent_tree->Branch("reco_hit_StartTick", &freco_hit_StartTick);
  fevent_tree->Branch("reco_hit_EndTick", &freco_hit_EndTick);
  fevent_tree->Branch("reco_hit_PeakTime", &freco_hit_PeakTime);
  fevent_tree->Branch("reco_hit_SigmaPeakTime", &freco_hit_SigmaPeakTime);
  fevent_tree->Branch("reco_hit_RMS", &freco_hit_RMS);
  fevent_tree->Branch("reco_hit_PeakAmplitude", &freco_hit_PeakAmplitude);
  fevent_tree->Branch("reco_hit_SigmaPeakAmplitude", &freco_hit_SigmaPeakAmplitude);
  fevent_tree->Branch("reco_hit_SummedADC", &freco_hit_SummedADC);
  fevent_tree->Branch("reco_hit_Integral", &freco_hit_Integral);
  fevent_tree->Branch("reco_hit_SigmaIntegral", &freco_hit_SigmaIntegral);
  fevent_tree->Branch("reco_hit_Multiplicity", &freco_hit_Multiplicity);
  fevent_tree->Branch("reco_hit_LocalIndex", &freco_hit_LocalIndex);
  fevent_tree->Branch("reco_hit_GoodnessOfFit", &freco_hit_GoodnessOfFit);
  fevent_tree->Branch("reco_hit_DegreesOfFreedom", &freco_hit_DegreesOfFreedom);
  fevent_tree->Branch("reco_hit_View", &freco_hit_View);
  fevent_tree->Branch("reco_hit_SignalType", &freco_hit_SignalType);
  fevent_tree->Branch("reco_hit_WireID_CryostatID", &freco_hit_WireID_CryostatID);
  fevent_tree->Branch("reco_hit_WireID_TPCID", &freco_hit_WireID_TPCID);
  fevent_tree->Branch("reco_hit_WireID_PlaneID", &freco_hit_WireID_PlaneID);
  fevent_tree->Branch("reco_hit_WireID_WireID", &freco_hit_WireID_WireID);

  fevent_tree->Branch("reco_track_producer_index", &freco_track_producer_index);
  fevent_tree->Branch("reco_track_NumberTrajectoryPoints", &freco_track_NumberTrajectoryPoints);
  fevent_tree->Branch("reco_track_NPoints", &freco_track_NPoints);
  fevent_tree->Branch("reco_track_FirstPoint", &freco_track_FirstPoint);
  fevent_tree->Branch("reco_track_LastPoint", &freco_track_LastPoint);
  fevent_tree->Branch("reco_track_FirstValidPoint", &freco_track_FirstValidPoint);
  fevent_tree->Branch("reco_track_LastValidPoint", &freco_track_LastValidPoint);
  fevent_tree->Branch("reco_track_CountValidPoints", &freco_track_CountValidPoints);
  fevent_tree->Branch("reco_track_X", &freco_track_X);
  fevent_tree->Branch("reco_track_Y", &freco_track_Y);
  fevent_tree->Branch("reco_track_Z", &freco_track_Z);
  fevent_tree->Branch("reco_track_Px", &freco_track_Px);
  fevent_tree->Branch("reco_track_Py", &freco_track_Py);
  fevent_tree->Branch("reco_track_Pz", &freco_track_Pz);
  fevent_tree->Branch("reco_track_HasMomentum", &freco_track_HasMomentum);
  fevent_tree->Branch("reco_track_Length", &freco_track_Length);
  fevent_tree->Branch("reco_track_Chi2", &freco_track_Chi2);
  fevent_tree->Branch("reco_track_Chi2PerNdof", &freco_track_Chi2PerNdof);
  fevent_tree->Branch("reco_track_Ndof", &freco_track_Ndof);
  fevent_tree->Branch("reco_track_ParticleId", &freco_track_ParticleId);
  fevent_tree->Branch("reco_track_Theta", &freco_track_Theta);
  fevent_tree->Branch("reco_track_Phi", &freco_track_Phi);
  fevent_tree->Branch("reco_track_ZenithAngle", &freco_track_ZenithAngle);
  fevent_tree->Branch("reco_track_AzimuthAngle", &freco_track_AzimuthAngle);
  fevent_tree->Branch("reco_track_mc_type", &freco_track_mc_type);
  fevent_tree->Branch("reco_track_mc_index", &freco_track_mc_index);
  fevent_tree->Branch("reco_track_charge_contribution", &freco_track_charge_contribution);
  fevent_tree->Branch("reco_track_charge_total", &freco_track_charge_total);

  fevent_tree->Branch("passed_swtrigger", &fpassed_swtrigger, "passed_swtrigger/I");
  fevent_tree->Branch("totalpe_sum", &ftotalpe_sum, "totalpe_sum/D");
  fevent_tree->Branch("totalpe_ibg_sum", &ftotalpe_ibg_sum, "totalpe_ibg_sum/D");
  fevent_tree->Branch("totalpe_bbg_sum", &ftotalpe_bbg_sum, "totalpe_bbg_sum/D");

  if(fmc) {

    fevent_tree->Branch("nu_pdg", &fnu_pdg, "nu_pdg/I");
    fevent_tree->Branch("nu_energy", &fnu_energy, "nu_energy/D");
    fevent_tree->Branch("lep_pdg", &flep_pdg, "lep_pdg/I");
    fevent_tree->Branch("lep_energy", &flep_energy, "lep_energy/D");
    fevent_tree->Branch("ccnc", &fccnc, "ccnc/I");
    fevent_tree->Branch("mode", &fmode, "mode/I");
    fevent_tree->Branch("interaction_type", &finteraction_type, "interaction_type/I");

    fevent_tree->Branch("true_nu_E", &ftrue_nu_E, "true_nu_E/D");

    fevent_tree->Branch("true_nuvertx", &ftrue_nuvertx, "true_nuvertx/D");
    fevent_tree->Branch("true_nuverty", &ftrue_nuverty, "true_nuverty/D");
    fevent_tree->Branch("true_nuvertz", &ftrue_nuvertz, "true_nuvertz/D");

    fevent_tree->Branch("true_nu_vtx_tpc_contained", &ftrue_nu_vtx_tpc_contained, "true_nu_vtx_tpc_contained/I"); 
    fevent_tree->Branch("true_nu_vtx_fid_contained", &ftrue_nu_vtx_fid_contained, "true_nu_vtx_fid_contained/I"); 

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

    fevent_tree->Branch("mctrack_Origin", &fmctrack_Origin);
    fevent_tree->Branch("mctrack_PdgCode", &fmctrack_PdgCode);
    fevent_tree->Branch("mctrack_TrackID", &fmctrack_TrackID);
    fevent_tree->Branch("mctrack_Process", &fmctrack_Process);
    fevent_tree->Branch("mctrack_X", &fmctrack_X);
    fevent_tree->Branch("mctrack_Y", &fmctrack_Y);
    fevent_tree->Branch("mctrack_Z", &fmctrack_Z);
    fevent_tree->Branch("mctrack_T", &fmctrack_T);
    fevent_tree->Branch("mctrack_Px", &fmctrack_Px);
    fevent_tree->Branch("mctrack_Py", &fmctrack_Py);
    fevent_tree->Branch("mctrack_Pz", &fmctrack_Pz);
    fevent_tree->Branch("mctrack_E", &fmctrack_E);
    fevent_tree->Branch("mctrack_dQdx", &fmctrack_dQdx);
    fevent_tree->Branch("mctrack_dEdx", &fmctrack_dEdx);
    fevent_tree->Branch("mctrack_MotherTrackID", &fmctrack_MotherTrackID);
    fevent_tree->Branch("mctrack_MotherPdgCode", &fmctrack_MotherPdgCode);
    fevent_tree->Branch("mctrack_MotherProcess", &fmctrack_MotherProcess);
    fevent_tree->Branch("mctrack_AncestorTrackID", &fmctrack_AncestorTrackID);
    fevent_tree->Branch("mctrack_AncestorPdgCode", &fmctrack_AncestorPdgCode);
    fevent_tree->Branch("mctrack_AncestorProcess", &fmctrack_AncestorProcess);
    fevent_tree->Branch("mctrack_contributed_charge", &fmctrack_contributed_charge);

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
    fevent_tree->Branch("mcshower_Charge", &fmcshower_Charge);
    fevent_tree->Branch("mcshower_dQdx", &fmcshower_dQdx);
    fevent_tree->Branch("mcshower_StartDir_X", &fmcshower_StartDir_X);
    fevent_tree->Branch("mcshower_StartDir_Y", &fmcshower_StartDir_Y);
    fevent_tree->Branch("mcshower_StartDir_Z", &fmcshower_StartDir_Z);
    fevent_tree->Branch("mcshower_contributed_charge", &fmcshower_contributed_charge);

    fevent_tree->Branch("is_delta_rad", &fis_delta_rad, "is_delta_rad/I");
    fevent_tree->Branch("delta_true_pdg", &fdelta_true_pdg, "delta_true_pdg/I");
    fevent_tree->Branch("delta_true_energy", &fdelta_true_energy, "delta_true_energy/D");
    fevent_tree->Branch("delta_photon_index", &fdelta_photon_index, "delta_photon_index");
    fevent_tree->Branch("delta_mcshower_index", &fdelta_mcshower_index, "delta_mcshower_index");
    fevent_tree->Branch("delta_proton_index", &fdelta_proton_index, "delta_proton_index");
    fevent_tree->Branch("delta_mctrack_index", &fdelta_mctrack_index, "delta_mctrack_index");
    fevent_tree->Branch("delta_photon_energy", &fdelta_photon_energy, "delta_photon_energy/D");
    fevent_tree->Branch("delta_proton_energy", &fdelta_proton_energy, "delta_proton_energy/D");
    fevent_tree->Branch("delta_mcshower_true_pdg", &fdelta_mcshower_true_pdg, "delta_mcshower_true_pdg/I");
    fevent_tree->Branch("delta_mcshower_true_energy", &fdelta_mcshower_true_energy, "delta_mcshower_true_energy/D");
    fevent_tree->Branch("delta_mcshower_detprofile_energy", &fdelta_mcshower_detprofile_energy, "delta_mcshower_detprofile_energy/D");
    fevent_tree->Branch("delta_mctrack_true_pdg", &fdelta_mctrack_true_pdg, "delta_mctrack_true_pdg/I");
    fevent_tree->Branch("delta_mctrack_true_energy", &fdelta_mctrack_true_energy, "delta_mctrack_true_energy/D");
    fevent_tree->Branch("delta_mctrack_true_length", &fdelta_mctrack_true_length, "delta_mctrack_true_length/D");

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

  frun_number = -1;
  fsubrun_number = -1;
  fevent_number = -1;

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
  freco_track_mc_type.clear();
  freco_track_mc_index.clear();
  freco_track_charge_contribution.clear();
  freco_track_charge_total.clear();

  fpassed_swtrigger = -1;
  ftotalpe_sum = -1;
  ftotalpe_ibg_sum = -1;
  ftotalpe_bbg_sum = -1;

  fnu_pdg = 0;
  fnu_energy = -1;
  flep_pdg = -1;
  flep_energy = 0;
  fccnc = -1;
  fmode = -1;
  finteraction_type = -1;

  ftrue_nuvertx = -10000;
  ftrue_nuverty = -10000;
  ftrue_nuvertz = -10000;

  ftrue_nu_E = -1;

  ftrue_nu_vtx_tpc_contained = -1;
  ftrue_nu_vtx_fid_contained = -1;

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

  fmctrack_Origin.clear();
  fmctrack_PdgCode.clear();
  fmctrack_TrackID.clear();
  fmctrack_Process.clear();
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

  fis_delta_rad = -1;
  fdelta_true_pdg = 0;
  fdelta_true_energy = -1;
  fdelta_photon_index = -1;
  fdelta_mcshower_index = -1;
  fdelta_proton_index = -1;
  fdelta_mctrack_index = -1;
  fdelta_photon_energy = -1;
  fdelta_proton_energy = -1;
  fdelta_mcshower_true_pdg = 0;
  fdelta_mcshower_true_energy = -1;
  fdelta_mcshower_detprofile_energy = -1;
  fdelta_mctrack_true_pdg = 0;
  fdelta_mctrack_true_energy = -1;
  fdelta_mctrack_true_length = -1;  

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


void FillLightEvent::FillWeights(art::Event const & e) {

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


void FillLightEvent::FillRecoHitVectors(art::Event const & e, int const producer_index) {

  std::cout << fhit_producers.at(producer_index) << "\n";
  art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(producer_index));  

  for(recob::Hit const & h : *ev_h) {

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
    
  }

}


void FillLightEvent::FillRecoHitVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < fhit_producers.size(); ++i) {
    std::cout << fhit_producers.at(i) << "\n";
    art::ValidHandle<std::vector<recob::Hit>> const & ev_h = e.getValidHandle<std::vector<recob::Hit>>(fhit_producers.at(i));
    size += ev_h->size();
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

  for(size_t i = 0; i < fhit_producers.size(); ++i) {
    FillRecoHitVectors(e, i);
  }  

}


void FillLightEvent::FillRecoTrackVectors(art::Event const & e, int const producer_index) {

  art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producers.at(producer_index));  
  RecoMCMatching const & rmcm = frmcm.at(producer_index);
  std::vector<RecoMCMatch> const & track_matches = rmcm.GetTrackMatches();

  for(size_t i = 0; i < ev_t->size(); ++i) {

    recob::Track const & t = ev_t->at(i);
    size_t const traj_size = t.NumberTrajectoryPoints();
    RecoMCMatch const & track_match = track_matches.at(i);

    freco_track_producer_index.push_back(i);
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

    std::vector<double> & reco_track_X_traj = freco_track_X.at(i);
    std::vector<double> & reco_track_Y_traj = freco_track_Y.at(i);
    std::vector<double> & reco_track_Z_traj = freco_track_Z.at(i);
    std::vector<double> & reco_track_Px_traj = freco_track_Px.at(i);
    std::vector<double> & reco_track_Py_traj = freco_track_Py.at(i);
    std::vector<double> & reco_track_Pz_traj = freco_track_Pz.at(i);
    std::vector<double> & reco_track_Length_traj = freco_track_Length.at(i);
    std::vector<double> & reco_track_Theta_traj = freco_track_Theta.at(i);
    std::vector<double> & reco_track_Phi_traj = freco_track_Phi.at(i);
    std::vector<double> & reco_track_ZenithAngle_traj = freco_track_ZenithAngle.at(i);
    std::vector<double> & reco_track_AzimuthAngle_traj = freco_track_AzimuthAngle.at(i);    

    reco_track_X_traj.reserve(traj_size);
    reco_track_Y_traj.reserve(traj_size);
    reco_track_Z_traj.reserve(traj_size);
    reco_track_Px_traj.reserve(traj_size);
    reco_track_Py_traj.reserve(traj_size);
    reco_track_Pz_traj.reserve(traj_size);
    reco_track_Length_traj.reserve(traj_size);
    reco_track_Theta_traj.reserve(traj_size);
    reco_track_Phi_traj.reserve(traj_size);
    reco_track_ZenithAngle_traj.reserve(traj_size);
    reco_track_AzimuthAngle_traj.reserve(traj_size);
    
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
      reco_track_Theta_traj.push_back(t.Theta(j));
      reco_track_Phi_traj.push_back(t.Phi(j));
      reco_track_ZenithAngle_traj.push_back(t.ZenithAngle(j));
      reco_track_AzimuthAngle_traj.push_back(t.AzimuthAngle(j));

    }

    std::vector<int> & reco_track_mc_type = freco_track_mc_type.at(i);
    std::vector<unsigned int> & reco_track_mc_index = freco_track_mc_index.at(i);
    std::vector<double> & reco_track_charge_contribution = freco_track_charge_contribution.at(i);    

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


void FillLightEvent::FillRecoTrackVectors(art::Event const & e) {

  size_t size = 0;
  for(size_t i = 0; i < ftrack_producers.size(); ++i) {
    art::ValidHandle<std::vector<recob::Track>> const & ev_t = e.getValidHandle<std::vector<recob::Track>>(ftrack_producers.at(i));
    size += ev_t->size();
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
  freco_track_Theta.resize(size, std::vector<double>());
  freco_track_Phi.resize(size, std::vector<double>());
  freco_track_ZenithAngle.resize(size, std::vector<double>());
  freco_track_AzimuthAngle.resize(size, std::vector<double>());
  freco_track_mc_type.resize(size, std::vector<int>());
  freco_track_mc_index.resize(size, std::vector<unsigned int>());
  freco_track_charge_contribution.resize(size, std::vector<double>());
  freco_track_charge_total.reserve(size);

  for(size_t i = 0; i < ftrack_producers.size(); ++i) {
    FillRecoTrackVectors(e, i);
  }

}


void FillLightEvent::FillGenieParticleVectors(simb::MCTruth const & mct) {

  fgenie_particle_TrackId.reserve(mct.NParticles());
  fgenie_particle_StatusCode.reserve(mct.NParticles());
  fgenie_particle_PdgCode.reserve(mct.NParticles());
  fgenie_particle_Mother.reserve(mct.NParticles());
  fgenie_particle_X.reserve(mct.NParticles());
  fgenie_particle_Y.reserve(mct.NParticles());
  fgenie_particle_Z.reserve(mct.NParticles());
  fgenie_particle_T.reserve(mct.NParticles());
  fgenie_particle_Px.reserve(mct.NParticles());
  fgenie_particle_Py.reserve(mct.NParticles());
  fgenie_particle_Pz.reserve(mct.NParticles());
  fgenie_particle_E.reserve(mct.NParticles());

  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    fgenie_particle_TrackId.push_back(mcp.TrackId());
    fgenie_particle_StatusCode.push_back(mcp.StatusCode());
    fgenie_particle_PdgCode.push_back(mcp.PdgCode());
    fgenie_particle_Mother.push_back(mcp.Mother());
    simb::MCTrajectory const & mctraj = mcp.Trajectory();
    fgenie_particle_X.push_back(mctraj.Position(0).X());
    fgenie_particle_Y.push_back(mctraj.Position(0).Y());
    fgenie_particle_Z.push_back(mctraj.Position(0).Z());
    fgenie_particle_T.push_back(mctraj.Position(0).T());
    fgenie_particle_Px.push_back(mctraj.Momentum(0).Px());
    fgenie_particle_Py.push_back(mctraj.Momentum(0).Py());
    fgenie_particle_Pz.push_back(mctraj.Momentum(0).Pz());
    fgenie_particle_E.push_back(mctraj.Momentum(0).Pz());    
  }
  
}
  
  
void FillLightEvent::FillMCTrackVectors(art::Event const & e) {

  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctr = e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");  
  size_t const size = ev_mctr->size();

  fmctrack_Origin.reserve(size);
  fmctrack_PdgCode.reserve(size);
  fmctrack_TrackID.reserve(size);
  fmctrack_Process.reserve(size);
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
  fmctrack_contributed_charge.resize(size, std::vector<double>());  

  for(size_t i = 0; i < size; ++i) {

    sim::MCTrack const & mctr = ev_mctr->at(i);
    size_t const mctr_size = mctr.size();

    fmctrack_Origin.push_back(mctr.Origin());
    fmctrack_PdgCode.push_back(mctr.PdgCode());
    fmctrack_TrackID.push_back(mctr.TrackID());
    fmctrack_Process.push_back(mctr.Process());

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

    std::vector<double> & mctrack_contributed_charge = fmctrack_contributed_charge.at(i);
    mctrack_contributed_charge.reserve(ftrack_producers.size());

    for(size_t j = 0; j < frmcm_size; ++j) {
      RecoMCMatching const & rmcm = frmcm.at(j);
      mctrack_contributed_charge.push_back(rmcm.GetMCTrackCharge().at(i));
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
  fmcshower_contributed_charge.resize(size, std::vector<double>());

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
    fmcshower_DaughterTrackID.push_back(mcs.DaughterTrackID());
    fmcshower_Charge.push_back(mcs.Charge());
    fmcshower_dQdx.push_back(mcs.Charge());
    fmcshower_StartDir_X.push_back(mcs.StartDir().X());
    fmcshower_StartDir_Y.push_back(mcs.StartDir().Y());
    fmcshower_StartDir_Z.push_back(mcs.StartDir().Z());

    std::vector<double> & mcshower_contributed_charge = fmcshower_contributed_charge.at(i);
    mcshower_contributed_charge.reserve(fshower_producers.size());

    for(size_t j = 0; j < frmcm_size; ++j) {
      RecoMCMatching const & rmcm = frmcm.at(j);
      mcshower_contributed_charge.push_back(rmcm.GetMCShowerCharge().at(i));
    }

  }

}


void FillLightEvent::FillTruth(art::Event const & e,
			       size_t & delta_rad_mct_index) {
  
  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");  
  
  if(ffs.Run(e, delta_rad_mct_index)) fis_delta_rad = 1;
  else fis_delta_rad = 0;
  if(fis_delta_rad == 1) {
    if(delta_rad_mct_index == SIZE_MAX) {
      std::cout << "delta_rad_mct_index == SIZE_MAX\n";
      return;
    }
    if(delta_rad_mct_index >= ev_mct->size()) {
      std::cout << "delta_rad_mct_index: " << delta_rad_mct_index
		<< " >= ev_mct->size(): " << ev_mct->size() << std::endl;
      return;
    }
  }
  if(fis_delta_rad == 0) delta_rad_mct_index = 0;

  simb::MCTruth const & mct = ev_mct->at(delta_rad_mct_index);

  FillGenieParticleVectors(mct);
  FillMCTrackVectors(e);
  FillMCShowerVectors(e);

  simb::MCNeutrino const & mcn = mct.GetNeutrino();

  fnu_pdg = mcn.Nu().PdgCode();
  fnu_energy = mcn.Nu().Trajectory().E(0);
  flep_pdg = mcn.Lepton().PdgCode();
  flep_energy = mcn.Lepton().Trajectory().E(0);
  fccnc = mcn.CCNC();
  fmode = mcn.Mode();
  finteraction_type = mcn.InteractionType();

  TLorentzVector const & true_nu_pos = mcn.Nu().Trajectory().Position(0);
  if(ftpc_volume.Contain(geoalgo::Point_t(true_nu_pos))) ftrue_nu_vtx_tpc_contained = 1;
  else ftrue_nu_vtx_tpc_contained = 0;
  if(ffiducial_volume.Contain(geoalgo::Point_t(true_nu_pos))) ftrue_nu_vtx_fid_contained = 1;
  else ftrue_nu_vtx_fid_contained = 0;

  ftrue_nuvertx = true_nu_pos.X();
  ftrue_nuverty = true_nu_pos.Y();
  ftrue_nuvertz = true_nu_pos.Z();

  ftrue_nu_E = mcn.Nu().Trajectory().E(0);

}


void FillLightEvent::GetDeltaMCShowerMCTrackIndices(art::Event const & e,
						    size_t const delta_rad_mct_index,
						    size_t & delta_photon_index,
						    size_t & delta_mcshower_index,
						    size_t & delta_proton_index,
						    size_t & delta_mctrack_index) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mctruth =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");  
  art::ValidHandle<std::vector<sim::MCTrack>> const & ev_mctrack =
    e.getValidHandle<std::vector<sim::MCTrack>>("mcreco");
  art::ValidHandle<std::vector<sim::MCShower>> const & ev_mcshower =
    e.getValidHandle<std::vector<sim::MCShower>>("mcreco");

  simb::MCTruth const & mct = ev_mctruth->at(delta_rad_mct_index);

  int delta_index = -1;
  for(int i = 0; i < mct.NParticles(); ++i) {
    simb::MCParticle const & mcp = mct.GetParticle(i);
    if(abs(mcp.PdgCode()) == 2214 || abs(mcp.PdgCode()) == 2114) {
      delta_index = i;
      fdelta_true_pdg = mcp.PdgCode();
      fdelta_true_energy = mcp.Trajectory().E(0);
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

  fdelta_photon_energy = mct.GetParticle(delta_photon_index).Trajectory().E(0);
  if(delta_proton_index != SIZE_MAX) 
    fdelta_proton_energy = mct.GetParticle(delta_proton_index).Trajectory().E(0);

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

  if(delta_mcshower_index != SIZE_MAX) {
    sim::MCShower const & mcs = ev_mcshower->at(delta_mcshower_index);
    fdelta_mcshower_true_pdg = mcs.PdgCode();
    fdelta_mcshower_true_energy = mcs.Start().E();
    fdelta_mcshower_detprofile_energy = mcs.DetProfile().E();
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

  if(delta_mctrack_index != SIZE_MAX) {
    sim::MCTrack const & mctr = ev_mctrack->at(delta_mctrack_index);
    fdelta_mctrack_true_pdg = mctr.PdgCode();
    fdelta_mctrack_true_energy = mctr.Start().E();
    fdelta_mctrack_true_length = geoalgo::Point_t(mctr.Start().Position()).Dist(mctr.End().Position());
  }

}


void FillLightEvent::analyze(art::Event const & e) {

  /*

  art::ValidHandle<std::vector<recob::Shower>> const & ev_s = e.getValidHandle<std::vector<recob::Shower>>(fshower_producer);
  */

  ResetEvent();

  frun_number = e.id().run();
  fsubrun_number = e.id().subRun();
  fevent_number = e.id().event();

  std::cout << "1\n";

  for(RecoMCMatching & rmcm : frmcm) rmcm.MatchWAssociations(e);

  std::cout << "2\n";

  FillRecoHitVectors(e);
  std::cout << "3\n";
  FillRecoTrackVectors(e);
  std::cout << "4\n";

  if(fmc) {

    art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct = e.getValidHandle<std::vector<simb::MCTruth>>("generator");  

    FillWeights(e);

    size_t delta_rad_mct_index = SIZE_MAX;

    if(ev_mct->front().GetNeutrino().Nu().Mother() == -1) FillTruth(e, delta_rad_mct_index);

    size_t delta_photon_index = SIZE_MAX;
    size_t delta_mcshower_index = SIZE_MAX;
    size_t delta_proton_index = SIZE_MAX;
    size_t delta_mctrack_index = SIZE_MAX;

    if(fis_delta_rad == 1) {
      GetDeltaMCShowerMCTrackIndices(e,
				     delta_rad_mct_index,
				     delta_photon_index,
				     delta_mcshower_index,
				     delta_proton_index,
				     delta_mctrack_index);
    }

    fdelta_photon_index = delta_photon_index;
    fdelta_mcshower_index = delta_mcshower_index;
    fdelta_proton_index = delta_proton_index;
    fdelta_mctrack_index = delta_mctrack_index;

  }  
  
  fevent_tree->Fill();

}


DEFINE_ART_MODULE(FillLightEvent)
