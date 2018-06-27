
/**
 * @file    GoodRunSelectionAna_module.cc
 * @brief   Module for monitoring good run selection variables
 * @Author: Aleena Rafique (aleena@ksu.edu)
 * 
 ******************************************************************************/

  #ifndef GOODRUNSELECTIONANA_H
  #define GOODRUNSELECTIONANA_H
  
  // Framework includes
#include "cetlib_except/exception.h"
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
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
  // LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
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
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"

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

  
  // C++ Includes
  #include <iostream>
  #include <cstring>
  #include <sys/stat.h>
  #include <vector>
  #include "TMatrixD.h"
  #include "TVector3.h"
  
  // root
  #include "TH1D.h"
  #include <TMath.h>
  #include "TVectorD.h"
  #include "TFile.h"
  #include "TGeoManager.h"
  #include "TF1.h"
  #include "TF2.h"
  #include "TGraph.h"
  #include "TMath.h"
  #include "TTree.h"
  #include "TH2F.h"
  #include "TProfile.h"
  #include <TStyle.h>
  #include <TCanvas.h>
  
  namespace GoodRunSelectionAna {  
   
    class GoodRunSelectionAna : public art::EDAnalyzer
    {  

    public:
 
     explicit GoodRunSelectionAna(fhicl::ParameterSet const &pset);
     virtual ~GoodRunSelectionAna();                        
     
      // The actual method for analyzing the events
     void analyze(art::Event const& evt) ;
     //void beginSubRun(const art::SubRun& sr);
     
     // Allow for fhicl parameters to possibly change during processing...
     void reconfigure(fhicl::ParameterSet const&)  ;
     
     // Called when job begins for definitions of histograms/tuples/etc.   
     void beginJob()  ;
     
     // Called when job completes to deal with output of stuff from beginJob
     void endJob();
    
     // Recover information from the start of a run (if processing across runs)
     void beginRun(const art::Run&);

     // Recover information at the end of a run (if processing across runs)
     void endRun(const art::Run&);     

    // Length of reconstructed track, trajectory by trajectory.
    double length(const recob::Track& track);

    TTree* fDataTree;
    
    int numevent = 0; //number of event counter

    int Run=-999;
    static const int n = 20000;
        
    int    Nentries[n]={0};
    int    Nevts[n]={0};
    float  Nrun[n]={0};
    float  Ntrack_pmtrack[n]={0};
    float  NtrackRMS_pmtrack[n]={0};
    float  Ntrack_panCos[n]={0};
    float  NtrackRMS_panCos[n]={0};
    float  Ntrack_panNu[n]={0};
    float  NtrackRMS_panNu[n]={0};
    float  Ntrack_panNuPMA[n]={0};
    float  NtrackRMS_panNuPMA[n]={0};
    float  Ntrklen_pmtrack[n]={0};
    float  NtrklenRMS_pmtrack[n]={0};
    float  Ntrklen_panCos[n]={0};
    float  NtrklenRMS_panCos[n]={0};
    float  Ntrklen_panNu[n]={0};
    float  NtrklenRMS_panNu[n]={0};
    float  Ntrklen_panNuPMA[n]={0};
    float  NtrklenRMS_panNuPMA[n]={0};
    float  Nhit[n]={0};
    float  NhitRMS[n]={0};
    float  NhitU[n]={0};
    float  NhitURMS[n]={0};
    float  NhitV[n]={0};
    float  NhitVRMS[n]={0};
    float  NhitY[n]={0};
    float  NhitYRMS[n]={0}; 
    float  Nhitph[n]={0};
    float  NhitphRMS[n]={0};
    float  NhitphU[n]={0};
    float  NhitphURMS[n]={0};
    float  NhitphV[n]={0};
    float  NhitphVRMS[n]={0};
    float  NhitphY[n]={0};
    float  NhitphYRMS[n]={0};
    float  Nhitcharge[n]={0};
    float  NhitchargeRMS[n]={0};
    float  NhitchargeU[n]={0};
    float  NhitchargeURMS[n]={0};
    float  NhitchargeV[n]={0};
    float  NhitchargeVRMS[n]={0};
    float  NhitchargeY[n]={0};
    float  NhitchargeYRMS[n]={0};
    float  Nflash_opFlashSat[n]={0};
    float  NflashRMS_opFlashSat[n]={0};
    float  Nflash50_opFlashSat[n]={0};
    float  Nflash50RMS_opFlashSat[n]={0};
    float  Nflash20_opFlashSat[n]={0};
    float  Nflash20RMS_opFlashSat[n]={0};
    float  Nflash0_20_opFlashSat[n]={0};
    float  Nflash0_20RMS_opFlashSat[n]={0};
    float  Nflashy_opFlashSat[n]={0};
    float  NflashyRMS_opFlashSat[n]={0};
    float  Nflashy50_opFlashSat[n]={0};
    float  Nflashy50RMS_opFlashSat[n]={0};
    float  Nflashy20_opFlashSat[n]={0};
    float  Nflashy20RMS_opFlashSat[n]={0};
    float  Nflashy0_20_opFlashSat[n]={0};
    float  Nflashy0_20RMS_opFlashSat[n]={0};
    float  Nflashz_opFlashSat[n]={0};
    float  NflashzRMS_opFlashSat[n]={0};
    float  Nflashz50_opFlashSat[n]={0};
    float  Nflashz50RMS_opFlashSat[n]={0};
    float  Nflashz20_opFlashSat[n]={0};
    float  Nflashz20RMS_opFlashSat[n]={0};
    float  Nflashz0_20_opFlashSat[n]={0};
    float  Nflashz0_20RMS_opFlashSat[n]={0};
    float  Nflashpe_opFlashSat[n]={0};
    float  NflashpeRMS_opFlashSat[n]={0};
    float  Nflash_simpleFlashBeam[n]={0};
    float  NflashRMS_simpleFlashBeam[n]={0};
    float  Nflash50_simpleFlashBeam[n]={0};
    float  Nflash50RMS_simpleFlashBeam[n]={0};
    float  Nflash20_simpleFlashBeam[n]={0};
    float  Nflash20RMS_simpleFlashBeam[n]={0};
    float  Nflash0_20_simpleFlashBeam[n]={0};
    float  Nflash0_20RMS_simpleFlashBeam[n]={0};
    float  Nflashy_simpleFlashBeam[n]={0};
    float  NflashyRMS_simpleFlashBeam[n]={0};
    float  Nflashy50_simpleFlashBeam[n]={0};
    float  Nflashy50RMS_simpleFlashBeam[n]={0};
    float  Nflashy20_simpleFlashBeam[n]={0};
    float  Nflashy20RMS_simpleFlashBeam[n]={0};
    float  Nflashy0_20_simpleFlashBeam[n]={0};
    float  Nflashy0_20RMS_simpleFlashBeam[n]={0};
    float  Nflashz_simpleFlashBeam[n]={0};
    float  NflashzRMS_simpleFlashBeam[n]={0};
    float  Nflashz50_simpleFlashBeam[n]={0};
    float  Nflashz50RMS_simpleFlashBeam[n]={0};
    float  Nflashz20_simpleFlashBeam[n]={0};
    float  Nflashz20RMS_simpleFlashBeam[n]={0};
    float  Nflashz0_20_simpleFlashBeam[n]={0};
    float  Nflashz0_20RMS_simpleFlashBeam[n]={0};
    float  Nflashpe_simpleFlashBeam[n]={0};
    float  NflashpeRMS_simpleFlashBeam[n]={0};
    float  Nvtx_pmtrk[n]={0};
    float  NvtxRMS_pmtrk[n]={0};
    float  Nvtx_panCos[n]={0};
    float  NvtxRMS_panCos[n]={0};
    float  Nvtx_panNu[n]={0};
    float  NvtxRMS_panNu[n]={0};
    float  Nvtx_panNuPMA[n]={0};
    float  NvtxRMS_panNuPMA[n]={0};
    
    float _frun;
    float _fnumevent; 
    float _f_mean_ntrack_pmtrack;
    float _f_mean_ntrack_pandoraNu;
    float _f_mean_ntrack_pandoraNuPMA;    
    float _f_mean_ntrack_pandoraCosmic;     
    float _f_mean_trklen_pmtrack;
    float _f_mean_trklen_pandoraNu;
    float _f_mean_trklen_pandoraNuPMA;    
    float _f_mean_trklen_pandoraCosmic;     
    float _f_mean_nhit_tot;
    float _f_mean_nhit_U;
    float _f_mean_nhit_V;     
    float _f_mean_nhit_Y; 
    float _f_mean_hitPH_tot;
    float _f_mean_hitPH_U;
    float _f_mean_hitPH_V;     
    float _f_mean_hitPH_Y; 
    float _f_mean_hitCharge_tot;
    float _f_mean_hitCharge_U;
    float _f_mean_hitCharge_V;     
    float _f_mean_hitCharge_Y;	  
    float _f_mean_nflashTot_opFlashSat;     
    float _f_mean_nflashPE50_opFlashSat;     
    float _f_mean_nflashPE20_opFlashSat;     
    float _f_mean_nflashPE0_20_opFlashSat;     
    float _f_mean_nflashTot_simpleFlashBeam;     
    float _f_mean_nflashPE50_simpleFlashBeam;     
    float _f_mean_nflashPE20_simpleFlashBeam;     
    float _f_mean_nflashPE0_20_simpleFlashBeam;     
    float _f_mean_flashycenter_opFlashSat;     
    float _f_mean_flashycenterPE50_opFlashSat;     
    float _f_mean_flashycenterPE20_opFlashSat;     
    float _f_mean_flashycenterPE0_20_opFlashSat;     
    float _f_mean_flashzcenter_opFlashSat;     
    float _f_mean_flashzcenterPE50_opFlashSat;     
    float _f_mean_flashzcenterPE20_opFlashSat;     
    float _f_mean_flashzcenterPE0_20_opFlashSat;     
    float _f_mean_flashycenter_simpleFlashBeam;     
    float _f_mean_flashycenterPE50_simpleFlashBeam;     
    float _f_mean_flashycenterPE20_simpleFlashBeam;	
    float _f_mean_flashycenterPE0_20_simpleFlashBeam;	
    float _f_mean_flashzcenter_simpleFlashBeam;     
    float _f_mean_flashzcenterPE50_simpleFlashBeam;     
    float _f_mean_flashzcenterPE20_simpleFlashBeam;	
    float _f_mean_flashzcenterPE0_20_simpleFlashBeam;	
    float _f_mean_flashPE_opFlashSat;     
    float _f_mean_flashPE_simpleFlashBeam;     
    float _f_mean_nvrtx_pmtrack;
    float _f_mean_nvrtx_pandoraNu;
    float _f_mean_nvrtx_pandoraNuPMA;
    float _f_mean_nvrtx_pandoraCosmic;
    
    float _f_rms_ntrack_pmtrack;
    float _f_rms_ntrack_pandoraNu;
    float _f_rms_ntrack_pandoraNuPMA;    
    float _f_rms_ntrack_pandoraCosmic;     
    float _f_rms_trklen_pmtrack;
    float _f_rms_trklen_pandoraNu;
    float _f_rms_trklen_pandoraNuPMA;    
    float _f_rms_trklen_pandoraCosmic;     
    float _f_rms_nhit_tot;
    float _f_rms_nhit_U;
    float _f_rms_nhit_V;     
    float _f_rms_nhit_Y; 
    float _f_rms_hitPH_tot;
    float _f_rms_hitPH_U;
    float _f_rms_hitPH_V;     
    float _f_rms_hitPH_Y; 
    float _f_rms_hitCharge_tot;
    float _f_rms_hitCharge_U;
    float _f_rms_hitCharge_V;     
    float _f_rms_hitCharge_Y;	  
    float _f_rms_nflashTot_opFlashSat;     
    float _f_rms_nflashPE50_opFlashSat;     
    float _f_rms_nflashPE20_opFlashSat;     
    float _f_rms_nflashPE0_20_opFlashSat;     
    float _f_rms_nflashTot_simpleFlashBeam;	
    float _f_rms_nflashPE50_simpleFlashBeam;     
    float _f_rms_nflashPE20_simpleFlashBeam;     
    float _f_rms_nflashPE0_20_simpleFlashBeam;     
    float _f_rms_flashycenter_opFlashSat;     
    float _f_rms_flashycenterPE50_opFlashSat;     
    float _f_rms_flashycenterPE20_opFlashSat;     
    float _f_rms_flashycenterPE0_20_opFlashSat;     
    float _f_rms_flashzcenter_opFlashSat;     
    float _f_rms_flashzcenterPE50_opFlashSat;     
    float _f_rms_flashzcenterPE20_opFlashSat;     
    float _f_rms_flashzcenterPE0_20_opFlashSat;     
    float _f_rms_flashycenter_simpleFlashBeam;     
    float _f_rms_flashycenterPE50_simpleFlashBeam;     
    float _f_rms_flashycenterPE20_simpleFlashBeam;	
    float _f_rms_flashycenterPE0_20_simpleFlashBeam;	
    float _f_rms_flashzcenter_simpleFlashBeam;     
    float _f_rms_flashzcenterPE50_simpleFlashBeam;     
    float _f_rms_flashzcenterPE20_simpleFlashBeam;	
    float _f_rms_flashzcenterPE0_20_simpleFlashBeam;	
    float _f_rms_flashPE_opFlashSat;     
    float _f_rms_flashPE_simpleFlashBeam;     
    float _f_rms_nvrtx_pmtrack;
    float _f_rms_nvrtx_pandoraNu;
    float _f_rms_nvrtx_pandoraNuPMA;
    float _f_rms_nvrtx_pandoraCosmic;	     
     
    TH1F* _run;      	
    TH1F* _ntrack_pmtrack ;
    TH1F* _ntrack_pandoraNu    ;
    TH1F* _ntrack_pandoraNuPMA;
    TH1F* _ntrack_pandoraCosmic ;
    TH1F* _trklen_pmtrack  ;
    TH1F* _trklen_pandoraNu     ;
    TH1F* _trklen_pandoraNuPMA  ;
    TH1F* _trklen_pandoraCosmic ;
    TH1F* _nhit                 ;
    TH1F* _nhitU                ;
    TH1F* _nhitV                ;
    TH1F* _nhitY                ;
    TH1F* _hit_ph               ;
    TH1F* _hit_phU              ;
    TH1F* _hit_phV              ;
    TH1F* _hit_phY              ;
    TH1F* _hit_charge           ;
    TH1F* _hit_chargeU          ;
    TH1F* _hit_chargeV          ;
    TH1F* _hit_chargeY          ;
    TH1F* _nflash_opFlashSat    ;
    TH1F* _nflash50_opFlashSat  ;
    TH1F* _nflash20_opFlashSat  ;
    TH1F* _nflash0_20_opFlashSat;
    TH1F* _nflash_simpleFlashBeam ;
    TH1F* _nflash50_simpleFlashBeam  ;
    TH1F* _nflash20_simpleFlashBeam  ;
    TH1F* _nflash0_20_simpleFlashBeam;
    TH1F* _flash_ycenter_opFlashSat  ;
    TH1F* _flash_ycenter50_opFlashSat;
    TH1F* _flash_ycenter20_opFlashSat;
    TH1F* _flash_ycenter0_20_opFlashSat;
    TH1F* _flash_zcenter_opFlashSat  ;
    TH1F* _flash_zcenter50_opFlashSat;
    TH1F* _flash_zcenter20_opFlashSat;
    TH1F* _flash_zcenter0_20_opFlashSat;
    TH1F* _flash_ycenter_simpleFlashBeam   ;
    TH1F* _flash_ycenter50_simpleFlashBeam ;
    TH1F* _flash_ycenter20_simpleFlashBeam ;
    TH1F* _flash_ycenter0_20_simpleFlashBeam ;
    TH1F* _flash_zcenter_simpleFlashBeam   ;
    TH1F* _flash_zcenter50_simpleFlashBeam ;
    TH1F* _flash_zcenter20_simpleFlashBeam ;
    TH1F* _flash_zcenter0_20_simpleFlashBeam ;
    TH1F* _flash_pe_opFlashSat     ;      
    TH1F* _flash_pe_simpleFlashBeam;      
    TH1F* _nvrtx_pmtrack;
    TH1F* _nvrtx_pandoraNu;
    TH1F* _nvrtx_pandoraNuPMA;
    TH1F* _nvrtx_pandoraCosmic;

      private:
  
    std::string              fHitsModuleLabel;
    std::string              fpmtrackModuleLabel;
    std::string              fPanCosTrackModuleLabel;
    std::string              fPanNuTrackModuleLabel;
    std::string              fPanNuPMATrackModuleLabel;
    std::string              fpmtrackVrtxModuleLabel;
    std::string              fPanCosVrtxModuleLabel;
    std::string              fPanNuVrtxModuleLabel;
    std::string              fPanNuPMAVrtxModuleLabel;
    std::string              fOpFlashModuleLabel;
    std::string              fSimpFlashBeamModuleLabel;

   };
  
   //-----------------------------------------------------------------------
    // Constructor
   GoodRunSelectionAna::GoodRunSelectionAna(fhicl::ParameterSet const& pset)
   : EDAnalyzer(pset)
  {
      this->reconfigure(pset);
  }
  
    //-----------------------------------------------------------------------
    // Destructor
    GoodRunSelectionAna::~GoodRunSelectionAna() 
   {
    }
 
   //-----------------------------------------------------------------------
   void GoodRunSelectionAna::beginJob()
   {

    art::ServiceHandle<art::TFileService> tfs;

    fDataTree = tfs->make<TTree>("fDataTree","Data Holder");

    fDataTree->Branch("_frun",&_frun,"_frun/F");
    fDataTree->Branch("_fnumevent",&_fnumevent,"_fnumevent/F"); 
   
    fDataTree->Branch("_f_mean_ntrack_pmtrack",&_f_mean_ntrack_pmtrack,"_f_mean_ntrack_pmtrack/F");   
    fDataTree->Branch("_f_mean_ntrack_pandoraNu",&_f_mean_ntrack_pandoraNu,"_f_mean_ntrack_pandoraNu/F");   
    fDataTree->Branch("_f_mean_ntrack_pandoraNuPMA",&_f_mean_ntrack_pandoraNuPMA,"_f_mean_ntrack_pandoraNuPMA/F");   
    fDataTree->Branch("_f_mean_ntrack_pandoraCosmic",&_f_mean_ntrack_pandoraCosmic,"_f_mean_ntrack_pandoraCosmic/F");     
    fDataTree->Branch("_f_mean_trklen_pmtrack",&_f_mean_trklen_pmtrack,"_f_mean_trklen_pmtrack/F");
    fDataTree->Branch("_f_mean_trklen_pandoraNu",&_f_mean_trklen_pandoraNu,"_f_mean_trklen_pandoraNu/F");
    fDataTree->Branch("_f_mean_trklen_pandoraNuPMA",&_f_mean_trklen_pandoraNuPMA,"_f_mean_trklen_pandoraNuPMA/F");
    fDataTree->Branch("_f_mean_trklen_pandoraCosmic",&_f_mean_trklen_pandoraCosmic,"_f_mean_trklen_pandoraCosmic/F");     
    fDataTree->Branch("_f_mean_nhit_tot",&_f_mean_nhit_tot,"_f_mean_nhit_tot/F");
    fDataTree->Branch("_f_mean_nhit_U",&_f_mean_nhit_U,"_f_mean_nhit_U/F");
    fDataTree->Branch("_f_mean_nhit_V",&_f_mean_nhit_V,"_f_mean_nhit_V/F");     
    fDataTree->Branch("_f_mean_nhit_Y",&_f_mean_nhit_Y,"_f_mean_nhit_Y/F"); 
    fDataTree->Branch("_f_mean_hitPH_tot",&_f_mean_hitPH_tot,"_f_mean_hitPH_tot/F");
    fDataTree->Branch("_f_mean_hitPH_U",&_f_mean_hitPH_U,"_f_mean_hitPH_U/F");
    fDataTree->Branch("_f_mean_hitPH_V",&_f_mean_hitPH_V,"_f_mean_hitPH_V/F");     
    fDataTree->Branch("_f_mean_hitPH_Y",&_f_mean_hitPH_Y,"_f_mean_hitPH_Y/F"); 
    fDataTree->Branch("_f_mean_hitCharge_tot",&_f_mean_hitCharge_tot,"_f_mean_hitCharge_tot/F");
    fDataTree->Branch("_f_mean_hitCharge_U",&_f_mean_hitCharge_U,"_f_mean_hitCharge_U/F");
    fDataTree->Branch("_f_mean_hitCharge_V",&_f_mean_hitCharge_V,"_f_mean_hitCharge_V/F");     
    fDataTree->Branch("_f_mean_hitCharge_Y",&_f_mean_hitCharge_Y,"_f_mean_hitCharge_Y/F");	  
    fDataTree->Branch("_f_mean_nflashTot_opFlashSat",&_f_mean_nflashTot_opFlashSat,"_f_mean_nflashTot_opFlashSat/F");     
    fDataTree->Branch("_f_mean_nflashPE50_opFlashSat",&_f_mean_nflashPE50_opFlashSat,"_f_mean_nflashPE50_opFlashSat/F");     
    fDataTree->Branch("_f_mean_nflashPE20_opFlashSat",&_f_mean_nflashPE20_opFlashSat,"_f_mean_nflashPE20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_nflashPE0_20_opFlashSat",&_f_mean_nflashPE0_20_opFlashSat,"_f_mean_nflashPE0_20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_nflashTot_simpleFlashBeam",&_f_mean_nflashTot_simpleFlashBeam,"_f_mean_nflashTot_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_nflashPE50_simpleFlashBeam",&_f_mean_nflashPE50_simpleFlashBeam,"_f_mean_nflashPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_nflashPE20_simpleFlashBeam",&_f_mean_nflashPE20_simpleFlashBeam,"_f_mean_nflashPE20_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_nflashPE0_20_simpleFlashBeam",&_f_mean_nflashPE0_20_simpleFlashBeam,"_f_mean_nflashPE0_20_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_flashycenter_opFlashSat",&_f_mean_flashycenter_opFlashSat,"_f_mean_flashycenter_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashycenterPE50_opFlashSat",&_f_mean_flashycenterPE50_opFlashSat,"_f_mean_flashycenterPE50_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashycenterPE20_opFlashSat",&_f_mean_flashycenterPE20_opFlashSat,"_f_mean_flashycenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashycenterPE0_20_opFlashSat",&_f_mean_flashycenterPE20_opFlashSat,"_f_mean_flashycenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashzcenter_opFlashSat",&_f_mean_flashzcenter_opFlashSat,"_f_mean_flashzcenter_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashzcenterPE50_opFlashSat",&_f_mean_flashzcenterPE50_opFlashSat,"_f_mean_flashzcenterPE50_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashzcenterPE20_opFlashSat",&_f_mean_flashzcenterPE20_opFlashSat,"_f_mean_flashzcenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashzcenterPE0_20_opFlashSat",&_f_mean_flashzcenterPE20_opFlashSat,"_f_mean_flashzcenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashycenter_simpleFlashBeam",&_f_mean_flashycenter_simpleFlashBeam,"_f_mean_flashycenter_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_flashycenterPE50_simpleFlashBeam",&_f_mean_flashycenterPE50_simpleFlashBeam,"_f_mean_flashycenterPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_flashycenterPE20_simpleFlashBeam",&_f_mean_flashycenterPE20_simpleFlashBeam,"_f_mean_flashycenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_mean_flashycenterPE0_20_simpleFlashBeam",&_f_mean_flashycenterPE20_simpleFlashBeam,"_f_mean_flashycenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_mean_flashzcenter_simpleFlashBeam",&_f_mean_flashzcenter_simpleFlashBeam,"_f_mean_flashzcenter_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_flashzcenterPE50_simpleFlashBeam",&_f_mean_flashzcenterPE50_simpleFlashBeam,"_f_mean_flashzcenterPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_flashzcenterPE20_simpleFlashBeam",&_f_mean_flashzcenterPE20_simpleFlashBeam,"_f_mean_flashzcenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_mean_flashzcenterPE0_20_simpleFlashBeam",&_f_mean_flashzcenterPE20_simpleFlashBeam,"_f_mean_flashzcenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_mean_flashPE_opFlashSat",&_f_mean_flashPE_opFlashSat,"_f_mean_flashPE_opFlashSat/F");     
    fDataTree->Branch("_f_mean_flashPE_simpleFlashBeam",&_f_mean_flashPE_simpleFlashBeam,"_f_mean_flashPE_simpleFlashBeam/F");     
    fDataTree->Branch("_f_mean_nvrtx_pmtrack",&_f_mean_nvrtx_pmtrack,"_f_mean_nvrtx_pmtrack/F");
    fDataTree->Branch("_f_mean_nvrtx_pandoraNu",&_f_mean_nvrtx_pandoraNu,"_f_mean_nvrtx_pandoraNu/F");
    fDataTree->Branch("_f_mean_nvrtx_pandoraNuPMA",&_f_mean_nvrtx_pandoraNuPMA,"_f_mean_nvrtx_pandoraNuPMA/F");
    fDataTree->Branch("_f_mean_nvrtx_pandoraCosmic",&_f_mean_nvrtx_pandoraCosmic,"_f_mean_nvrtx_pandoraCosmic/F");

    fDataTree->Branch("_f_rms_ntrack_pmtrack",&_f_rms_ntrack_pmtrack,"_f_rms_ntrack_pmtrack/F");
    fDataTree->Branch("_f_rms_ntrack_pandoraNu",&_f_rms_ntrack_pandoraNu,"_f_rms_ntrack_pandoraNu/F");
    fDataTree->Branch("_f_rms_ntrack_pandoraNuPMA",&_f_rms_ntrack_pandoraNuPMA,"_f_rms_ntrack_pandoraNuPMA/F");
    fDataTree->Branch("_f_rms_ntrack_pandoraCosmic",&_f_rms_ntrack_pandoraCosmic,"_f_rms_ntrack_pandoraCosmic/F");     
    fDataTree->Branch("_f_rms_trklen_pmtrack",&_f_rms_trklen_pmtrack,"_f_rms_trklen_pmtrack/F");
    fDataTree->Branch("_f_rms_trklen_pandoraNu",&_f_rms_trklen_pandoraNu,"_f_rms_trklen_pandoraNu/F");
    fDataTree->Branch("_f_rms_trklen_pandoraNuPMA",&_f_rms_trklen_pandoraNuPMA,"_f_rms_trklen_pandoraNuPMA/F");
    fDataTree->Branch("_f_rms_trklen_pandoraCosmic",&_f_rms_trklen_pandoraCosmic,"_f_rms_trklen_pandoraCosmic/F");     
    fDataTree->Branch("_f_rms_nhit_tot",&_f_rms_nhit_tot,"_f_rms_nhit_tot/F");
    fDataTree->Branch("_f_rms_nhit_U",&_f_rms_nhit_U,"_f_rms_nhit_U/F");
    fDataTree->Branch("_f_rms_nhit_V",&_f_rms_nhit_V,"_f_rms_nhit_V/F");     
    fDataTree->Branch("_f_rms_nhit_Y",&_f_rms_nhit_Y,"_f_rms_nhit_Y/F"); 
    fDataTree->Branch("_f_rms_hitPH_tot",&_f_rms_hitPH_tot,"_f_rms_hitPH_tot/F");
    fDataTree->Branch("_f_rms_hitPH_U",&_f_rms_hitPH_U,"_f_rms_hitPH_U/F");
    fDataTree->Branch("_f_rms_hitPH_V",&_f_rms_hitPH_V,"_f_rms_hitPH_V/F");     
    fDataTree->Branch("_f_rms_hitPH_Y",&_f_rms_hitPH_Y,"_f_rms_hitPH_Y/F"); 
    fDataTree->Branch("_f_rms_hitCharge_tot",&_f_rms_hitCharge_tot,"_f_rms_hitCharge_tot/F");
    fDataTree->Branch("_f_rms_hitCharge_U",&_f_rms_hitCharge_U,"_f_rms_hitCharge_U/F");
    fDataTree->Branch("_f_rms_hitCharge_V",&_f_rms_hitCharge_V,"_f_rms_hitCharge_V/F");     
    fDataTree->Branch("_f_rms_hitCharge_Y",&_f_rms_hitCharge_Y,"_f_rms_hitCharge_Y/F");	  
    fDataTree->Branch("_f_rms_nflashTot_opFlashSat",&_f_rms_nflashTot_opFlashSat,"_f_rms_nflashTot_opFlashSat/F");     
    fDataTree->Branch("_f_rms_nflashPE50_opFlashSat",&_f_rms_nflashPE50_opFlashSat,"_f_rms_nflashPE50_opFlashSat/F");     
    fDataTree->Branch("_f_rms_nflashPE20_opFlashSat",&_f_rms_nflashPE20_opFlashSat,"_f_rms_nflashPE20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_nflashPE0_20_opFlashSat",&_f_rms_nflashPE0_20_opFlashSat,"_f_rms_nflashPE0_20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_nflashTot_simpleFlashBeam",&_f_rms_nflashTot_simpleFlashBeam,"_f_rms_nflashTot_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_nflashPE50_simpleFlashBeam",&_f_rms_nflashPE50_simpleFlashBeam,"_f_rms_nflashPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_nflashPE20_simpleFlashBeam",&_f_rms_nflashPE20_simpleFlashBeam,"_f_rms_nflashPE20_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_nflashPE0_20_simpleFlashBeam",&_f_rms_nflashPE0_20_simpleFlashBeam,"_f_rms_nflashPE0_20_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_flashycenter_opFlashSat",&_f_rms_flashycenter_opFlashSat,"_f_rms_flashycenter_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashycenterPE50_opFlashSat",&_f_rms_flashycenterPE50_opFlashSat,"_f_rms_flashycenterPE50_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashycenterPE20_opFlashSat",&_f_rms_flashycenterPE20_opFlashSat,"_f_rms_flashycenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashycenterPE0_20_opFlashSat",&_f_rms_flashycenterPE20_opFlashSat,"_f_rms_flashycenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashzcenter_opFlashSat",&_f_rms_flashzcenter_opFlashSat,"_f_rms_flashzcenter_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashzcenterPE50_opFlashSat",&_f_rms_flashzcenterPE50_opFlashSat,"_f_rms_flashzcenterPE50_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashzcenterPE20_opFlashSat",&_f_rms_flashzcenterPE20_opFlashSat,"_f_rms_flashzcenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashzcenterPE0_20_opFlashSat",&_f_rms_flashzcenterPE20_opFlashSat,"_f_rms_flashzcenterPE20_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashycenter_simpleFlashBeam",&_f_rms_flashycenter_simpleFlashBeam,"_f_rms_flashycenter_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_flashycenterPE50_simpleFlashBeam",&_f_rms_flashycenterPE50_simpleFlashBeam,"_f_rms_flashycenterPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_flashycenterPE20_simpleFlashBeam",&_f_rms_flashycenterPE20_simpleFlashBeam,"_f_rms_flashycenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_rms_flashycenterPE0_20_simpleFlashBeam",&_f_rms_flashycenterPE20_simpleFlashBeam,"_f_rms_flashycenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_rms_flashzcenter_simpleFlashBeam",&_f_rms_flashzcenter_simpleFlashBeam,"_f_rms_flashzcenter_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_flashzcenterPE50_simpleFlashBeam",&_f_rms_flashzcenterPE50_simpleFlashBeam,"_f_rms_flashzcenterPE50_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_flashzcenterPE20_simpleFlashBeam",&_f_rms_flashzcenterPE20_simpleFlashBeam,"_f_rms_flashzcenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_rms_flashzcenterPE0_20_simpleFlashBeam",&_f_rms_flashzcenterPE20_simpleFlashBeam,"_f_rms_flashzcenterPE20_simpleFlashBeam/F");	
    fDataTree->Branch("_f_rms_flashPE_opFlashSat",&_f_rms_flashPE_opFlashSat,"_f_rms_flashPE_opFlashSat/F");     
    fDataTree->Branch("_f_rms_flashPE_simpleFlashBeam",&_f_rms_flashPE_simpleFlashBeam,"_f_rms_flashPE_simpleFlashBeam/F");     
    fDataTree->Branch("_f_rms_nvrtx_pmtrack",&_f_rms_nvrtx_pmtrack,"_f_rms_nvrtx_pmtrack/F");
    fDataTree->Branch("_f_rms_nvrtx_pandoraNu",&_f_rms_nvrtx_pandoraNu,"_f_rms_nvrtx_pandoraNu/F");
    fDataTree->Branch("_f_rms_nvrtx_pandoraNuPMA",&_f_rms_nvrtx_pandoraNuPMA,"_f_rms_nvrtx_pandoraNuPMA/F");
    fDataTree->Branch("_f_rms_nvrtx_pandoraCosmic",&_f_rms_nvrtx_pandoraCosmic,"_f_rms_nvrtx_pandoraCosmic/F");
                            
   }
   
   //-----------------------------------------------------------------------
   void GoodRunSelectionAna::endJob()
  {
  for(int i=0;i<n;i++)
  {
    if(Nentries[i]!=0 && Nrun[i]!=0)
    {
        _frun =  Nrun[i]/Nentries[i];
	_fnumevent = Nevts[i];
  	_f_mean_ntrack_pmtrack             = Ntrack_pmtrack[i]/Nentries[i];
  	_f_rms_ntrack_pmtrack              = NtrackRMS_pmtrack[i]/Nentries[i];
	_f_mean_ntrack_pandoraCosmic       = Ntrack_panCos[i]/Nentries[i];
	_f_rms_ntrack_pandoraCosmic        = NtrackRMS_panCos[i]/Nentries[i];
	_f_mean_ntrack_pandoraNu           = Ntrack_panNu[i]/Nentries[i];
	_f_rms_ntrack_pandoraNu            = NtrackRMS_panNu[i]/Nentries[i];
	_f_mean_ntrack_pandoraNuPMA        = Ntrack_panNuPMA[i]/Nentries[i];
	_f_rms_ntrack_pandoraNuPMA         = NtrackRMS_panNuPMA[i]/Nentries[i];	
	_f_mean_trklen_pmtrack	           = Ntrklen_pmtrack[i]/Nentries[i];
	_f_rms_trklen_pmtrack	           = NtrklenRMS_pmtrack[i]/Nentries[i];
	_f_mean_trklen_pandoraCosmic       = Ntrklen_panCos[i]/Nentries[i];
	_f_rms_trklen_pandoraCosmic        = NtrklenRMS_panCos[i]/Nentries[i];
	_f_mean_trklen_pandoraNu           = Ntrklen_panNu[i]/Nentries[i];
	_f_rms_trklen_pandoraNu            = NtrklenRMS_panNu[i]/Nentries[i];
	_f_mean_trklen_pandoraNuPMA        = Ntrklen_panNuPMA[i]/Nentries[i];
	_f_rms_trklen_pandoraNuPMA         = NtrklenRMS_panNuPMA[i]/Nentries[i];
	_f_mean_nhit_tot                   = Nhit[i]/Nentries[i];
	_f_rms_nhit_tot                    = NhitRMS[i]/Nentries[i];
	_f_mean_nhit_U                     = NhitU[i]/Nentries[i];
	_f_rms_nhit_U                      = NhitURMS[i]/Nentries[i];
	_f_mean_nhit_V                     = NhitV[i]/Nentries[i];
	_f_rms_nhit_V                      = NhitVRMS[i]/Nentries[i];
	_f_mean_nhit_Y                     = NhitY[i]/Nentries[i];
	_f_rms_nhit_Y                      = NhitYRMS[i]/Nentries[i];
	_f_mean_hitPH_tot                  = Nhitph[i]/Nentries[i];
	_f_rms_hitPH_tot                   = NhitphRMS[i]/Nentries[i];
	_f_mean_hitPH_U                    = NhitphU[i]/Nentries[i];
	_f_rms_hitPH_U                     = NhitphURMS[i]/Nentries[i];
	_f_mean_hitPH_V                    = NhitphV[i]/Nentries[i];
	_f_rms_hitPH_V                     = NhitphVRMS[i]/Nentries[i];
	_f_mean_hitPH_Y                    = NhitphY[i]/Nentries[i];
	_f_rms_hitPH_Y                     = NhitphYRMS[i]/Nentries[i];      
	_f_mean_hitCharge_tot              = Nhitcharge[i]/Nentries[i];
	_f_rms_hitCharge_tot               = NhitchargeRMS[i]/Nentries[i];
	_f_mean_hitCharge_U                = NhitchargeU[i]/Nentries[i];
	_f_rms_hitCharge_U                 = NhitchargeURMS[i]/Nentries[i];
	_f_mean_hitCharge_V                = NhitchargeV[i]/Nentries[i];
	_f_rms_hitCharge_V                 = NhitchargeVRMS[i]/Nentries[i];
	_f_mean_hitCharge_Y                = NhitchargeY[i]/Nentries[i];
	_f_rms_hitCharge_Y                 = NhitchargeYRMS[i]/Nentries[i];
	_f_mean_nflashTot_opFlashSat       = Nflash_opFlashSat[i]/Nentries[i];
	_f_rms_nflashTot_opFlashSat        = NflashRMS_opFlashSat[i]/Nentries[i];
	_f_mean_nflashPE50_opFlashSat      = Nflash50_opFlashSat[i]/Nentries[i];
	_f_rms_nflashPE50_opFlashSat       = Nflash50RMS_opFlashSat[i]/Nentries[i];
	_f_mean_nflashPE20_opFlashSat      = Nflash20_opFlashSat[i]/Nentries[i];
	_f_rms_nflashPE20_opFlashSat       = Nflash20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_nflashPE0_20_opFlashSat    = Nflash0_20_opFlashSat[i]/Nentries[i];
	_f_rms_nflashPE0_20_opFlashSat     = Nflash0_20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashycenter_opFlashSat               = Nflashy_opFlashSat[i]/Nentries[i];
	_f_rms_flashycenter_opFlashSat                = NflashyRMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashycenterPE50_opFlashSat           = Nflashy50_opFlashSat[i]/Nentries[i];
	_f_rms_flashycenterPE50_opFlashSat            = Nflashy50RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashycenterPE20_opFlashSat           = Nflashy20_opFlashSat[i]/Nentries[i];
	_f_rms_flashycenterPE20_opFlashSat            = Nflashy20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashycenterPE0_20_opFlashSat         = Nflashy0_20_opFlashSat[i]/Nentries[i];
	_f_rms_flashycenterPE0_20_opFlashSat          = Nflashy0_20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashzcenter_opFlashSat               = Nflashz_opFlashSat[i]/Nentries[i];
	_f_rms_flashzcenter_opFlashSat                = NflashzRMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashzcenterPE50_opFlashSat           = Nflashz50_opFlashSat[i]/Nentries[i];
	_f_rms_flashzcenterPE50_opFlashSat           = Nflashz50RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashzcenterPE20_opFlashSat           = Nflashz20_opFlashSat[i]/Nentries[i];
	_f_rms_flashzcenterPE20_opFlashSat           = Nflashz20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashzcenterPE0_20_opFlashSat         = Nflashz0_20_opFlashSat[i]/Nentries[i];
	_f_rms_flashzcenterPE0_20_opFlashSat         = Nflashz0_20RMS_opFlashSat[i]/Nentries[i];
	_f_mean_flashPE_opFlashSat                    = Nflashpe_opFlashSat[i]/Nentries[i];
	_f_rms_flashPE_opFlashSat                     = NflashpeRMS_opFlashSat[i]/Nentries[i];
	_f_mean_nflashTot_simpleFlashBeam             = Nflash_simpleFlashBeam[i]/Nentries[i];
	_f_rms_nflashTot_simpleFlashBeam              = NflashRMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_nflashPE50_simpleFlashBeam            = Nflash50_simpleFlashBeam[i]/Nentries[i];
	_f_rms_nflashPE50_simpleFlashBeam             = Nflash50RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_nflashPE20_simpleFlashBeam            = Nflash20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_nflashPE20_simpleFlashBeam             = Nflash20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_nflashPE0_20_simpleFlashBeam          = Nflash0_20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_nflashPE0_20_simpleFlashBeam           = Nflash0_20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashycenter_simpleFlashBeam          = Nflashy_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashycenter_simpleFlashBeam           = NflashyRMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashycenterPE50_simpleFlashBeam      = Nflashy50_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashycenterPE50_simpleFlashBeam       = Nflashy50RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashycenterPE20_simpleFlashBeam      = Nflashy20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashycenterPE20_simpleFlashBeam       = Nflashy20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashycenterPE0_20_simpleFlashBeam    = Nflashy0_20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashycenterPE0_20_simpleFlashBeam     = Nflashy0_20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashzcenter_simpleFlashBeam          = Nflashz_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashzcenter_simpleFlashBeam           = NflashzRMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashzcenterPE50_simpleFlashBeam      = Nflashz50_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashzcenterPE50_simpleFlashBeam      = Nflashz50RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashzcenterPE20_simpleFlashBeam      = Nflashz20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashzcenterPE20_simpleFlashBeam      = Nflashz20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashzcenterPE0_20_simpleFlashBeam    = Nflashz0_20_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashzcenterPE0_20_simpleFlashBeam    = Nflashz0_20RMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_flashPE_simpleFlashBeam               = Nflashpe_simpleFlashBeam[i]/Nentries[i];
	_f_rms_flashPE_simpleFlashBeam               = NflashpeRMS_simpleFlashBeam[i]/Nentries[i];
	_f_mean_nvrtx_pmtrack                         = Nvtx_pmtrk[i]/Nentries[i];
	_f_rms_nvrtx_pmtrack                          = NvtxRMS_pmtrk[i]/Nentries[i];
	_f_mean_nvrtx_pandoraCosmic                   = Nvtx_panCos[i]/Nentries[i];
	_f_rms_nvrtx_pandoraCosmic                    = NvtxRMS_panCos[i]/Nentries[i];
	_f_mean_nvrtx_pandoraNu                       = Nvtx_panNu[i]/Nentries[i];
	_f_rms_nvrtx_pandoraNu                        = NvtxRMS_panNu[i]/Nentries[i];
	_f_mean_nvrtx_pandoraNuPMA                    = Nvtx_panNuPMA[i]/Nentries[i];
	_f_rms_nvrtx_pandoraNuPMA                     = NvtxRMS_panNuPMA[i]/Nentries[i];
	
	fDataTree->Fill();

	}
      }	

  }
  
   //-----------------------------------------------------------------------
   void GoodRunSelectionAna::beginRun(const art::Run& run)
  {
 _run                               = new TH1F("_run","run number",20000,0,20000);  
 _ntrack_pmtrack                    = new TH1F("_ntrack_pmtrack","Number of tracks (pmtrack) per event",100,0,100);
 _ntrack_pandoraNu                  = new TH1F("_ntrack_pandoraNu","Number of tracks (pandoraNu) per event",100,0,100);
 _ntrack_pandoraNuPMA               = new TH1F("_ntrack_pandoraNuPMA","Number of tracks (pandoraNu) per event",100,0,100);
 _ntrack_pandoraCosmic              = new TH1F("_ntrack_pandoraCosmic","Number of tracks (pandoraCosmic) per event",100,0,100);
 _trklen_pmtrack                    = new TH1F("_trklen_pmtrack","track length (pmtrack) per event",100,0,700000);
 _trklen_pandoraNu                  = new TH1F("_trklen_pandoraNu","track length (pandoraNu) per event",100,0,700000);
 _trklen_pandoraNuPMA               = new TH1F("_trklen_pandoraNuPMA","track length (pandoraNu) per event",100,0,700000);
 _trklen_pandoraCosmic              = new TH1F("_trklen_pandoraCosmic","track length (pandoraCosmic) per event",100,0,700000);
 _nhit                              = new TH1F("_nhit","Number of hits per event ",100,0,50000);
 _nhitU                             = new TH1F("_nhitU","Number of U-plane hits per event ",100,0,50000);
 _nhitV                             = new TH1F("_nhitV","Number of V-plane hits per event ",100,0,50000);
 _nhitY                             = new TH1F("_nhitY","Number of Y-plane hits per event ",100,0,50000);
 _hit_ph                            = new TH1F("_hit_ph","Number of hit ph per event ",100,0,500);
 _hit_phU                           = new TH1F("_hit_phU","Number of U-plane hit ph per event ",100,0,500);
 _hit_phV                           = new TH1F("_hit_phV","Number of V-plane hit ph per event ",100,0,500);
 _hit_phY                           = new TH1F("_hit_phY","Number of Y-plane hit ph per event ",100,0,500);
 _hit_charge                        = new TH1F("_hit_charge","Number of hit charge per event ",100,0,1600000);
 _hit_chargeU                       = new TH1F("_hit_chargeU","Number of U-plane hit charge per event ",100,0,1600000);
 _hit_chargeV                       = new TH1F("_hit_chargeV","Number of V-plane hit charge per event ",100,0,1600000);
 _hit_chargeY                       = new TH1F("_hit_chargeY","Number of Y-plane hit charge per event ",100,0,1600000);
 _nflash_opFlashSat                 = new TH1F("_flash_opFlashSat","Number of optical sat flashes per event ",100,0,20000);
 _nflash50_opFlashSat               = new TH1F("_nflash50_opFlashSat","Number of > 50 optical sat flashes per event ",100,0,20000);
 _nflash20_opFlashSat               = new TH1F("_nflash20_opFlashSat","Number of > 20 optical sat flashes per event ",100,0,20000);
 _nflash0_20_opFlashSat             = new TH1F("_nflash0_20_opFlashSat","Number of < 20 optical sat flashes per event ",100,0,20000);
 _nflash_simpleFlashBeam            = new TH1F("_nflash_simpleFlashBeam","Number of simpleFlashBeam flashes per event ",100,0,20000);
 _nflash50_simpleFlashBeam          = new TH1F("_nflash50_simpleFlashBeam","Number of > 50 simpleFlashBeam flashes per event ",100,0,20000);
 _nflash20_simpleFlashBeam          = new TH1F("_nflash20_simpleFlashBeam","Number of > 20 simpleFlashBeam flashes per event ",100,0,20000);
 _nflash0_20_simpleFlashBeam        = new TH1F("_nflash0_20_simpleFlashBeam","Number of < 20 simpleFlashBeam flashes per event ",100,0,20000);
 _flash_ycenter_opFlashSat          = new TH1F("_flash_ycenter_opFlashSat","Number of y center of flashes per event ",100,-200,200);
 _flash_ycenter50_opFlashSat        = new TH1F("_flash_ycenter50_opFlashSat","y center of >50 flashes per event ",100,-200,200);
 _flash_ycenter20_opFlashSat        = new TH1F("_flash_ycenter20_opFlashSat","y center of >20 flashes per event ",100,-200,200);
 _flash_ycenter0_20_opFlashSat      = new TH1F("_flash_ycenter0_20_opFlashSat","y center of <20 flashes per event ",100,-200,200);
 _flash_zcenter_opFlashSat          = new TH1F("_flash_zcenter_opFlashSat","Number of z center of flashes per event",100,0,1500);
 _flash_zcenter50_opFlashSat        = new TH1F("_flash_zcenter50_opFlashSat","z center of >50 flashes per event",100,0,1500);
 _flash_zcenter20_opFlashSat        = new TH1F("_flash_zcenter20_opFlashSat","z center of >20 flashes per event",100,0,1500);
 _flash_zcenter0_20_opFlashSat      = new TH1F("_flash_zcenter0_20_opFlashSat","z center of <20 flashes per event",100,0,1500);
 _flash_ycenter_simpleFlashBeam     = new TH1F("_flash_ycenter_simpleFlashBeam","Number of y center of flashes per event ",100,-200,200);
 _flash_ycenter50_simpleFlashBeam   = new TH1F("_flash_ycenter50_simpleFlashBeam","y center of >50 flashes per event ",100,-200,200);
 _flash_ycenter20_simpleFlashBeam   = new TH1F("_flash_ycenter20_simpleFlashBeam","y center of >20 flashes per event ",100,-200,200);
 _flash_ycenter0_20_simpleFlashBeam = new TH1F("_flash_ycenter0_20_simpleFlashBeam","y center of <20 flashes per event ",100,-200,200);
 _flash_zcenter_simpleFlashBeam     = new TH1F("_flash_zcenter_simpleFlashBeam","Number of z center of flashes per event",100,0,1500);
 _flash_zcenter50_simpleFlashBeam   = new TH1F("_flash_zcenter50_simpleFlashBeam","z center of >50 flashes per event",100,0,1500);
 _flash_zcenter20_simpleFlashBeam   = new TH1F("_flash_zcenter20_simpleFlashBeam","z center of >20 flashes per event",100,0,1500);
 _flash_zcenter0_20_simpleFlashBeam = new TH1F("_flash_zcenter0_20_simpleFlashBeam","z center of <20 flashes per event",100,0,1500);
 _flash_pe_opFlashSat               = new TH1F("_flash_pe_opFlashSat","Number of flashes pe per event ",100,0,7000);      
 _flash_pe_simpleFlashBeam          = new TH1F("_flash_pe_simpleFlashBeam","Number of flashes pe per event ",100,0,7000);      
 _nvrtx_pmtrack                     = new TH1F("_nvrtx_pmtrack","number of vertex (pmtrack) per event",100,0,600);
 _nvrtx_pandoraNu                   = new TH1F("_nvrtx_pandoraNu","number of vertex (pandoraNu) per event",100,0,600);
 _nvrtx_pandoraNuPMA                = new TH1F("_nvrtx_pandoraNuPMA","number of vertex (pandoraNu) per event",100,0,600);
 _nvrtx_pandoraCosmic               = new TH1F("_nvrtx_pandoraCosmic","number of vertex (pandoraCosmic) per event",100,0,600);
  }
  
   //-----------------------------------------------------------------------
   void GoodRunSelectionAna::endRun(const art::Run& run)
  {
        Run                    =_run->GetMean(1);
	Nevts[Run]             += numevent;
	numevent = 0; //resetting value
	Nentries[Run]            +=1;
	Nrun[Run]                +=Run;	
        std::cout<<"\n"<<"\n"<<"run num: "<<Run<<" number of events: "<<Nevts[Run]<<"\n"<<"\n";
	Ntrack_pmtrack[Run]          +=_ntrack_pmtrack->GetMean(1);
	NtrackRMS_pmtrack[Run]       +=_ntrack_pmtrack->GetRMS(1);
	Ntrack_panNu[Run]       +=_ntrack_pandoraNu->GetMean(1);
	NtrackRMS_panNu[Run]    +=_ntrack_pandoraNu->GetRMS(1);
	Ntrack_panNuPMA[Run]        +=_ntrack_pandoraNuPMA->GetMean(1);
	NtrackRMS_panNuPMA[Run]     +=_ntrack_pandoraNuPMA->GetRMS(1);
	Ntrack_panCos[Run]    +=_ntrack_pandoraCosmic->GetMean(1);
	NtrackRMS_panCos[Run]    +=_ntrack_pandoraCosmic->GetRMS(1);
	Ntrklen_pmtrack[Run]         +=_trklen_pmtrack->GetMean(1);
	NtrklenRMS_pmtrack[Run]      +=_trklen_pmtrack->GetRMS(1);
	Ntrklen_panNu[Run]       +=_trklen_pandoraNu->GetMean(1);
	NtrklenRMS_panNu[Run]    +=_trklen_pandoraNu->GetRMS(1);
	Ntrklen_panNuPMA[Run]       +=_trklen_pandoraNuPMA->GetMean(1);
	NtrklenRMS_panNuPMA[Run]    +=_trklen_pandoraNuPMA->GetRMS(1);
	Ntrklen_panCos[Run]      +=_trklen_pandoraCosmic->GetMean(1);
	NtrklenRMS_panCos[Run]   +=_trklen_pandoraCosmic->GetRMS(1);
	Nhit[Run]                +=_nhit->GetMean(1);
	NhitRMS[Run]             +=_nhit->GetRMS(1);
	NhitU[Run]               +=_nhitU->GetMean(1);
	NhitURMS[Run]            +=_nhitU->GetRMS(1);
	NhitV[Run]               +=_nhitV->GetMean(1);
	NhitVRMS[Run]               +=_nhitV->GetRMS(1);
	NhitY[Run]               +=_nhitY->GetMean(1);
	NhitYRMS[Run]               +=_nhitY->GetRMS(1);
	Nhitph[Run]              +=_hit_ph->GetMean(1);
	NhitphRMS[Run]              +=_hit_ph->GetRMS(1);
	NhitphU[Run]              +=_hit_phU->GetMean(1);
	NhitphURMS[Run]              +=_hit_phU->GetRMS(1);
	NhitphV[Run]             +=_hit_phV->GetMean(1);
	NhitphVRMS[Run]             +=_hit_phV->GetRMS(1);
	NhitphY[Run]             +=_hit_phY->GetMean(1);
	NhitphYRMS[Run]             +=_hit_phY->GetRMS(1);
	Nhitcharge[Run]          +=_hit_charge->GetMean(1);
	NhitchargeRMS[Run]          +=_hit_charge->GetRMS(1);
	NhitchargeU[Run]          +=_hit_chargeU->GetMean(1);
	NhitchargeURMS[Run]          +=_hit_chargeU->GetRMS(1);
	NhitchargeV[Run]          +=_hit_chargeV->GetMean(1);
	NhitchargeVRMS[Run]          +=_hit_chargeV->GetRMS(1);
	NhitchargeY[Run]          +=_hit_chargeY->GetMean(1);
	NhitchargeYRMS[Run]          +=_hit_chargeY->GetRMS(1);
	Nflash_opFlashSat[Run]              +=_nflash_opFlashSat->GetMean(1);
	NflashRMS_opFlashSat[Run]              +=_nflash_opFlashSat->GetRMS(1);
	Nflash50_opFlashSat[Run]            +=_nflash50_opFlashSat->GetMean(1);
	Nflash50RMS_opFlashSat[Run]            +=_nflash50_opFlashSat->GetRMS(1);
	Nflash20_opFlashSat[Run]            +=_nflash20_opFlashSat->GetMean(1);
	Nflash20RMS_opFlashSat[Run]            +=_nflash20_opFlashSat->GetRMS(1);
	Nflash0_20_opFlashSat[Run]            +=_nflash0_20_opFlashSat->GetMean(1);
	Nflash0_20RMS_opFlashSat[Run]            +=_nflash0_20_opFlashSat->GetRMS(1);
	Nflash_simpleFlashBeam[Run]              +=_nflash_simpleFlashBeam->GetMean(1);
	NflashRMS_simpleFlashBeam[Run]              +=_nflash_simpleFlashBeam->GetRMS(1);
	Nflash50_simpleFlashBeam[Run]            +=_nflash50_simpleFlashBeam->GetMean(1);
	Nflash50RMS_simpleFlashBeam[Run]            +=_nflash50_simpleFlashBeam->GetRMS(1);
	Nflash20_simpleFlashBeam[Run]            +=_nflash20_simpleFlashBeam->GetMean(1);
	Nflash20RMS_simpleFlashBeam[Run]            +=_nflash20_simpleFlashBeam->GetRMS(1);
	Nflash0_20_simpleFlashBeam[Run]            +=_nflash0_20_simpleFlashBeam->GetMean(1);
	Nflash0_20RMS_simpleFlashBeam[Run]            +=_nflash0_20_simpleFlashBeam->GetRMS(1);		
	Nflashy_opFlashSat[Run]             +=_flash_ycenter_opFlashSat->GetMean(1);
	NflashyRMS_opFlashSat[Run]             +=_flash_ycenter_opFlashSat->GetRMS(1);
	Nflashy50_opFlashSat[Run]           +=_flash_ycenter50_opFlashSat->GetMean(1);
	Nflashy50RMS_opFlashSat[Run]           +=_flash_ycenter50_opFlashSat->GetRMS(1);
	Nflashy20_opFlashSat[Run]           +=_flash_ycenter20_opFlashSat->GetMean(1);
	Nflashy20RMS_opFlashSat[Run]           +=_flash_ycenter20_opFlashSat->GetRMS(1);
	Nflashy0_20_opFlashSat[Run]           +=_flash_ycenter0_20_opFlashSat->GetMean(1);
	Nflashy0_20RMS_opFlashSat[Run]           +=_flash_ycenter0_20_opFlashSat->GetRMS(1);
	Nflashz_opFlashSat[Run]             +=_flash_zcenter_opFlashSat->GetMean(1);
	NflashzRMS_opFlashSat[Run]             +=_flash_zcenter_opFlashSat->GetRMS(1);
	Nflashz50_opFlashSat[Run]           +=_flash_zcenter50_opFlashSat->GetMean(1);
	Nflashz50RMS_opFlashSat[Run]           +=_flash_zcenter50_opFlashSat->GetRMS(1);
	Nflashz20_opFlashSat[Run]           +=_flash_zcenter20_opFlashSat->GetMean(1);
	Nflashz20RMS_opFlashSat[Run]           +=_flash_zcenter20_opFlashSat->GetRMS(1);
	Nflashz0_20_opFlashSat[Run]           +=_flash_zcenter0_20_opFlashSat->GetMean(1);
	Nflashz0_20RMS_opFlashSat[Run]           +=_flash_zcenter0_20_opFlashSat->GetRMS(1);
	Nflashy_simpleFlashBeam[Run]             +=_flash_ycenter_simpleFlashBeam->GetMean(1);
	NflashyRMS_simpleFlashBeam[Run]             +=_flash_ycenter_simpleFlashBeam->GetRMS(1);
	Nflashy50_simpleFlashBeam[Run]           +=_flash_ycenter50_simpleFlashBeam->GetMean(1);
	Nflashy50RMS_simpleFlashBeam[Run]           +=_flash_ycenter50_simpleFlashBeam->GetRMS(1);
	Nflashy20_simpleFlashBeam[Run]           +=_flash_ycenter20_simpleFlashBeam->GetMean(1);
	Nflashy20RMS_simpleFlashBeam[Run]           +=_flash_ycenter20_simpleFlashBeam->GetRMS(1);
	Nflashy0_20_simpleFlashBeam[Run]           +=_flash_ycenter0_20_simpleFlashBeam->GetMean(1);
	Nflashy0_20RMS_simpleFlashBeam[Run]           +=_flash_ycenter0_20_simpleFlashBeam->GetRMS(1);
	Nflashz_simpleFlashBeam[Run]             +=_flash_zcenter_simpleFlashBeam->GetMean(1);
	NflashzRMS_simpleFlashBeam[Run]             +=_flash_zcenter_simpleFlashBeam->GetRMS(1);
	Nflashz50_simpleFlashBeam[Run]           +=_flash_zcenter50_simpleFlashBeam->GetMean(1);
	Nflashz50RMS_simpleFlashBeam[Run]           +=_flash_zcenter50_simpleFlashBeam->GetRMS(1);
	Nflashz20_simpleFlashBeam[Run]           +=_flash_zcenter20_simpleFlashBeam->GetMean(1);
	Nflashz20RMS_simpleFlashBeam[Run]           +=_flash_zcenter20_simpleFlashBeam->GetRMS(1);
	Nflashz0_20_simpleFlashBeam[Run]           +=_flash_zcenter0_20_simpleFlashBeam->GetMean(1);
	Nflashz0_20RMS_simpleFlashBeam[Run]           +=_flash_zcenter0_20_simpleFlashBeam->GetRMS(1);
	Nflashpe_opFlashSat[Run]            +=_flash_pe_opFlashSat->GetMean(1);
	NflashpeRMS_opFlashSat[Run]            +=_flash_pe_opFlashSat->GetRMS(1);		
	Nflashpe_simpleFlashBeam[Run]            +=_flash_pe_simpleFlashBeam->GetMean(1);
	NflashpeRMS_simpleFlashBeam[Run]            +=_flash_pe_simpleFlashBeam->GetRMS(1);		
	Nvtx_pmtrk[Run]          +=_nvrtx_pmtrack->GetMean(1);
	NvtxRMS_pmtrk[Run]       +=_nvrtx_pmtrack->GetRMS(1);
	Nvtx_panNu[Run]         +=_nvrtx_pandoraNu->GetMean(1);
	NvtxRMS_panNu[Run]         +=_nvrtx_pandoraNu->GetRMS(1);	
	Nvtx_panNuPMA[Run]       +=_nvrtx_pandoraNuPMA->GetMean(1);
	NvtxRMS_panNuPMA[Run]       +=_nvrtx_pandoraNuPMA->GetRMS(1);	
	Nvtx_panCos[Run]         +=_nvrtx_pandoraCosmic->GetMean(1);
	NvtxRMS_panCos[Run]         +=_nvrtx_pandoraCosmic->GetRMS(1);

  
  }  
    //-----------------------------------------------------------------------
   void GoodRunSelectionAna::reconfigure(fhicl::ParameterSet const& pset)
   {

    fHitsModuleLabel            = pset.get< std::string > ("HitsModuleLabel","gaushit");
    fpmtrackModuleLabel         = pset.get< std::string > ("pmtrackModuleLabel","pmtrack");
    fPanCosTrackModuleLabel     = pset.get< std::string > ("PanCosTrackModuleLabel","pandoraCosmic");
    fPanNuTrackModuleLabel      = pset.get< std::string > ("PanNuTrackModuleLabel","pandoraNu");
    fPanNuPMATrackModuleLabel   = pset.get< std::string > ("PanNuPMATrackModuleLabel","pandoraNuPMA");
    fpmtrackVrtxModuleLabel     = pset.get< std::string > ("pmtrackVrtxModuleLabel","pmtrack");
    fPanNuVrtxModuleLabel       = pset.get< std::string > ("PanNuVrtxModuleLabel","pandoraNu");
    fPanNuPMAVrtxModuleLabel    = pset.get< std::string > ("PanNuPMAVrtxModuleLabel","pandoraNuPMA");
    fPanCosVrtxModuleLabel      = pset.get< std::string > ("PanCosVrtxModuleLabel","pandoraCosmic");
    fOpFlashModuleLabel         = pset.get< std::string > ("OpFlashModuleLabel","opflash");
    fSimpFlashBeamModuleLabel   = pset.get< std::string > ("SimpFlashBeamModuleLabel","simpleFlashBeam");

   }
 //========================================================================	
   
  // Length of reconstructed track, trajectory by trajectory.
  double GoodRunSelectionAna::length(const recob::Track& track)
  {
    double result = 0.;
    TVector3 disp = track.LocationAtPoint(0);
    int n = track.NumberTrajectoryPoints();

    for(int i = 1; i < n; ++i) 
    {
      const TVector3& pos = track.LocationAtPoint(i);
      disp -= pos;
      result += disp.Mag();
      disp = pos;
    }
    return result;
  } 
 //========================================================================	

   void GoodRunSelectionAna::analyze(art::Event const& evt) 
   {
   
    _run ->Fill(evt.run());
    
    numevent +=1;

   /////////////Track related quantities///////////////////

    art::Handle<std::vector<recob::Track> > pmtrkCol;	
    evt.getByLabel(fpmtrackModuleLabel, pmtrkCol);
    std::vector<recob::Track> const& pmtrkVector(*pmtrkCol);
    _ntrack_pmtrack->Fill(pmtrkVector.size());
    for (size_t it=0; it < pmtrkVector.size(); it++)
    {
      auto pmtrk = pmtrkVector.at(it);
      _trklen_pmtrack->Fill(length(pmtrk));
    } 
    
    art::Handle<std::vector<recob::Track> > panNuTrkCol;	
    evt.getByLabel(fPanNuTrackModuleLabel, panNuTrkCol);
    std::vector<recob::Track> const& panNuTrkVector(*panNuTrkCol);
    _ntrack_pandoraNu->Fill(panNuTrkVector.size());  
    for (size_t it=0; it < panNuTrkVector.size(); it++)
    {
      auto panNutrk = panNuTrkVector.at(it);
      _trklen_pandoraNu->Fill(length(panNutrk));
    }  
        
    art::Handle<std::vector<recob::Track> > panNuPMATrkCol;	
    evt.getByLabel(fPanNuPMATrackModuleLabel, panNuPMATrkCol);
    std::vector<recob::Track> const& panNuPMATrkVector(*panNuPMATrkCol);
    _ntrack_pandoraNuPMA->Fill(panNuPMATrkVector.size());    
    for (size_t it=0; it < panNuPMATrkVector.size(); it++)
    {
      auto panNuPMAtrk = panNuPMATrkVector.at(it);
      _trklen_pandoraNuPMA->Fill(length(panNuPMAtrk));
    }     
    art::Handle<std::vector<recob::Track> > panCosTrkCol;	
    evt.getByLabel(fPanCosTrackModuleLabel, panCosTrkCol);
    std::vector<recob::Track> const& panCosTrkVector(*panCosTrkCol);
    _ntrack_pandoraCosmic->Fill(panCosTrkVector.size());   
    for (size_t it=0; it < panCosTrkVector.size(); it++)
    {
      auto panCostrk = panCosTrkVector.at(it);
      _trklen_pandoraCosmic->Fill(length(panCostrk));
    }     
   
   ////////////Hit related quantities////////////////
    int hit_all=0, hitU=0, hitV=0, hitY=0;
    art::Handle< std::vector<recob::Hit>> hitListCol;
    evt.getByLabel(fHitsModuleLabel,hitListCol);
    std::vector<recob::Hit> const& hitlistVector(*hitListCol);
    _nhit->Fill(hitlistVector.size());

    for (size_t it=0; it < hitlistVector.size(); it++)
    {    
      auto gaushit = hitlistVector.at(it);
      hit_all+=1;
      _hit_ph->Fill(gaushit.PeakAmplitude());
      _hit_charge->Fill(gaushit.Integral());
      
      if(gaushit.WireID().Plane==0)
      {
        hitU+=1;
        _hit_phU->Fill(gaushit.PeakAmplitude());
        _hit_chargeU->Fill(gaushit.Integral());	
      }
      if(gaushit.WireID().Plane==1)
      {
        hitV+=1;
        _hit_phV->Fill(gaushit.PeakAmplitude());
        _hit_chargeV->Fill(gaushit.Integral());	
      }
      if(gaushit.WireID().Plane==2)
      {
        hitY+=1;
        _hit_phY->Fill(gaushit.PeakAmplitude());
        _hit_chargeY->Fill(gaushit.Integral());	
      }
    }
    _nhitU->Fill(hitU);
    _nhitV->Fill(hitV);
    _nhitY->Fill(hitY);

   ///////////Flash related quantities///////////////

    int opflash_all=0, opflash50=0, opflash20=0, opflash0_20=0;   
    art::Handle< std::vector<recob::OpFlash>> opflashlistCol;
    evt.getByLabel(fOpFlashModuleLabel,opflashlistCol);
    std::vector<recob::OpFlash> const& opflashlistVector(*opflashlistCol);
    _nflash_opFlashSat->Fill(opflashlistVector.size());
    
    for (size_t it=0; it < opflashlistVector.size(); it++)
    {    
      auto opflash = opflashlistVector.at(it);
      opflash_all+=1;
      _flash_ycenter_opFlashSat->Fill(opflash.YCenter()); 
      _flash_zcenter_opFlashSat->Fill(opflash.ZCenter());
      _flash_pe_opFlashSat->Fill(opflash.TotalPE()); 
      
      if(opflash.TotalPE() > 50)
      {
        opflash50+=1;
        _flash_ycenter50_opFlashSat->Fill(opflash.YCenter()); 
        _flash_zcenter50_opFlashSat->Fill(opflash.ZCenter());
      }
      if(opflash.TotalPE() > 20)
      {
        opflash20+=1;
        _flash_ycenter20_opFlashSat->Fill(opflash.YCenter()); 
        _flash_zcenter20_opFlashSat->Fill(opflash.ZCenter());
      }      
      if(opflash.TotalPE() < 20)
      {
        opflash0_20+=1;
        _flash_ycenter0_20_opFlashSat->Fill(opflash.YCenter()); 
        _flash_zcenter0_20_opFlashSat->Fill(opflash.ZCenter());
      }
    }
    _nflash50_opFlashSat->Fill(opflash50);
    _nflash20_opFlashSat->Fill(opflash20);
    _nflash0_20_opFlashSat->Fill(opflash0_20);
        
    
    int simpflash_all=0, simpflash50=0, simpflash20=0, simpflash0_20=0;   
    art::Handle< std::vector<recob::OpFlash>> simpflashlistCol;
    evt.getByLabel(fSimpFlashBeamModuleLabel,simpflashlistCol);
    std::vector<recob::OpFlash> const& simpflashlistVector(*simpflashlistCol);
    _nflash_simpleFlashBeam->Fill(simpflashlistVector.size());
    
    for (size_t it=0; it < simpflashlistVector.size(); it++)
    {    
      auto simpflash = simpflashlistVector.at(it);
      simpflash_all+=1;
      _flash_ycenter_simpleFlashBeam->Fill(simpflash.YCenter()); 
      _flash_zcenter_simpleFlashBeam->Fill(simpflash.ZCenter());
      _flash_pe_simpleFlashBeam->Fill(simpflash.TotalPE()); 

      if(simpflash.TotalPE() > 50)
      {
        simpflash50+=1;
        _flash_ycenter50_simpleFlashBeam->Fill(simpflash.YCenter()); 
        _flash_zcenter50_simpleFlashBeam->Fill(simpflash.ZCenter());
      }
      if(simpflash.TotalPE() > 20)
      {
        simpflash20+=1;
        _flash_ycenter20_simpleFlashBeam->Fill(simpflash.YCenter()); 
        _flash_zcenter20_simpleFlashBeam->Fill(simpflash.ZCenter());
      }      
      if(simpflash.TotalPE() < 20)
      {
        simpflash0_20+=1;
        _flash_ycenter0_20_simpleFlashBeam->Fill(simpflash.YCenter()); 
        _flash_zcenter0_20_simpleFlashBeam->Fill(simpflash.ZCenter());
      }
    }
    _nflash50_simpleFlashBeam->Fill(simpflash50);
    _nflash20_simpleFlashBeam->Fill(simpflash20);
    _nflash0_20_simpleFlashBeam->Fill(simpflash0_20);
                
   ///////////Vertex related quantities/////////////

    art::Handle<std::vector<recob::Vertex> > pmtrackVrtxCol;	
    evt.getByLabel(fpmtrackVrtxModuleLabel, pmtrackVrtxCol);
    std::vector<recob::Vertex> const& pmtrackVrtxVector(*pmtrackVrtxCol);
    _nvrtx_pmtrack->Fill(pmtrackVrtxVector.size());  

    art::Handle<std::vector<recob::Vertex> > panNuVrtxCol;	
    evt.getByLabel(fPanNuVrtxModuleLabel, panNuVrtxCol);
    std::vector<recob::Vertex> const& panNuVrtxVector(*panNuVrtxCol);
    _nvrtx_pandoraNu->Fill(panNuVrtxVector.size());  
    
    art::Handle<std::vector<recob::Vertex> > panNuPMAVrtxCol;	
    evt.getByLabel(fPanNuPMAVrtxModuleLabel, panNuPMAVrtxCol);
    std::vector<recob::Vertex> const& panNuPMAVrtxVector(*panNuPMAVrtxCol);
    _nvrtx_pandoraNuPMA->Fill(panNuPMAVrtxVector.size());  
    
    art::Handle<std::vector<recob::Vertex> > panCosVrtxCol;	
    evt.getByLabel(fPanCosVrtxModuleLabel, panCosVrtxCol);
    std::vector<recob::Vertex> const& panCosVrtxVector(*panCosVrtxCol);
    _nvrtx_pandoraCosmic->Fill(panCosVrtxVector.size());      
    
                			
 return;
 } // end GoodRunSelectionAna::analyze()
   	
   DEFINE_ART_MODULE(GoodRunSelectionAna)
 
} // namespace GoodRunSelectionAna_module
 
  #endif //GOODRUNSELECTIONANA_H
  
