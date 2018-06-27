#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <vector>
#include <string>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "uboone/DQMTools/DQMFlashAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "TH1F.h"

namespace dqm {
  class DQMFlashModule;
}

class dqm::DQMFlashModule : public art::EDAnalyzer {
public:
  explicit DQMFlashModule(fhicl::ParameterSet const & p);
  virtual ~DQMFlashModule();

  void analyze(art::Event const & e) override;

  void beginJob() override;

  /* import from GoodRunSelectionAna_module.cc for optically pre-filtered data */
     // Allow for fhicl parameters to possibly change during processing...
     void reconfigure(fhicl::ParameterSet const&)  ;
     
     // Recover information from the start of a run (if processing across runs)
     void beginRun(const art::Run&) override;

     // Recover information at the end of a run (if processing across runs)
     void endRun(const art::Run&) override;

     // Called when job completes to deal with output of stuff from beginJob
     void endJob() override;

private:

  TH1F *NFlashes;

  TH1F *MeanFlashLight;
  TH1F *FlashLightVar;
  TH1F *FlashLightSkew;

  DQMFlashAlg analysisAlg;

  /* import from GoodRunSelectionAna_module.cc for optically pre-filtered data */
    TTree* fDataTree;
    
    int numevent = 0; //number of event counter

    int Run=-999;
    static const int n = 20000;
        
    int    Nentries[n]={0};
    int    Nevts[n]={0};
    float  Nrun[n]={0};
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
    
    float _frun;
    float _fnumevent; 
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

    TH1F* _run;      	
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
    
    std::string              fOpFlashModuleLabel;
    std::string              fSimpFlashBeamModuleLabel;
};


dqm::DQMFlashModule::DQMFlashModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ 
  this->reconfigure(p); 
}

   void dqm::DQMFlashModule::reconfigure(fhicl::ParameterSet const& pset)
   {

    fOpFlashModuleLabel         = pset.get< std::string > ("OpFlashModuleLabel","opflash");
    fSimpFlashBeamModuleLabel   = pset.get< std::string > ("SimpFlashBeamModuleLabel","simpleFlashBeam");

   }

dqm::DQMFlashModule::~DQMFlashModule()
{
  // Clean up dynamic memory and other resources here.
}

void dqm::DQMFlashModule::analyze(art::Event const & e)
{
  //get the flash data
  art::Handle< std::vector<recob::OpFlash> > flashHandle;
  std::string _flash_producer_name = "simpleFlashBeam";
  e.getByLabel(_flash_producer_name,flashHandle);
  if(!flashHandle.isValid()) {
    //std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::OpFlash from "
    //<< _flash_producer_name << std::endl;
    //throw std::exception();
    return;
  }
  std::vector<recob::OpFlash> const& flashVector(*flashHandle);
 
  analysisAlg.AnalyzeFlashes(flashVector, NFlashes, MeanFlashLight, FlashLightVar, FlashLightSkew);

  /* import from GoodRunSelectionAna_module.cc for optically pre-filtered data */

    _run ->Fill(e.run());
    
    numevent +=1;

   /////////////Track related quantities///////////////////


   ///////////Flash related quantities///////////////

    int opflash_all=0, opflash50=0, opflash20=0, opflash0_20=0;   
    art::Handle< std::vector<recob::OpFlash>> opflashlistCol;
    e.getByLabel(fOpFlashModuleLabel,opflashlistCol);
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
    e.getByLabel(fSimpFlashBeamModuleLabel,simpflashlistCol);
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
                
}

void dqm::DQMFlashModule::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  NFlashes       = tfs->make<TH1F>("NFlashes"      , "# of flashes"        , 300, 0, 300);
  MeanFlashLight = tfs->make<TH1F>("MeanFlashLight", "Mean light per flash", 100, 0, 500);
  FlashLightVar  = tfs->make<TH1F>("FlashLightVar" , "Variance on light per flash", 1E3, 0, 1E6);
  FlashLightSkew = tfs->make<TH1F>("FlashLightSkew", "Skewness on light per flash", 100, -500, 500);
 
  /* import from GoodRunSelectionAna_module.cc for optically pre-filtered data */
    fDataTree = tfs->make<TTree>("fDataTree","Data Holder");

    fDataTree->Branch("_frun",&_frun,"_frun/F");
    fDataTree->Branch("_fnumevent",&_fnumevent,"_fnumevent/F");

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
}


   void dqm::DQMFlashModule::beginRun(const art::Run& run)
  {
 _run                               = new TH1F("_run","run number",20000,0,20000);  
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
  }
   
  void dqm::DQMFlashModule::endRun(const art::Run& run)
  {
        Run                    =_run->GetMean(1);
	Nevts[Run]             += numevent;
	numevent = 0; //resetting value
	Nentries[Run]            +=1;
	Nrun[Run]                +=Run;	
        std::cout<<"\n"<<"\n"<<"run num: "<<Run<<" number of events: "<<Nevts[Run]<<"\n"<<"\n";
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

  
  }

   void dqm::DQMFlashModule::endJob()
  {
  for(int i=0;i<n;i++)
  {
    if(Nentries[i]!=0 && Nrun[i]!=0)
    {
        _frun =  Nrun[i]/Nentries[i];
	_fnumevent = Nevts[i];
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
	
	fDataTree->Fill();

	}
      }	

  }

DEFINE_ART_MODULE(dqm::DQMFlashModule)
