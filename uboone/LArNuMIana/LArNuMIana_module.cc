// LArNuMIana_module.cc
// D. Davis <douglasdavis@utexas.edu>

#ifndef LArNuMIana_Module
#define LArNuMIana_Module

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCFlux.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "MCCheater/BackTracker.h"
#include "SummaryData/POTSummary.h"
#include "Simulation/ParticleList.h"
#include "Simulation/sim.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace LArNuMIana {

  class LArNuMIana : public art::EDAnalyzer {
    
  public:
    
    explicit LArNuMIana(fhicl::ParameterSet const& pset);
    virtual ~LArNuMIana();
    
    void beginJob();
    void beginSubRun(const art::SubRun& r);
    void analyze (const art::Event& evt);
    
  private:
    
    std::string fSimulationProducerLabel;
    std::string fGenieGenModuleLabel;    
    std::string fHitProducerLabel;       
    std::string fClusterProducerLabel;   
    std::string fPOTModuleLabel;

    TH1D  *fPOT;

    TTree *fNuMIEventNtuple;
    int    fEvent;
    int    fRun;
    int    fSubRun;
    bool   fCCint;
    bool   fCCQEint;
    bool   fNCint;
    bool   fNCQEint;
    int    fNuPdgCode;
    double fNuVx;
    double fNuVy;
    double fNuVz;
    double fNuE;
    double fNuPx;
    double fNuPy;
    double fNuPz;
    int    fLeptonPdgCode;
    double fLeptonVx;
    double fLeptonVy;
    double fLeptonVz;
    double fLeptonPx;
    double fLeptonPy;
    double fLeptonPz;
    double fLeptonE;
    double fLeptonThetaYZ;
    double fLeptonThetaXZ;
    double fLeptonThetaYZ2;
    double fLeptonThetaXZ2;

    int    flux_run;
    int    flux_evtno;
    double flux_tpx;
    double flux_tpy;
    double flux_tpz;
    int    flux_tptype;
    double flux_vx;
    double flux_vy;
    double flux_vz;
    int    flux_ndecay;
    int    flux_ppmedium;

    std::vector<std::string> fProcess;
    std::vector<std::string> fProdVolume;
    std::vector<std::string> fEndVolume;
    std::vector<std::string> fProdMaterial;
    std::vector<std::string> fEndMaterial;
    std::vector<int>         fParentID;
    std::vector<int>         fPdgCode;
    std::vector<int>         fTrackID;
    std::vector<double>      fStartVx;
    std::vector<double>      fStartVy;
    std::vector<double>      fStartVz;
    std::vector<double>      fStartPx;
    std::vector<double>      fStartPy;
    std::vector<double>      fStartPz;
    std::vector<double>      fStartE;
    std::vector<bool>        fInTPC;

    art::ServiceHandle<geo::Geometry> fGeometry;      
    double                            fElectronsToGeV;

  };
  
  LArNuMIana::LArNuMIana(fhicl::ParameterSet const& p) : EDAnalyzer(p)
  {
    fGenieGenModuleLabel     = p.get< std::string >("GenieGenModuleLabel");
    fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
    fPOTModuleLabel          = p.get< std::string >("POTModuleLabel","generator");
  }
  
  LArNuMIana::~LArNuMIana() {}
  
  void LArNuMIana::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    fPOT             = tfs->make<TH1D>("POT",";POT;",10000,0,1e22);

    fNuMIEventNtuple = tfs->make<TTree>("LArNuMIanaSimulation","LArNuMIanaSimulation");

    fNuMIEventNtuple->Branch("NuPdgCode",     &fNuPdgCode,     "NuPdgCode/I");
    fNuMIEventNtuple->Branch("NuPx",          &fNuPx,          "NuPx/D");
    fNuMIEventNtuple->Branch("NuPy",          &fNuPy,          "NuPy/D");
    fNuMIEventNtuple->Branch("NuPz",          &fNuPz,          "NuPz/D");
    fNuMIEventNtuple->Branch("NuVx",          &fNuVx,          "NuVx/D");
    fNuMIEventNtuple->Branch("NuVy",          &fNuVy,          "NuVy/D");
    fNuMIEventNtuple->Branch("NuVz",          &fNuVz,          "NuVz/D");
    fNuMIEventNtuple->Branch("NuE",           &fNuE,           "NuE/D");
    fNuMIEventNtuple->Branch("LeptonPdgCode", &fLeptonPdgCode, "LeptonPdgCode/I");
    fNuMIEventNtuple->Branch("LeptonVx",      &fLeptonVx,      "LeptonVx/D");
    fNuMIEventNtuple->Branch("LeptonVy",      &fLeptonVy,      "LeptonVy/D");
    fNuMIEventNtuple->Branch("LeptonVz",      &fLeptonVz,      "LeptonVz/D");
    fNuMIEventNtuple->Branch("LeptonPx",      &fLeptonPx,      "LeptonPx/D");
    fNuMIEventNtuple->Branch("LeptonPy",      &fLeptonPy,      "LeptonPy/D");
    fNuMIEventNtuple->Branch("LeptonPz",      &fLeptonPz,      "LeptonPz/D");
    fNuMIEventNtuple->Branch("LeptonE",       &fLeptonE,       "LeptonE/D");
    fNuMIEventNtuple->Branch("LeptonThetaXZ2",&fLeptonThetaXZ2,"LeptonThetaXZ2/D");
    fNuMIEventNtuple->Branch("LeptonThetaYZ2",&fLeptonThetaYZ2,"LeptonThetaYZ2/D");
    fNuMIEventNtuple->Branch("LeptonThetaXZ", &fLeptonThetaXZ, "LeptonThetaXZ/D");
    fNuMIEventNtuple->Branch("LeptonThetaYZ", &fLeptonThetaYZ, "LeptonThetaYZ/D");
    fNuMIEventNtuple->Branch("CCint",         &fCCint,         "CCint/O");
    fNuMIEventNtuple->Branch("CCQEint",       &fCCQEint,       "CCQEint/O");
    fNuMIEventNtuple->Branch("NCint",         &fNCint,         "NCint/O");
    fNuMIEventNtuple->Branch("NCQEint",       &fNCQEint,       "NCQEint/O");
    fNuMIEventNtuple->Branch("Event",         &fEvent,         "Event/I");
    fNuMIEventNtuple->Branch("SubRun",        &fSubRun,        "SubRun/I");
    fNuMIEventNtuple->Branch("Run",           &fRun,           "Run/I");
    fNuMIEventNtuple->Branch("Process",       &fProcess);
    fNuMIEventNtuple->Branch("ProdMaterial",  &fProdMaterial);
    fNuMIEventNtuple->Branch("EndMaterial",   &fEndMaterial);
    fNuMIEventNtuple->Branch("ProdVolume",    &fProdVolume);
    fNuMIEventNtuple->Branch("EndVolume",     &fEndVolume);
    fNuMIEventNtuple->Branch("TrackID",       &fTrackID);
    fNuMIEventNtuple->Branch("ParentID",      &fParentID);
    fNuMIEventNtuple->Branch("PdgCode",       &fPdgCode);
    fNuMIEventNtuple->Branch("StartVx",       &fStartVx);
    fNuMIEventNtuple->Branch("StartVy",       &fStartVy);
    fNuMIEventNtuple->Branch("StartVz",       &fStartVz);
    fNuMIEventNtuple->Branch("StartPx",       &fStartPx);
    fNuMIEventNtuple->Branch("StartPy",       &fStartPy);
    fNuMIEventNtuple->Branch("StartPz",       &fStartPz);
    fNuMIEventNtuple->Branch("StartE",        &fStartE);
    fNuMIEventNtuple->Branch("InTPC",         &fInTPC);

    fNuMIEventNtuple->Branch("flux_run",     &flux_run,     "flux_run/I");
    fNuMIEventNtuple->Branch("flux_evtno",   &flux_evtno,   "flux_evtno/I");
    fNuMIEventNtuple->Branch("flux_tpx",     &flux_tpx,     "flux_tpx/D");
    fNuMIEventNtuple->Branch("flux_tpy",     &flux_tpy,     "flux_tpy/D");
    fNuMIEventNtuple->Branch("flux_tpz",     &flux_tpz,     "flux_tpz/D");
    fNuMIEventNtuple->Branch("flux_tptype",  &flux_tptype,  "flux_tptype/I");
    fNuMIEventNtuple->Branch("flux_vx",      &flux_vx,      "flux_vx/D");
    fNuMIEventNtuple->Branch("flux_vy",      &flux_vy,      "flux_vy/D");
    fNuMIEventNtuple->Branch("flux_vz",      &flux_vz,      "flux_vz/D");
    fNuMIEventNtuple->Branch("flux_ndecay",  &flux_ndecay,  "flux_ndecay/I");
    fNuMIEventNtuple->Branch("flux_ppmedium",&flux_ppmedium,"flux_ppmedium/I");

  }
   
  void LArNuMIana::beginSubRun(const art::SubRun& r)
  {
    art::Handle< sumdata::POTSummary > potListHandle;
    if ( r.getByLabel(fPOTModuleLabel,potListHandle) )
      fPOT->Fill(potListHandle->totpot);
    else
      std::cout << "POT BAD" << std::endl;
  }

  void LArNuMIana::analyze(const art::Event& event) 
  {
    art::ServiceHandle<geo::Geometry> geom;
      
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector< art::Ptr<simb::MCTruth> > mclist;
    if ( event.getByLabel(fGenieGenModuleLabel,mctruthListHandle) )
      art::fill_ptr_vector(mclist,mctruthListHandle);
    art::Ptr<simb::MCTruth> mctruth = mclist[0];

    art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
    std::vector< art::Ptr<simb::MCFlux> > mcfluxvec;
    if ( event.getByLabel(fGenieGenModuleLabel,mcfluxListHandle) )
      art::fill_ptr_vector(mcfluxvec,mcfluxListHandle);
    art::Ptr<simb::MCFlux> mcflux = mcfluxvec[0];
    
    if ( mctruth->NeutrinoSet() ) {
      fEvent  = event.id().event(); 
      fRun    = event.run();
      fSubRun = event.subRun();
      
      fNuPdgCode = mctruth->GetNeutrino().Nu().PdgCode();
      fNuVx      = mctruth->GetNeutrino().Nu().Vx();
      fNuVy      = mctruth->GetNeutrino().Nu().Vy();
      fNuVz      = mctruth->GetNeutrino().Nu().Vz();
      fNuPx      = mctruth->GetNeutrino().Nu().Px();
      fNuPy      = mctruth->GetNeutrino().Nu().Py();
      fNuPz      = mctruth->GetNeutrino().Nu().Pz();
      fNuE       = mctruth->GetNeutrino().Nu().E();

      fLeptonPdgCode  = mctruth->GetNeutrino().Lepton().PdgCode();
      fLeptonPx       = mctruth->GetNeutrino().Lepton().Px();
      fLeptonPy       = mctruth->GetNeutrino().Lepton().Py();
      fLeptonPz       = mctruth->GetNeutrino().Lepton().Pz();
      fLeptonE        = mctruth->GetNeutrino().Lepton().E();
      fLeptonVx       = mctruth->GetNeutrino().Lepton().Vx();
      fLeptonVy       = mctruth->GetNeutrino().Lepton().Vy();
      fLeptonVz       = mctruth->GetNeutrino().Lepton().Vz();
      fLeptonThetaXZ2 = std::atan2(fLeptonPy,fLeptonPz);
      fLeptonThetaYZ2 = std::atan2(fLeptonPx,fLeptonPz);
      fLeptonThetaXZ  = std::atan(fLeptonPy/fLeptonPz);
      fLeptonThetaYZ  = std::atan(fLeptonPx/fLeptonPz);
      
      flux_run      = mcflux->frun;
      flux_evtno    = mcflux->fevtno;
      flux_tpx      = mcflux->ftpx;
      flux_tpy      = mcflux->ftpy;
      flux_tpz      = mcflux->ftpz;
      flux_tptype   = mcflux->ftptype;
      flux_vx       = mcflux->fvx;
      flux_vy       = mcflux->fvy;
      flux_vz       = mcflux->fvz;
      flux_ndecay   = mcflux->fndecay;
      flux_ppmedium = mcflux->fppmedium;

      fCCint       = false;
      fNCint       = false;
      fCCQEint     = false;
      fNCQEint     = false;

      if ( mctruth->GetNeutrino().CCNC() == simb::kCC ) {
	fCCint = true;
	if ( mctruth->GetNeutrino().Mode() == 0 )
	  fCCQEint = true;
      }
      
      if ( mctruth->GetNeutrino().CCNC() == simb::kNC ) {
	fNCint = true;
	if ( mctruth->GetNeutrino().Mode() == 0 )
	  fNCQEint = true;
      }
      

      art::Handle< std::vector<simb::MCParticle> > particleHandle;
      event.getByLabel(fSimulationProducerLabel, particleHandle);
      //      std::map< int, const simb::MCParticle* > particleMap;
      
      int n_traj_points;      
      std::string tpc_string = "TPC";
      std::string vol_string;
      bool in_tpc = false;

      for ( auto const& particle : (*particleHandle) ) {
	//	particleMap[particle.TrackId()] = &particle; 

	n_traj_points = particle.NumberTrajectoryPoints();
      
	in_tpc = false;
	for ( int ii = 0; ii < n_traj_points; ++ii ) {
	  vol_string = geom->VolumeName(particle.Position(ii).Vect());
	  //std::cout << vol_string      << std::endl;
	  if ( vol_string.find(tpc_string) != std::string::npos )
	    in_tpc = true;
	}

	fInTPC.push_back(in_tpc);
		
	fProcess.push_back(particle.Process());
	fProdMaterial.push_back(geom->MaterialName(particle.Position(0).Vect()));
	fProdVolume.push_back(geom->VolumeName(particle.Position(0).Vect()));
	fEndMaterial.push_back(geom->MaterialName(particle.EndPosition().Vect()));
	fEndVolume.push_back(geom->VolumeName(particle.EndPosition().Vect()));

	fTrackID.push_back(particle.TrackId());
	fPdgCode.push_back(particle.PdgCode());
	fStartVx.push_back(particle.Vx(0));
	fStartVy.push_back(particle.Vy(0));
	fStartVz.push_back(particle.Vz(0));
	fStartPx.push_back(particle.Px(0));
	fStartPy.push_back(particle.Py(0));
	fStartPz.push_back(particle.Pz(0));
	fStartE.push_back(particle.E(0));

      } // loop over all particles in the event. 
      
    } // if neutrino set
    
    fNuMIEventNtuple->Fill();

    fProcess.clear();
    fProdMaterial.clear();
    fEndMaterial.clear();
    fProdVolume.clear();
    fEndVolume.clear();
    fTrackID.clear();
    fParentID.clear();
    fPdgCode.clear();
    fStartVx.clear();
    fStartVy.clear();
    fStartVz.clear();
    fStartPx.clear();
    fStartPy.clear();
    fStartPz.clear();
    fStartE.clear(); 
    fInTPC.clear();
  }

  DEFINE_ART_MODULE(LArNuMIana)

}

#endif
