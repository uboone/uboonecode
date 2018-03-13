////////////////////////////////////////////////////////////////////////
// Class:       TestModule
// Plugin Type: analyzer (art v2_05_01)
// File:        TestModule_module.cc
//
// Generated at Mon Mar 12 16:28:34 2018 by Adam Lister using cetskelgen
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

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "TTree.h"

class TestModule;


class TestModule : public art::EDAnalyzer {
  public:
    explicit TestModule(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TestModule(TestModule const &) = delete;
    TestModule(TestModule &&) = delete;
    TestModule & operator = (TestModule const &) = delete;
    TestModule & operator = (TestModule &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // selected optional functions
    void beginJob() override;

    int    evtno;
    double vtxx, vtxy, vtxz;
    double px, py, pz, E;
    int pdg;
    int ptype;
    double wgt;
    double dist;
    double nenergyn;
    double nwtnear;
    double nimpwt;

  private:
    TTree* flux;
    art::ServiceHandle< art::TFileService > tfs;

};


TestModule::TestModule(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{}

void TestModule::beginJob()
{

  flux = tfs->make<TTree>("flux", "flux");
  flux->Branch("evtno", &evtno);
  flux->Branch("vtxx" ,  &vtxx);
  flux->Branch("vtxy" ,  &vtxy);
  flux->Branch("vtxz" ,  &vtxz);
  flux->Branch("px"   ,  &px);
  flux->Branch("py"   ,  &py);
  flux->Branch("pz"   ,  &pz);
  flux->Branch("E"    ,  &E);
  flux->Branch("pdg"  ,  &pdg);
  flux->Branch("ptype",  &ptype);
  flux->Branch("wgt"  ,  &wgt);
  flux->Branch("dist" ,  &dist);
  flux->Branch("nenergyn", &nenergyn);
  flux->Branch("nwtnear", &nwtnear);
  flux->Branch("nimpwt", &nimpwt);

}

void TestModule::analyze(art::Event const & e)
{

  art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  e.getByLabel("generator", mcFluxHandle);
  if (!mcFluxHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;
  art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
  if (mcFluxVec.size() == 0){
    std::cout << ">> No MCFlux information" << std::endl;
    return;
  }

  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
  if (mcTruthVec.size() == 0){
    std::cout << ">> No MCTruth information" << std::endl;
    return;
  }

  const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
  const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const simb::MCParticle& nu = mcTruth->GetNeutrino().Nu();

  // possibly the wrong variables, but let's see for now...
  evtno = mcFlux->fevtno;
  vtxx  = mcFlux->fgenx;
  vtxy  = mcFlux->fgeny;
  vtxz  = mcFlux->fgenz;
  px    = nu.Px(); 
  py    = nu.Py(); 
  pz    = nu.Pz(); 
  E     = nu.E();
  pdg   = mcFlux->fntype;
  ptype = mcFlux->fptype;
  wgt   = 1; //mcFlux->fnwtnear;
  dist  = mcFlux->fdk2gen;
  nenergyn = mcFlux->fnenergyn;
  nwtnear = mcFlux->fnwtnear;
  nimpwt = mcFlux->fnimpwt;


  std::cout << "fgenx: " << mcFlux->fgenx << std::endl;
  std::cout << "fgeny: " << mcFlux->fgeny << std::endl;
  std::cout << "fgenz: " << mcFlux->fgenz << std::endl;
  std::cout << "mcpVx: " << nu.Vx() << std::endl;
  std::cout << "mcpVy: " << nu.Vy() << std::endl;
  std::cout << "mcpVz: " << nu.Vz() << std::endl;
  std::cout << "px   : " << px << std::endl;
  std::cout << "py   : " << py << std::endl;
  std::cout << "pz   : " << pz << std::endl;

  flux->Fill();

}

DEFINE_ART_MODULE(TestModule)
