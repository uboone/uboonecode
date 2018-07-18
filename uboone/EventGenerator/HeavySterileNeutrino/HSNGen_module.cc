////////////////////////////////////////////////////////////////////////
/// \file  HSNGen_module.cc
/// \brief Generator for heavy sterile neutrinos based on pre-generated HSN fluxes.
///
/// Generator for heavy sterile neutrino decays inside the MicroBooNE detector.
/// The code is largely based on InFlight generator, an event generator for sterile decays at
/// SBL facilities independent of LArSoft, written by Mark Ross-Lonergan and Peter Ballett.
/// Credit for the code goes to the two original authors.
///
/// Any malfunction, instability or bug is to be attributed solely to
/// Salvatore Davide Porzio, responsible for re-writing the code in its current form
/// and porting it to the LArSoft code.
///
/// \author  salvatore.porzio@postgrad.manchester.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_HSNGen_H
#define EVGEN_HSNGen_H

// ROOT includes
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TString.h"
#include "TSystem.h" //need BaseName and DirName

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcore/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcoreobj/SummaryData/RunData.h"

// HSNGen includes
#include "DataObjects/FourMomentum.h"
#include "DataObjects/SterileNeutrino.h"
#include "DataObjects/Flux.h"
#include "DataObjects/Observables.h"
#include "DataObjects/Channel.h"
#include "Helpers/Helper.h"
#include "Helpers/Settings.h"

// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <sqlite3.h> 
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoissonQ.h"
// #include "ifdh.h"  //to handle flux files

namespace hsngen
{
  /// A module to check the results from the Monte Carlo generator
  class HSNGen : public art::EDProducer
  {
  public:
    explicit HSNGen(fhicl::ParameterSet const& pset);
    virtual ~HSNGen();                       
    void produce(art::Event& evt);  
    void beginJob();
    void endJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
    void ClearData();
  private:
    // Fcl settings
    bool fPrintHepEvt;
    double fSterileMass;
    int fDecayChannel;
    std::string fFluxFile;
    double fDistance;
    double fGlobalTimeOffset;
    double fBeamWindow;
    std::vector<double> fBoundariesX;
    std::vector<double> fBoundariesY;
    std::vector<double> fBoundariesZ;
    std::vector<double> fGeneratedTimeWindow;

    // Analysis variables
    Settings gSett;
    std::vector<double> gModelParams;
    twoIP_channel *gChan;
    FluxFile gFlux;
    int gFakeRunNumber;

    // Diagnostic tree
    TTree *tTree;
    int run, subrun, event;
    std::vector<int> pdgCode;
    std::vector<double> Vx, Vy, Vz, T, Px, Py, Pz, E, P, Pt;
    double OpeningAngle, InvariantMass;

    // Auxiliary functions
    void CompressSettings(Settings &set);
  }; // END class HSNGen

  HSNGen::HSNGen(fhicl::ParameterSet const& p)
  : fPrintHepEvt(p.get<bool>("PrintHepEvt")),
    fSterileMass(p.get<double>("SterileMass")),
    fDecayChannel(p.get<double>("DecayChannel")),
    fFluxFile(p.get<std::string>("FluxFile")),
    fDistance(p.get<double>("Distance")),
    fGlobalTimeOffset(p.get<double>("GlobalTimeOffset")),
    fBeamWindow(p.get<double>("BeamWindow")),
    fBoundariesX(p.get<std::vector<double>>("BoundariesX")),
    fBoundariesY(p.get<std::vector<double>>("BoundariesY")),
    fBoundariesZ(p.get<std::vector<double>>("BoundariesZ")),
    fGeneratedTimeWindow(p.get<std::vector<double>>("GeneratedTimeWindow"))
  {
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });
    this->reconfigure(p);

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >(); 
    return;
  }
  
  HSNGen::~HSNGen()
  {}

  void HSNGen::reconfigure(fhicl::ParameterSet const& p)
  {

    return;
  }

  void HSNGen::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    tTree = tfs->make<TTree>("Data","");
    tTree->Branch("run",&run,"run/I");
    tTree->Branch("subrun",&subrun,"subrun/I");
    tTree->Branch("event",&event,"event/I");
    tTree->Branch("Vx",&Vx);
    tTree->Branch("Vy",&Vy);
    tTree->Branch("Vz",&Vz);
    tTree->Branch("T",&T);
    tTree->Branch("Px",&Px);
    tTree->Branch("Py",&Py);
    tTree->Branch("Pz",&Pz);
    tTree->Branch("E",&E);
    tTree->Branch("P",&P);
    tTree->Branch("Pt",&Pt);
    tTree->Branch("OpeningAngle",&OpeningAngle);
    tTree->Branch("InvariantMass",&InvariantMass);

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &gEngine = rng->getEngine("gen");
    CompressSettings(gSett);
    FillModel(gEngine, gChan, gModelParams, gSett);
    gFlux = FluxFile(fFluxFile, fSterileMass);
    gFakeRunNumber = 0; // Used for the Hepevt format output
    return;
  }

  void HSNGen::endJob()
  {
    delete gChan;
    return;
  }


  void HSNGen::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    run.put(std::move(runcol));
    return;
  }

  void HSNGen::ClearData()
  {
    run = -1;
    subrun = -1;
    event = -1;
    OpeningAngle = -999;
    InvariantMass = -999;
    pdgCode.clear();
    Vx.clear();
    Vy.clear();
    Vz.clear();
    T.clear();
    Px.clear();
    Py.clear();
    Pz.clear();
    E.clear();
    P.clear();
    Pt.clear();
  } // END function ClearData

  void HSNGen::produce(art::Event& evt)
  {
    // Clear diagnostic vectors
    ClearData();
    // Declare geometry and MCTruth objects
    art::ServiceHandle<geo::Geometry> geom;
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    truth.SetOrigin(simb::kUnknown);

    // Generate observables characterizing the event
    double neutrinoTime = -1;
    Observables obs;
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &gEngine = rng->getEngine("gen");

    while (neutrinoTime <= fGeneratedTimeWindow[0] || neutrinoTime >=fGeneratedTimeWindow[1])
    {
      GenerateObservables(gEngine, gChan, gFlux, gSett, obs);
      neutrinoTime = obs.time;
    }

    // Generate MCParts from observables
    simb::MCParticle p1(0,obs.pdg1,"primary",-1,obs.mass1,1);
    simb::MCParticle p2(1,obs.pdg2,"primary",-1,obs.mass2,1);
    TLorentzVector pos(obs.xPos, obs.yPos, obs.zPos, obs.time);
    TLorentzVector mom1(obs.P1[0],obs.P1[1],obs.P1[2],obs.E1);
    TLorentzVector mom2(obs.P2[0],obs.P2[1],obs.P2[2],obs.E2);
    p1.AddTrajectoryPoint(pos,mom1);
    p2.AddTrajectoryPoint(pos,mom2);
    truth.Add(p1);
    truth.Add(p2);
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));

    // Fill tree variables
    pdgCode.push_back(p1.PdgCode());
    Vx.push_back(p1.Vx());
    Vy.push_back(p1.Vy());
    Vz.push_back(p1.Vz());
    T.push_back(p1.T());
    Px.push_back(p1.Px());
    Py.push_back(p1.Py());
    Pz.push_back(p1.Pz());
    E.push_back(p1.E());
    P.push_back(p1.P());
    Pt.push_back(p1.Pt());
    pdgCode.push_back(p2.PdgCode());
    Vx.push_back(p2.Vx());
    Vy.push_back(p2.Vy());
    Vz.push_back(p2.Vz());
    T.push_back(p2.T());
    Px.push_back(p2.Px());
    Py.push_back(p2.Py());
    Pz.push_back(p2.Pz());
    E.push_back(p2.E());
    P.push_back(p2.P());
    Pt.push_back(p2.Pt());
    double dotProduct = Px[0]*Px[1] + Py[0]*Py[1] + Pz[0]*Pz[1];
    OpeningAngle = dotProduct / float(P[0]*P[1]);
    double eTerm = pow((E[0] + E[1]),2.);
    double pTerm = pow(P[0],2.) + pow(P[1],2.) + 2.*dotProduct;
    InvariantMass = sqrt(eTerm - pTerm);
    tTree->Fill();

    gFakeRunNumber++;
    return;
  } // END function produce

  // Compress fcl settings to a struct in order to make it easier
  // to pass them to other functions
  void HSNGen::CompressSettings(Settings &set)
  {
    set.printHepEvt = fPrintHepEvt;
    set.sterileMass = fSterileMass;
    set.decayChannel = fDecayChannel;
    set.fluxFile = fFluxFile;
    set.distance = fDistance;
    set.globalTimeOffset = fGlobalTimeOffset;
    set.beamWindow = fBeamWindow;
    set.boundariesX = fBoundariesX;
    set.boundariesY = fBoundariesY;
    set.boundariesZ = fBoundariesZ;
    set.generatedTimeWindow = fGeneratedTimeWindow;
    return;
  }

  DEFINE_ART_MODULE(HSNGen)

} // END namespace hsngen

#endif
