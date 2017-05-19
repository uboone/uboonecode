////////////////////////////////////////////////////////////////////////
// Class:       FitterTree
// Plugin Type: analyzer (art v2_06_03)
// File:        FitterTree_module.cc
//
// Generate a skimmed TTree for output to a fitter.
//
// Generated at Thu May 18 09:55:36 2017 by Andrew Mastbaum using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include <map>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"

class FitterTree : public art::EDAnalyzer {
public:
  explicit FitterTree(fhicl::ParameterSet const & p);

  FitterTree(FitterTree const &) = delete;
  FitterTree(FitterTree &&) = delete;
  FitterTree& operator = (FitterTree const &) = delete;
  FitterTree& operator = (FitterTree &&) = delete;

  void analyze(art::Event const & e) override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p) override;

  /**
   * \struct IntData
   * \brief A simple struct to store interaction data for filling the output
   * ROOT tree.
   */
  struct IntData {
    // Interaction
    int type;  //!< MCNeutrino.InteractionType
    int target;  //!< MCNeutrino.Target
    evwgh::MCEventWeight mcweight;  //!< MCEventWeight

    // Neutrino
    int nuPDG;  //!< MCNeutrino.Nu.PdgCode
    TLorentzVector nuVtx;  //!< MCNeutrino.Nu.Position
    TLorentzVector nuMom;  //!< MCNeutrino.Nu.Momentum
    double nuEnergy;  //!< MCNeutrino.Nu.E
    double nuQ2;  //!< MCNeutrino.Qsqr
    double nuPt;  //!< MCNeutrino.Pt
    double nuTheta;  //!< MCNeutrino.Theta

    // Lepton
    int leptonPDG;  //!< MCNeutrino.Lepton.PdgCode
    TLorentzVector leptonVtx;  //!< MCNeutrino.Lepton.Position
    TLorentzVector leptonMom;  //!< MCNeutrino.Lepton.Momentum
    double leptonEnergy;  //!< MCNeutrino.Lepton.E
  };

private:
  std::string fGeneratorLabel;  // Label for the truth generator
  std::string fWeightLabel;  // Label for the MC weight generator
  TFile* fFile;  //!< Output ROOT file
  TTree* fTree;  //!< Output ROOT tree
  IntData fEvent;  //!< Event data for the output tree
};


FitterTree::FitterTree(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
  this->reconfigure(p);
}


void FitterTree::reconfigure(fhicl::ParameterSet const& p) {
  art::ServiceHandle<art::TFileService> tfs;
  fFile = &(tfs->file());
  assert(fFile);

  std::string fGeneratorLabel = \
    p.get<std::string>("GeneratorLabel", "generator");

  std::string fWeightLabel = \
    p.get<std::string>("WeightLabel", "eventweight");

  std::string treeName = p.get<std::string>("TreeName", "fittree");

  fFile->cd();
  fTree = new TTree(treeName.c_str(), "Fitter Tree");
  fTree->Branch("type", &fEvent.type);
  fTree->Branch("target", &fEvent.target);
  fTree->Branch("mcweight", &fEvent.mcweight);
  fTree->Branch("nuPDG", &fEvent.nuPDG);
  fTree->Branch("nuVtx", &fEvent.nuVtx);
  fTree->Branch("nuMom", &fEvent.nuMom);
  fTree->Branch("nuEnergy", &fEvent.nuEnergy);
  fTree->Branch("nuQ2", &fEvent.nuQ2);
  fTree->Branch("nuPt", &fEvent.nuPt);
  fTree->Branch("nuTheta", &fEvent.nuTheta);
  fTree->Branch("leptonPDG", &fEvent.leptonPDG);
  fTree->Branch("leptonVtx", &fEvent.leptonVtx);
  fTree->Branch("leptonMom", &fEvent.leptonMom);
  fTree->Branch("leptonEnergy", &fEvent.leptonEnergy);
}


void FitterTree::analyze(art::Event const & e) {
  // Grab the MCTruth information
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel("generator", truthHandle);
  if (!truthHandle.isValid()) {
    std::cout << "truth handle not valid" << std::endl;
  }

  const std::vector<simb::MCTruth>* truthVector = truthHandle.product();
  if (truthVector->empty()) {
    std::cout << "truth vector is empty" << std::endl;
  }

  simb::MCTruth mcTruth = truthVector->at(0);

  // Grab the MCEventWeights
  art::Handle< std::vector< evwgh::MCEventWeight > > wghHandle;
  e.getByLabel("eventweight", wghHandle);

  if (!wghHandle.isValid()) {
    std::cout << "eventweight handle not valid" << std::endl;
    return;
  }

  evwgh::MCEventWeight mcWeight;
  const std::vector<evwgh::MCEventWeight>* wghv = wghHandle.product();
  if (wghv->empty()) {
    std::cout << "eventweight vector is empty" << std::endl;
  }
  else {
    mcWeight = wghv->at(0);
  }

  // Fill interaction truth
  simb::MCNeutrino const& nu = mcTruth.GetNeutrino();
  fEvent.type = nu.InteractionType();
  fEvent.target = nu.Target();
  fEvent.mcweight = mcWeight;
  fEvent.nuPDG = nu.Nu().PdgCode();
  fEvent.nuVtx = nu.Nu().Position();
  fEvent.nuMom = nu.Nu().Momentum();
  fEvent.nuEnergy = nu.Nu().E();
  fEvent.nuQ2 = nu.QSqr();
  fEvent.nuPt = nu.Pt();
  fEvent.nuTheta = nu.Theta();
  fEvent.leptonPDG = nu.Lepton().PdgCode();
  fEvent.leptonVtx = nu.Lepton().Position();
  fEvent.leptonMom = nu.Lepton().Momentum();
  fEvent.leptonEnergy = nu.Lepton().E();

  fFile->cd();
  fTree->Fill();
}


void FitterTree::endJob() {
  fFile->cd();
  fTree->Write();
}

DEFINE_ART_MODULE(FitterTree)

