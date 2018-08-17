/**
 * \class ReinteractionWeightCalc
 * \brief Hadron reinteraction event reweighting
 * \author A. Mastbaum <mastbaum@uchicago.edu>, 2018/07
 *
 * Reweight events based on hadron reinteraction probabilities.
 */

#include <map>
#include <string>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "WeightCalcCreator.h"
#include "WeightCalc.h"

namespace evwgh {

class ReinteractionWeightCalc : public WeightCalc {
public:
  ReinteractionWeightCalc() {}

  void Configure(fhicl::ParameterSet const& p);

  std::vector<std::vector<double> > GetWeight(art::Event& e);

  /**
   * \struct ParticleDef
   * \brief A reweightable particle definition
   */
  class ParticleDef {
  public:
    ParticleDef() {}
    ParticleDef(std::string _name, int _pdg, float _sigma, TFile* probFile)
        : name(_name), pdg(_pdg), par_sigma(_sigma) {
      pint = dynamic_cast<TH1F*>(probFile->Get(name.c_str()));
      assert(pint);

      // Reconstitute the cross section vs. KE
      char name[100];
      snprintf(name, 100, "_xs_%s", name);
      xs = dynamic_cast<TH1F*>(pint->Clone(name));

      for (int j=1; j<xs->GetNbinsX()+1; j++) {
        float p1 = pint->GetBinContent(j);
        float p2 = pint->GetBinContent(j-1);
        float v = 0;

        // TODO: enforced monotonicity
        if (p1 > p2 && p1 < 1) {
          v = -1.0 * log((1.0 - p1) / (1.0 - p2));
        }

        xs->SetBinContent(j, v);
      }
    }

    std::string name;  //!< String name
    int pdg;  //!< PDG code
    float par_sigma;  //!< Variation sigma set by user
    TH1F* pint;  //!< Interaction probability as a function of KE
    TH1F* xs;  //!< Derived effective cross section
    std::vector<double> sigmas;  //!< Sigmas for universes
  };

private:
  std::string fMCParticleProducer;  //!< Label for the MCParticle producer
  std::string fMCTruthProducer;  //!< Label for the MCTruth producer
  CLHEP::RandGaussQ* fGaussRandom;  //!< Random number generator
  TFile* fProbFile;  //!< File with interaction probabilities, uncertainties
  std::map<int, ParticleDef> fParticles;  //!< Particles to reweight
  unsigned fNsims;  //!< Number of multisims

  DECLARE_WEIGHTCALC(ReinteractionWeightCalc)
};


void ReinteractionWeightCalc::Configure(fhicl::ParameterSet const& p) {
  // Get configuration for this function
  fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>(GetName());
  fMCParticleProducer = pset.get<std::string>("MCParticleProducer", "largeant");
  fMCTruthProducer = pset.get<std::string>("MCTruthProducer", "generator");
  std::vector<std::string> pars = pset.get< std::vector<std::string> >("parameter_list");	
  std::vector<float> sigmas = pset.get<std::vector<float> >("parameter_sigma");	
  std::string mode = pset.get<std::string>("mode");
  std::string probFileName = pset.get<std::string>("ProbFileName", "systematics/reint/interaction_probabilities.root");
  fNsims = pset.get<int> ("number_of_multisims", 0);

  // Prepare random generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));    

  // Load interaction probabilities
  cet::search_path sp("FW_SEARCH_PATH");
  std::string probFilePath = sp.find_file(probFileName);
  fProbFile = TFile::Open(probFilePath.c_str());
  assert(fProbFile && fProbFile->IsOpen());

  // Build parameter list
  for (size_t i=0; i<pars.size(); i++) {
    if (pars[i] == "p") {
      fParticles[2212] = ParticleDef("p",   2212, sigmas[i], fProbFile);
    }
    else if (pars[i] == "pip") {
      fParticles[211]  = ParticleDef("pip",  211, sigmas[i], fProbFile);
    }
    else if (pars[i] == "pim") {
      fParticles[-211] = ParticleDef("pim", -211, sigmas[i], fProbFile);
    }
    else {
      std::cerr << "Unknown particle type: " << pars[i] << std::endl;
      assert(false);
    }
  };


  // Set up universes
  for (auto& it : fParticles) {
    if (mode == "pm1sigma") {
      // pm1sigma mode: 0 = +1sigma, 1 = -1sigma
      it.second.sigmas.push_back( 1.0);
      it.second.sigmas.push_back(-1.0);
      fNsims = 2;
    }
    else if (mode == "multisim") {
      // multisim mode: Scale factors sampled within the given uncertainty
      for (unsigned j=0; j<fNsims; j++) {
        double r = fGaussRandom->shoot(&rng->getEngine(GetName()), 0.0, 1.0);
        it.second.sigmas.push_back(it.second.par_sigma * r);
      }
    }
    else {
      // Anything else mode: Set scale to user-defined scale factor
      it.second.sigmas.push_back(it.second.par_sigma);
      fNsims = 1;
    }
  }
}


std::vector<std::vector<double> >
ReinteractionWeightCalc::GetWeight(art::Event& e) {
  // Get MCParticles for each MCTruth in this event
  art::Handle<std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer, truthHandle);
  const art::FindManyP<simb::MCParticle> truthParticles(truthHandle, e, fMCParticleProducer);
  assert(truthParticles.isValid());

  // Initialize the vector of event weights
  std::vector<std::vector<double> > weight(truthHandle->size());

  // Loop over sets of MCTruth-associated particles
  for (size_t itruth=0; itruth<truthParticles.size(); itruth++) {

    // Initialize weight vector for this MCTruth
    weight[itruth].resize(fNsims, 1.0);

    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);

    for (size_t i=0; i<mcparticles.size(); i++) {
      const simb::MCParticle& p = *mcparticles.at(i);
      int pdg = p.PdgCode();
      double ke = p.E() - p.Mass();
      std::string endProc = p.EndProcess();
      bool interacted = (endProc.find("Inelastic") != std::string::npos);

      // Reweight particles under consideration
      if (fParticles.find(pdg) != fParticles.end()) {
        ParticleDef& def = fParticles[pdg];
        int kebin = def.pint->FindBin(ke);

        // Loop through universes
        for (size_t j=0; j<weight[0].size(); j++) {
          // Integrate a modified cross section to find a survival probability
          float sprob = 1.0;
          for (int k=0; k<kebin; k++) {
            float wbin = 1.0 - (1.0 - 0.3 /*def.huncert->GetBinContent(k)*/) * def.sigmas[j];
            float xs = wbin * def.xs->GetBinContent(k);
            sprob *= exp(-1.0 * xs);
          }

          float w;
          if (interacted) {
            w = (1.0 - sprob) / def.pint->GetBinContent(kebin);
          }
          else {
            w = sprob / (1.0 - def.pint->GetBinContent(kebin));
          }

          // Total weight is the product of track weights in the event
          weight[itruth][j] *= std::max((float)0.0, w);
        }
      }
    }
  }

  return weight;
}

REGISTER_WEIGHTCALC(ReinteractionWeightCalc)

}  // namespace evwgh

