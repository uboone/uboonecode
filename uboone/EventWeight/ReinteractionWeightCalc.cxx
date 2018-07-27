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
#include "WeightCalcCreator.h"
#include "WeightCalc.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "CLHEP/Random/RandGaussQ.h"

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
      TDirectory* dir = dynamic_cast<TDirectory*>(probFile->Get(name.c_str()));
      assert(dir);
      pint = dynamic_cast<TH1F*>(dir->Get("pint"));
      assert(pint);
    }

    std::string name;  //!< String name
    int pdg;  //!< PDG code
    float par_sigma;  //!< Variation sigma set by user
    TH1F* pint;  //!< Interaction probability as a function of momentum
    std::vector<double> sigmas;  //!< Sigmas for universes
  };

private:
  std::string fMCParticleProducer;  //!< Label for the MCParticle producer
  std::string fMCTruthProducer;  //!< Label for the MCTruth producer
  CLHEP::RandGaussQ* fGaussRandom;  //!< Random number generator
  TFile* fProbFile;  //!< File with interaction probabilities, uncertainties
  std::map<int, ParticleDef> particles;  //!< Particles to reweight
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
  fNsims = pset.get<int> ("number_of_multisims", 0);

  // Prepare random generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));    

  // Load interaction probabilities
  fProbFile = TFile::Open("pint.root");  // FIXME
  assert(fProbFile);

  // Build parameter list
  for (size_t i=0; i<pars.size(); i++) {
    if (pars[i] == "p") {
      particles[2212] = ParticleDef("p",   2212, sigmas[i], fProbFile);
    }
    else if (pars[i] == "pip") {
      particles[211]  = ParticleDef("pip",  211, sigmas[i], fProbFile);
    }
    else if (pars[i] == "pim") {
      particles[-211] = ParticleDef("pim", -211, sigmas[i], fProbFile);
    }
    else {
      std::cerr << "Unknown particle type: " << pars[i] << std::endl;
      assert(false);
    }
  };

  // Set up universes
  for (auto& it : particles) {
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
  for (size_t itruth=0; itruth<truthHandle->size(); itruth++) {
    // Initialize weight vector for this MCTruth
    weight[itruth].resize(fNsims, 1.0);

    // Loop over MCParticles in the event
    auto const& mcparticles = truthParticles.at(itruth);

    for (size_t i=0; i<mcparticles.size(); i++) {
      int pdg = mcparticles.at(i)->PdgCode();
      double momentum = mcparticles.at(i)->Momentum(0).Vect().Mag();
      std::string endProc = mcparticles.at(i)->EndProcess();
      bool interacted = (endProc.find("Inelastic") != std::string::npos);

      // Reweight particles under consideration
      if (particles.find(pdg) != particles.end()) {
        ParticleDef& def = particles[pdg];
        int pbin = def.pint->FindBin(momentum);

        // Central value
        float pstop = TMath::Exp(-def.pint->Integral(0, pbin));

        // Loop through universes
        for (size_t j=0; j<weight[0].size(); j++) {
          // Integrate the modified interaction probability
          float integral = 0.0;
          for (int k=0; k<pbin; k++) {
            float s = 1.0 - (1.0 - def.pint->GetBinError(k)) * def.sigmas[j];
            integral += def.pint->GetBinContent(k) * s;
          }

          double pstopPrime = TMath::Exp(-integral);

          // Compute event weight as ratio of interaction probabilities
          float w = interacted ? (1.0 - pstopPrime) / (1.0 - pstop) : (pstopPrime / pstop);

          // Total weight is the product of track weights in the event
          weight[itruth][j] *= w;
        }
      }
    }
  }

  return weight;
}

REGISTER_WEIGHTCALC(ReinteractionWeightCalc)

}  // namespace evwgh

