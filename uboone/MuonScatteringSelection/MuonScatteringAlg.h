#ifndef MUONSCATTERINGALG_H
#define MUONSCATTERINGALG_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <vector>
#include <string>
#include <exception>

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "TTree.h"
#include "TH1F.h"

namespace andy{

  class MuonScatteringAlgException : public std::exception{
    virtual const char* what() const throw(){
      return "MuonScatteringAlg Exception";
    }
  } muonscatteringexception;

  class MuonScatteringAlg{

  public:

    MuonScatteringAlg();
    //~MuonScatteringAlg();

    bool inFV(double xyz[3]);
    bool fullyContained(recob::Track t);
    double closestApproach(recob::Track t, double xyz[3]);
    double absDistance(double * a, double * b);
    TVector3 TPCentryPoint(const simb::MCParticle & part);
//    simb::MCParticle GetMCParticleFromRecoTrack(recob::Track t, art::Event const & e);

  private:
    const double FVx; //(256.35);
    const double FVy; //(233.);
    const double FVz; //(1036.8);
    const double borderx; //(10);
    const double bordery; //(10);
    const double borderz; //(40);
    //friend class MuonScatteringAlgTest;

  };

}//end namespace hit


#endif
