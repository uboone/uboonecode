////////////////////////////////////////////////////////////////////////
// Class:       DQMTrackModule
// Module Type: analyzer
// File:        DQMTrackModule_module.cc
//
////////////////////////////////////////////////////////////////////////

/*!
 * Title:   DQMTrackModule
 * Author:  pawel.guzowski@manchester.ac.uk
 * Inputs:  recob::Wire (calibrated), recob::Hit, Assns<recob::Wire, recob::Hit>
 * Outputs: validation histograms
 *
 * Description:
 * This module is intended to be yet another hit analyzer module. Its intention is
 *   (1) to compare hit-finding modules against each other, and eventually
 *   (2) to compare those to truth
 */

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

#include "lardataobj/RecoBase/Track.h"

#include "TH1F.h"

namespace dqm {
  class DQMTrackModule;
}

class dqm::DQMTrackModule : public art::EDAnalyzer {
public:
  explicit DQMTrackModule(fhicl::ParameterSet const & p);
  virtual ~DQMTrackModule();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void endJob() override;

private:
  // Declare member data here.
  std::vector<std::string> fTrackModuleLabels;

  TH1F *NTracks;

  unsigned long sumNTracks;
  unsigned long sumNTracks2; // squared
  unsigned long sumNTracks3; // cubed
  unsigned long nEvents;

};


dqm::DQMTrackModule::DQMTrackModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ this->reconfigure(p); }

dqm::DQMTrackModule::~DQMTrackModule()
{
  // Clean up dynamic memory and other resources here.
}


void dqm::DQMTrackModule::analyze(art::Event const & e)
{
  //get the track data
  //
  //
  bool noTracks = true;
  art::Handle<std::vector<recob::Track> > trackHandle;

  for(auto trackModuleLabel: fTrackModuleLabels) {
    if(e.getByLabel(trackModuleLabel, trackHandle)) {
      unsigned long n = trackHandle->size();
      NTracks->Fill(n);
      sumNTracks += n;
      sumNTracks2 += n*n;
      sumNTracks3 += n*n*n;
      nEvents += 1;
      noTracks = false;
      break;
    }

  }

  if(noTracks) {
    // raise error that no reconstructed track containers were found
  }

}

void dqm::DQMTrackModule::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  NTracks = tfs->make<TH1F>("NTracks", "# tracks ", 150, 0, 150);

  sumNTracks = 0;
  sumNTracks2 = 0;
  sumNTracks3 = 0;
  nEvents = 0;
  
}

void dqm::DQMTrackModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

  fTrackModuleLabels  = p.get< std::vector<std::string> >("TrackModuleLabels");

}

void dqm::DQMTrackModule::endJob()
{
  double nev = (double) nEvents;
  double mean = sumNTracks / nev;
  double variance = (sumNTracks2 - nev * mean * mean) / nev;
  double stdev = std::sqrt(variance);
  double skew = (sumNTracks3 - 3. * mean * sumNTracks2 + 2. * nev * mean * mean * mean) / stdev / stdev / stdev;
  std::cout << "Mean: " << mean << " Variance: " << variance << " Skew: " << skew << std::endl;
}

DEFINE_ART_MODULE(dqm::DQMTrackModule)
