////////////////////////////////////////////////////////////////////////
// Class:       QuietEventAna
// Plugin Type: analyzer (art v2_05_00)
// File:        QuietEventAna_module.cc
//
// Generated at Thu Feb 16 01:04:06 2017 by Taritree Wongjirad using cetskelgen
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

#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/Wire.h"

#include "TTree.h"

class QuietEventAna;


class QuietEventAna : public art::EDAnalyzer {
public:
  explicit QuietEventAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  QuietEventAna(QuietEventAna const &) = delete;
  QuietEventAna(QuietEventAna &&) = delete;
  QuietEventAna & operator = (QuietEventAna const &) = delete;
  QuietEventAna & operator = (QuietEventAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
  int m_Run;
  int m_SubRun;
  int m_Event;
  float m_ADCSum;

  std::string fWireModuleLabel;
  float fThreshold;
  TTree* fTree;

};


QuietEventAna::QuietEventAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  fWireModuleLabel = p.get<std::string>("WireModuleLabel");
  fThreshold       = p.get<float>("Threshold");

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("QuietEvents","quiet event tree");
  fTree->Branch( "run", &m_Run, "run/I" );
  fTree->Branch( "subrun", &m_SubRun, "subrun/I" );
  fTree->Branch( "event", &m_Event, "event/I" );
  fTree->Branch( "adcsum", &m_ADCSum, "adcsum/F" );

}

void QuietEventAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  // We do a dumb thing, we grab the Y-plane and sum the total ADC values.
  // We'll store the info in a ROOT file. Fun.
  m_Run    = e.run();
  m_Event  = e.event();
  m_SubRun = e.subRun();

  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireModuleLabel,wireHandle);

  std::vector<recob::Wire> const& wireVector(*wireHandle);

  m_ADCSum = 0.;
  for ( auto const& wire : wireVector ) {
    if ( wire.View()!=geo::kZ ) continue;

    std::vector<float> signal = wire.Signal();

    for ( auto &val : signal ) {
      if ( val>fThreshold ) {
	m_ADCSum += val;
      }
    }

  }

  fTree->Fill();

}

DEFINE_ART_MODULE(QuietEventAna)
