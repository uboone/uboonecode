// Analyzer module that will dump GENIE events from an art-format ROOT file
// into a ROOT file that resembles one produced by gevgen.
//
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <memory>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "cetlib/maybe_ref.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"

// GENIE includes
#ifdef GENIE_PRE_R3
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpWriter.h"
#else
// Use these includes for GENIE v3
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpWriter.h"
#endif

namespace sim {
  class GENIEExtractor;
}

class sim::GENIEExtractor : public art::EDAnalyzer {

public:

  explicit GENIEExtractor(const fhicl::ParameterSet& p);
  virtual ~GENIEExtractor();

  void analyze(const art::Event& evt) override;

protected:

  std::string fGenieLabel;
  std::string fOutputFilename;
  int fEventCount;

  genie::NtpWriter fWriter;
};

sim::GENIEExtractor::GENIEExtractor(const fhicl::ParameterSet& pset) : art::EDAnalyzer(pset),
  fEventCount(0)
{
  //#ifndef GENIE_PRE_R3
  //// Make sure the right tune is set up for GENIE v3+
  //std::string genie_tune = pset.get<std::string>("tune_name", "${GENIE_XSEC_TUNE}");
  //evgb::SetEventGeneratorListAndTune("", genie_tune);
  //#endif

  fGenieLabel = pset.get<std::string>("genie_label", "generator");
  fOutputFilename = pset.get<std::string>("output_filename", "genie_dump.ghep.root");

  fWriter.CustomizeFilename( fOutputFilename );
  fWriter.Initialize();
}

sim::GENIEExtractor::~GENIEExtractor()
{
  fWriter.Save();
}

void sim::GENIEExtractor::analyze(const art::Event& evt) {

  // Find each of the MCTruth/GTruth pairs created by GENIE in the current event
  const auto& mc_truths_handle = evt.getValidHandle< std::vector<simb::MCTruth> >(
    fGenieLabel );
  art::FindOne<simb::GTruth> g_truths_find(mc_truths_handle, evt, fGenieLabel);

  const auto& mc_truths = *mc_truths_handle;

  for ( size_t k = 0u; k < mc_truths.size(); ++k ) {

    // Use each pair of MCTruth and GTruth objects to reconstruct the
    // GENIE event record
    const simb::MCTruth& mct = mc_truths.at( k );
    const simb::GTruth& gt = g_truths_find.at( k ).ref();

    std::unique_ptr<genie::EventRecord> ev_rec( evgb::RetrieveGHEP(mct, gt) );

    // Write it to the output file in GENIE's native format
    fWriter.AddEventRecord( fEventCount, ev_rec.get() );
    ++fEventCount;
  }
  //MF_LOG_INFO("sim::GENIEExtractor") << "";
}

DEFINE_ART_MODULE(sim::GENIEExtractor)
