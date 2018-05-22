//
// Name: TriggerTiming_module.cc
//
// Purpose: Module TriggerTiming.
//
// Configuration parameters.
//
//  TriggerModuleLabel: Trigger module label.
//  PMTModuleLabel:     PMT (OpDetWaveform) module label.
//  PMTModuleInstance:  PMT (OpDetWaveform) module instance.
//
// Created: 18-Oct-2017  H. Greenlee
//

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "uboone/TriggerSim/UBTriggerTypes.h"

#include "TH1F.h"

namespace raw {

  class TriggerTiming : public art::EDAnalyzer
  {
  public:

    // Constructors, destructor

    explicit TriggerTiming(fhicl::ParameterSet const& pset);
    virtual ~TriggerTiming();

    // Overrides.

    void analyze(const art::Event& evt);
    void endJob();

  private:

    // Fcl parameters.

    std::string fTriggerModuleLabel;
    std::string fPMTModuleLabel;
    std::string fPMTModuleInstance;
    unsigned int fPMTMinLength;

    // Histograms.

    TH1F* fHdtbnb;    // BNB delta-t.
    TH1F* fHdtnumi;   // NUMI delta-t.
    TH1F* fHdtext;    // EXT delta-t.
    TH1F* fHdtmucs;   // MUCS delta-t.

    // Statistics.

    int fNumEvent;
  };

  DEFINE_ART_MODULE(TriggerTiming)

  TriggerTiming::TriggerTiming(const fhicl::ParameterSet& pset) :
  //
  // Purpose: Constructor.
  //
  // Arguments: pset - Module parameters.
  //
    EDAnalyzer(pset),
    fTriggerModuleLabel(pset.get<std::string>("TriggerModuleLabel")),
    fPMTModuleLabel(pset.get<std::string>("PMTModuleLabel")),
    fPMTModuleInstance(pset.get<std::string>("PMTModuleInstance")),
    fPMTMinLength(pset.get<unsigned int>("PMTMinLength")),
    fNumEvent(0)
  {
    // Book histograms.

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("TriggerTiming", "Trigger Timing");
    fHdtbnb = dir.make<TH1F>("dtbnb", "BNB Delta-t", 50, -0.5, 49.5);
    fHdtnumi = dir.make<TH1F>("dtnumi", "NUMI Delta-t", 50, -0.5, 49.5);
    fHdtext = dir.make<TH1F>("dtext", "EXT Delta-t", 50, -0.5, 49.5);
    fHdtmucs = dir.make<TH1F>("dtmucs", "MUCS Delta-t", 50, -0.5, 49.5);    

    // Report.

    mf::LogInfo("TriggerTiming") 
      << "TriggerTiming configured with the following parameters:\n"
      << "  TriggerModuleLabel = " << fTriggerModuleLabel << "\n"
      << "  PMTModuleLabel = " << fPMTModuleLabel << "\n"
      << "  PMTModuleInstance = " << fPMTModuleInstance << "\n"
      << "  PMTMinLength = " << fPMTMinLength;
  }

  TriggerTiming::~TriggerTiming()
  //
  // Purpose: Destructor.
  //
  {}

  void TriggerTiming::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Get services.

    art::ServiceHandle<art::TFileService> tfs;
    //const detinfo::DetectorClocks* dc = lar::providerFrom<detinfo::DetectorClocksService>();

    // Fetch raw::Trigger data product.

    art::Handle<std::vector<raw::Trigger> > htrig;
    evt.getByLabel(fTriggerModuleLabel, htrig);
    if(htrig.isValid() and htrig->size() > 0) {
      const raw::Trigger& trig = htrig->front();
      unsigned int tmask = trig.TriggerBits();
      double trigger_time = trig.TriggerTime();
      //std::cout << "Trigger bit mask = 0x" << std::hex << tmask << std::dec << std::endl;
      //std::cout << "Trigger time = " << trigger_time << std::endl;      

      // Fetch OpDetWaveform data product.

      art::Handle<std::vector<raw::OpDetWaveform> > hpmt;
      evt.getByLabel(fPMTModuleLabel, fPMTModuleInstance, hpmt);
      if(hpmt.isValid()) {
	for(auto const& pmt : *hpmt) {
	  unsigned int channel = pmt.ChannelNumber() % 100;
	  unsigned int chlen = pmt.size();
	  if(channel < 32 && chlen >= fPMTMinLength) {
	    double timestamp = pmt.TimeStamp();
	    //std::cout << "PMT channel = " << channel << std::endl;
	    //std::cout << "PMT length = " << chlen << std::endl;
	    //std::cout << "PMT time stamp = " << timestamp << std::endl;      
	    //std::cout << "PMT-trigger time difference = " << timestamp - trigger_time << std::endl;
	    // Fill histograms.

	    double dtbins = 64. * (timestamp - trigger_time);
	    if(tmask & (1 << trigger::kTriggerBNB))
	      fHdtbnb->Fill(dtbins);
	    if(tmask & (1 << trigger::kTriggerNuMI))
	      fHdtnumi->Fill(dtbins);
	    if(tmask & (1 << trigger::kTriggerEXT))
	      fHdtext->Fill(dtbins);
	    if(tmask & (1 << trigger::kSpare))
	      fHdtmucs->Fill(dtbins);
	  }
	}
      }
    }    
  }

  void TriggerTiming::endJob()
  //
  // Purpose: End of job.
  //
  {
    // Print summary.

    mf::LogInfo("TriggerTiming") 
      << "TriggerTiming statistics:\n"
      << "  Number of events = " << fNumEvent;
  }
}
