////////////////////////////////////////////////////////////////////////
// Class:       DQMChannelNoiseModule
// Plugin Type: analyzer (art v2_07_03)
// File:        DQMChannelNoiseModule_module.cc
//
// Generated at Mon Jul 24 14:54:06 2017 by pguzowski using cetskelgen
// from cetlib version v3_00_01.
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
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include "TH1D.h"
#include "TH2D.h"

namespace dqm {
  class DQMChannelNoiseModule;
}

const uint32_t num_channels = 8256;


class dqm::DQMChannelNoiseModule : public art::EDAnalyzer {
public:
  explicit DQMChannelNoiseModule(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DQMChannelNoiseModule(DQMChannelNoiseModule const &) = delete;
  DQMChannelNoiseModule(DQMChannelNoiseModule &&) = delete;
  DQMChannelNoiseModule & operator = (DQMChannelNoiseModule const &) = delete;
  DQMChannelNoiseModule & operator = (DQMChannelNoiseModule &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  std::vector<unsigned short> fADCs[num_channels];
  unsigned long long fTotalTicks[num_channels];
  unsigned long long fTotalADC[num_channels];
  unsigned long long fTotalADCSquared[num_channels];

  std::string fDigitModuleLabel;
  std::vector<uint32_t> fChannelsToHistogram;
  uint32_t fMaxChannel;

  std::vector<TH1D *> fAdcHists;
  TH1D *fMeanHist, *fVarHist, *fRMSHist, *fLengthHist, *fPlaneHist;
};


dqm::DQMChannelNoiseModule::DQMChannelNoiseModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), fMaxChannel(0)  // ,
 // More initializers here.
{ this->reconfigure(p);
  for(uint32_t i =0; i < num_channels; ++i) {
    fTotalTicks[i] = fTotalADC[i] = fTotalADCSquared[i] = 0;
    //fADCs[i].resize(100000);
  }
}

void dqm::DQMChannelNoiseModule::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  //
  /*
  if(!(fLengthHist->GetEntries()>0.)) {
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    for(uint32_t chan = 0; chan < num_channels; ++chan) {
      std::vector<geo::WireID> wireIDs = geom->ChannelToWire(chan);
      if(wireIDs.empty()) {
        continue;
      }
      geo::WireGeo const& wire = geom->Wire(wireIDs.front());
      double wirelength = wire.HalfL() * 2;
      fLengthHist->Fill(chan, wirelength);
      fPlaneHist->Fill(chan, geom->View(wireIDs.front().planeID()));
    }
    if(!(fLengthHist->GetEntries()>0.)) {
      fLengthHist->Fill(0.,0.); // so that GetEntries() >0 on later events, even if there are no wires.
    }
  }
  */
  art::Handle< std::vector<raw::RawDigit> > rawdigit;
  if (! e.getByLabel(fDigitModuleLabel, rawdigit)) {
    std::cout << "WARNING: no label " << fDigitModuleLabel << std::endl;
    return;
  }
  std::vector< art::Ptr<raw::RawDigit> >  wires;
  art::fill_ptr_vector(wires, rawdigit);

  for (auto const& wire: wires) {
    std::vector<Short_t> uncompressed(wire->Samples());
    raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());
    if(fMaxChannel < wire->Channel()) fMaxChannel=wire->Channel();
    if(wire->Channel() >= num_channels) {
      std::cerr << "Wire channel " << wire->Channel() << " > num_channels! " << std::endl;
      continue;
    }
    for (size_t j = 0; j < wire->Samples(); ++j){
      fTotalTicks[wire->Channel()] += 1;
      fTotalADC[wire->Channel()] += uncompressed[j];
      fTotalADCSquared[wire->Channel()] += uncompressed[j]*uncompressed[j];
      //fADCs[wire->Channel()].push_back(uncompressed[j]);
      //if(fChannelToHistogram >= 0 && wire->Channel() == fChannelToHistogram) {
      int ii = 0;
      for(auto i : fChannelsToHistogram) {
        if(wire->Channel()==i) {
        fAdcHists[ii]->Fill(uncompressed[j]);
        }
        ii++;
      }
    }
  }
}

void dqm::DQMChannelNoiseModule::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  //fAdcHist = tfs->make<TH2D>("hADC", "adc ", num_channels ,0,num_channels, 65536, 0, 65536);
  fMeanHist = tfs->make<TH1D>("hADCMean","adc mean", num_channels, 0, num_channels);
  fVarHist = tfs->make<TH1D>("hADCStdDev","adc std dev", num_channels, 0, num_channels);
  /*
  fLengthHist = tfs->make<TH1D>("hLength","wire length",num_channels,0,num_channels);
  fPlaneHist = tfs->make<TH1D>("hPlane","wire plane",num_channels,0,num_channels);
  */
  for(auto i : fChannelsToHistogram) {
    fAdcHists.push_back(tfs->make<TH1D>(Form("hADC%d",i),Form("adc chan %d",i),65536,0,65536));
  }
}

void dqm::DQMChannelNoiseModule::endJob()
{
  // Implementation of optional member function here.

  /*
  unsigned long long maxTotalTicks = 0;
  uint32_t maxChan = -1;
  for(auto i: fTotalTicks) {
    if(i.second > maxTotalTicks) {
      maxTotalTicks = i.second;
      maxChan = i.first;
    }
  }
  std::cout << "Max channel: " << maxChan << std::endl;
  */
  std::cout << "Max Channels no: " << fMaxChannel << std::endl;

  for(uint32_t i = 0; i < num_channels; ++i) {
    if(fTotalTicks[i]>0) {
      double N = fTotalTicks[i];
      double mean = fTotalADC[i]/N;
      fMeanHist->Fill(i,mean);
      double var = (fTotalADCSquared[i]-fTotalADC[i]*fTotalADC[i]/N)/N;
      fVarHist->Fill(i,sqrt(var));
    }
    else {
      fMeanHist->Fill(i,-1.);
      fVarHist->Fill(i,-1.);
    }
  }
  for(auto h : fAdcHists) {
    h->GetXaxis()->SetRange(h->FindFirstBinAbove(0),h->FindLastBinAbove(0));
  }
}

void dqm::DQMChannelNoiseModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fDigitModuleLabel = p.get<std::string>("DigitModuleLabel");
  fChannelsToHistogram = p.get<std::vector<uint32_t> >("ChannelsToHistogram");
}

DEFINE_ART_MODULE(dqm::DQMChannelNoiseModule)
