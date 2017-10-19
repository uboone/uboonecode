////////////////////////////////////////////////////////////////////////
// Class:       DQMTPCBurstModule
// Plugin Type: analyzer (art v2_07_03)
// File:        DQMTPCBurstModule_module.cc
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
  class DQMTPCBurstModule;
}

const uint32_t num_channels = 8256;


class dqm::DQMTPCBurstModule : public art::EDAnalyzer {
public:
  explicit DQMTPCBurstModule(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DQMTPCBurstModule(DQMTPCBurstModule const &) = delete;
  DQMTPCBurstModule(DQMTPCBurstModule &&) = delete;
  DQMTPCBurstModule & operator = (DQMTPCBurstModule const &) = delete;
  DQMTPCBurstModule & operator = (DQMTPCBurstModule &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  std::string fDigitModuleLabel;

  std::vector<int> fChannelToPlaneMap;
  TH1D *fUberUVYMaxHist;

  double calculate_median(const std::vector<double>& v);
};


dqm::DQMTPCBurstModule::DQMTPCBurstModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ this->reconfigure(p);
}

double dqm::DQMTPCBurstModule::calculate_median(const std::vector<double>& vec) {
  if(vec.empty()) return 0.;
  std::vector<double> vec2(vec);
  std::sort(vec2.begin(), vec2.end());
  if(vec2.size() % 2 == 1) return vec2.at(vec2.size()/2);
  else return 0.5*(vec2.at(vec2.size()/2) + vec2.at(vec2.size()/2+1));
}

void dqm::DQMTPCBurstModule::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  //
  if(fChannelToPlaneMap.empty()) {
    fChannelToPlaneMap.resize(num_channels);
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    for(uint32_t chan = 0; chan < num_channels; ++chan) {
      int plane = -1;
      {
        std::vector<geo::WireID> wireIDs = geom->ChannelToWire(chan);
        if(wireIDs.empty()) {
          continue;
        }
        plane = geom->View(wireIDs.front().planeID());
      }
      if(plane < 0 || plane > 2) continue;
      fChannelToPlaneMap[chan] = plane;
    }
  }
  //
  //
  art::Handle< std::vector<raw::RawDigit> > rawdigit;
  if (! e.getByLabel(fDigitModuleLabel, rawdigit)) {
    std::cout << "WARNING: no label " << fDigitModuleLabel << std::endl;
    return;
  }
  std::vector< art::Ptr<raw::RawDigit> >  wires;
  art::fill_ptr_vector(wires, rawdigit);


  unsigned int maxtick = 0;

  std::vector<std::vector<Short_t> > ucwires(num_channels);
  std::vector<double> uberU, uberV, uberY;

  for (auto const& wire: wires) {
    std::vector<Short_t> uncompressed(wire->Samples());
    if(wire->Samples() > maxtick) {
      maxtick = wire->Samples();
      uberU.resize(maxtick);
      uberV.resize(maxtick);
      uberY.resize(maxtick);
    }
    //std::vector<Short_t> uncompressed(wire->Samples());
    raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());
    int plane = fChannelToPlaneMap[wire->Channel()];
    for(size_t j = 0; j < maxtick; ++j) {
      if(plane == 0) {
        uberU[j] += uncompressed[j];
      }
      else if(plane == 1) {
        uberV[j] += uncompressed[j];
      }
      else if(plane == 2) {
        uberY[j] += uncompressed[j];
      }
    }
  }
  double medianU = calculate_median(uberU);
  double medianV = calculate_median(uberV);
  double medianY = calculate_median(uberY);

  for(size_t j = 0; j < maxtick; ++j) {
    uberU[j] = std::fabs(uberU[j] - medianU)/1e4;
    uberV[j] = std::fabs(uberV[j] - medianV)/1e4;
    uberY[j] = std::fabs(uberY[j] - medianY)/1e4;
  }

  double max_uvy500 = 0.;
  for(size_t i = 0; i < maxtick - 500; i++) {
    double uvy500 = 0.;
    double sumU = 0.;
    double sumV = 0.;
    double sumY = 0.;
    for(size_t j = i; j < maxtick && j < i+500; ++j) {
      //uvy500 += (uberU[j] * uberV[j] * uberY[j]) / 1e8;
      sumU += uberU[j];
      sumV += uberV[j];
      sumY += uberY[j];
    }
    uvy500 = (sumU * sumV * sumY) / 1e8;
    if(uvy500 > max_uvy500) {
      max_uvy500 = uvy500;
    }
  }

  fUberUVYMaxHist->Fill(max_uvy500);
}

void dqm::DQMTPCBurstModule::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  //fAdcHist = tfs->make<TH2D>("hADC", "adc ", num_channels ,0,num_channels, 65536, 0, 65536);
  //fUberYMaxHist = tfs->make<TH1D>("hUberYMax", "uber y-plane maximum", 100, 1, -1);
  //fUberUMaxHist = tfs->make<TH1D>("hUberUMax", "uber u-plane maximum", 100, 1, -1);
  //fUberVMaxHist = tfs->make<TH1D>("hUberVMax", "uber v-plane maximum", 100, 1, -1);
  fUberUVYMaxHist = tfs->make<TH1D>("hUberUVY500", "UVY500 metric", 100, 1, -1);
}

void dqm::DQMTPCBurstModule::endJob()
{
  // Implementation of optional member function here.

}

void dqm::DQMTPCBurstModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fDigitModuleLabel = p.get<std::string>("DigitModuleLabel");
}

DEFINE_ART_MODULE(dqm::DQMTPCBurstModule)
