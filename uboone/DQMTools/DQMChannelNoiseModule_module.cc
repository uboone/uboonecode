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

#include "TH2F.h"

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

  std::string fDigitModuleLabel;

  TH2F *fVarHist;
};


dqm::DQMChannelNoiseModule::DQMChannelNoiseModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
 // More initializers here.
{ this->reconfigure(p);
}

void dqm::DQMChannelNoiseModule::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  //
  art::Handle< std::vector<raw::RawDigit> > rawdigit;
  if (! e.getByLabel(fDigitModuleLabel, rawdigit)) {
    std::cout << "WARNING: no label " << fDigitModuleLabel << std::endl;
    return;
  }
  std::vector< art::Ptr<raw::RawDigit> >  wires;
  art::fill_ptr_vector(wires, rawdigit);

  for (auto const& wire: wires) {
    std::map<Short_t, unsigned int> entries;
    std::vector<Short_t> uncompressed(wire->Samples());
    raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());
    for (size_t j = 0; j < wire->Samples(); ++j){
      entries[uncompressed[j]]+=1;
    }
    // RMS based on 16 - 84 % quantiles
    float num_seen = 0.;
    Short_t low16 = -1;
    Short_t upp16 = 0x7FFF;
    for(auto entry : entries) {
      if(upp16==0x7FFF && num_seen / (float)wire->Samples() >= 0.84) {
        upp16 = entry.first;
      }
      num_seen += entry.second;
      if(low16 < 0 && num_seen / (float) wire->Samples() >= 0.16) {
        low16 = entry.first;
      }
    }
    long long total = 0, totSq = 0;
    double N=0;
    for(auto entry : entries) {
      if(entry.first >= low16 && entry.first < upp16) {
        long long adc = entry.first;
        long long nent = entry.second;
        total += nent * adc;
        totSq += nent * adc * adc;
        N += nent;
      }
    }
    float wire_rms = ((N>0) ? ((totSq - total/N*total) / N) : 0.);
    fVarHist->Fill(wire->Channel(), std::log10(wire_rms+1e-6)); // 1e-6 to offset any log(0) errors

    /*
    if(!(std::log10(wire_rms+1e-6)>=-6.&&std::log10(wire_rms+1e-6)<=3.)) {
      for(auto entry : entries) {
        std::cout << "entries["<<entry.first<<"] = "<<entry.second<<std::endl;
      }
      std::cout << "Chan: " << wire->Channel() << std::endl;
      std::cout << "NS " << num_seen << " " << low16 << " " << upp16 << std::endl;
      std::cout << "final " << total << " sq: " << totSq << " n:" << N << " rms "<<wire_rms << std::endl ;
      throw std::exception();
    }
    */
  }
}

void dqm::DQMChannelNoiseModule::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // cuts will be based on these numbers
  const int nImportantBins = 7;
  double important_bins[nImportantBins] = { 0.9, 1., 1.3, 1.9, 5., 8., 11. };
  double important_bins_log[nImportantBins];
  for(int i = 0; i < nImportantBins; ++i) important_bins_log[i] = std::log10(important_bins[i]);
  const int nbins = 40;
  double binsYlow = -1.5;
  double binsYup = 3.;
  std::vector<double> binsY;
  for(int i = 0; i < nbins+1-nImportantBins; ++i) {
    binsY.push_back(binsYlow + i * (binsYup - binsYlow) / (double)(nbins-nImportantBins));
  }
  for(int i = 0; i < nImportantBins; ++i) {
    binsY.push_back(important_bins_log[i]);
  }
  std::sort(binsY.begin(), binsY.end());

  fVarHist = tfs->make<TH2F>("hADCRMS","adc rms;Channel number;Log_{10} RMS", num_channels, 0, num_channels,nbins,binsY.data());
}

void dqm::DQMChannelNoiseModule::endJob()
{
  // Implementation of optional member function here.

}

void dqm::DQMChannelNoiseModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fDigitModuleLabel = p.get<std::string>("DigitModuleLabel");
}


DEFINE_ART_MODULE(dqm::DQMChannelNoiseModule)
