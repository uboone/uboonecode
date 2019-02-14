////////////////////////////////////////////////////////////////////////
// Class:       SNMichelAna
// Plugin Type: analyzer (art v2_11_02)
// File:        SNMichelAna_module.cc
//
// Generated at Thu Oct 18 16:31:00 2018 by Jose Ignacio Crespo Anadon using 
// cetskelgen -v -d . -e beginJob -e endJob analyzer SNMichelAna
// from cetlib version v3_03_01.
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

// Additional framework
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

//"larsoft" object includes
#include "lardataobj/RecoBase/Cluster.h"

// ROOT
#include "TH1F.h"

class SNMichelAna;


class SNMichelAna : public art::EDAnalyzer {
public:
  explicit SNMichelAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SNMichelAna(SNMichelAna const &) = delete;
  SNMichelAna(SNMichelAna &&) = delete;
  SNMichelAna & operator = (SNMichelAna const &) = delete;
  SNMichelAna & operator = (SNMichelAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TH1F* fheSpectrum;
  TH1F* fhgSpectrum;
  TH1F* fhtotSpectrum;

};


SNMichelAna::SNMichelAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
}

void SNMichelAna::beginJob()
{
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us.
  art::ServiceHandle<art::TFileService> tfs;

  // Implementation of optional member function here.
  fheSpectrum = tfs->make<TH1F>("heSpectrum","e cluster spectrum; MeV*; Entries/2 MeV", 40, 0, 80);
  fhgSpectrum = tfs->make<TH1F>("hgSpectrum","#gamma cluster spectrum; MeV*; Entries/2 MeV", 40, 0, 80);
  fhtotSpectrum = tfs->make<TH1F>("htotSpectrum","Total clusters spectrum; MeV*; Entries/2 MeV", 40, 0, 80);
}

void SNMichelAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  float fADC2MeV = 7.69e-3;
  int fSamplesOverlapPre = 1600;
  int fSamplesOverlapPost = 1600;
  int fTotalSamplesPerRecord = 6400;

  art::InputTag ecluster_tag { "michel", "electron", "MichelReco" };
  art::InputTag gcluster_tag { "michel", "photon", "MichelReco" };

  auto const& ecluster_handle = e.getValidHandle< std::vector<recob::Cluster> >(ecluster_tag);
  auto const& gcluster_handle = e.getValidHandle< std::vector<recob::Cluster> >(gcluster_tag);

  if(ecluster_handle->size() < 1) return;
  std::cout << "*** MICHEL ELECTRON EVENT ***" << std::endl;

  // Loop over clusters
  std::cout << "\tThere are " << gcluster_handle->size() << " gamma clusters in this event." << std::endl;;
  float radCharge = 0;
  // If there are several Michel e in the frame this is not correct
  for( auto const& cluster : *gcluster_handle){
    // Just count the Michel e in the central frame
    if( cluster.StartTick() < fSamplesOverlapPre || 
	cluster.StartTick() >= fTotalSamplesPerRecord - fSamplesOverlapPost ) continue;
    radCharge += cluster.Integral();
  }
  std::cout << "Gamma energies: " << radCharge << " ADC = " << fADC2MeV*radCharge << " MeV" << std::endl;
  fhgSpectrum->Fill( fADC2MeV*radCharge );

  std::cout << "\tThere are " << ecluster_handle->size() << " electron clusters in this event." << std::endl;;
  if( ecluster_handle->size() > 1 ) std::cout << "WARNING: POTENTIAL AMBIGUITY" << std::endl;

  for( auto const& cluster : *ecluster_handle){
    // Just count the Michel e in the central frame
    if( cluster.StartTick() < fSamplesOverlapPre || 
	cluster.StartTick() >= fTotalSamplesPerRecord - fSamplesOverlapPost ) continue;
    fheSpectrum->Fill( fADC2MeV*cluster.Integral() );
    std::cout << "Electron energy: " << cluster.Integral() << " ADC = " << fADC2MeV*cluster.Integral() << " MeV" << std::endl;
    fhtotSpectrum->Fill( fADC2MeV*radCharge + fADC2MeV*cluster.Integral() );
    std::cout << "Total energy: " << cluster.Integral() + radCharge << " ADC = " << fADC2MeV*(radCharge + cluster.Integral()) << " MeV" << std::endl;
  }

}

void SNMichelAna::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SNMichelAna)
