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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

//"larsoft" object includes
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

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

  std::string fHitProducer; // Module label that created the hits
  // std::string fMichelProcess; // Process name that created the Michel clusters
  std::string fMichelProducer; // Module label that created the Michel clusters

  float fADC2MeV; // Calibration factor
  int fSamplesOverlapPre; // Number of samples from event N - 1 in event N
  int fSamplesOverlapPost; // Number of samples from event N + 1 in event N
  int fTotalSamplesPerRecord; // Number of samples in event

  // Michel cluster's energy spectrum
  TH1F* fheSpectrum; // electron
  TH1F* fhgSpectrum; // photon(s)
  TH1F* fhtotSpectrum; // electron + photon(s)

  // Michel cluster's hit multiplicity
  TH1F* fheHitMult;
  TH1F* fhgHitMult;
  TH1F* fhtotHitMult;

  // Michel cluster's hit spectrum
  TH1F* fheHitSpectrum;
  TH1F* fhgHitSpectrum;
  TH1F* fhtotHitSpectrum;

  // Event-wide hit multiplicity
  TH1F* fhEventHitMult;

  // Event-wide hit spectrum
  TH1F* fhEventHitSpectrum;

  geo::View_t fPlane; // Wire plane analyzed
};


SNMichelAna::SNMichelAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  // fMichelProcess = p.get<std::string>("MichelProcess");
  fHitProducer = p.get<std::string>("HitProducer");
  fMichelProducer = p.get<std::string>("MichelProducer");
  fADC2MeV = p.get<float>("ADC2MeV");
  fSamplesOverlapPre = p.get<float>("SamplesOverlapPre", 1600.);
  fSamplesOverlapPost = p.get<float>("SamplesOverlapPost", 1600.);
  fTotalSamplesPerRecord = p.get<float>("TotalSamplesPerRecord", 6400.);

}

void SNMichelAna::beginJob()
{
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us.
  art::ServiceHandle<art::TFileService> tfs;

  // Implementation of optional member function here.
  fheSpectrum = tfs->make<TH1F>("heSpectrum","e cluster spectrum; Energy (MeV*); Entries/2 MeV", 40, 0, 80);
  fhgSpectrum = tfs->make<TH1F>("hgSpectrum","#gamma cluster spectrum; Energy (MeV*); Entries/2 MeV", 40, 0, 80);
  fhtotSpectrum = tfs->make<TH1F>("htotSpectrum","Total clusters spectrum; Energy (MeV*); Entries/2 MeV", 40, 0, 80);

  fheHitMult = tfs->make<TH1F>("heHitMult","e cluster hit mult.; Hit mult.; Entries", 100, 0, 100);
  fhgHitMult = tfs->make<TH1F>("hgHitMult","#gamma cluster hit mult.; Hit mult.; Entries", 100, 0, 100);
  fhtotHitMult = tfs->make<TH1F>("htotHitMult","Total clusters hit mult.; Hit mult.; Entries", 100, 0, 100);

  fheHitSpectrum = tfs->make<TH1F>("heHitSpectrum","e cluster hits spectrum; Hit integral (ADC); Entries/1 ADC", 1000, 0, 1000);
  fhgHitSpectrum = tfs->make<TH1F>("hgHitSpectrum","#gamma cluster hits spectrum; Hit integral (ADC); Entries/1 ADC", 1000, 0, 1000);
  fhtotHitSpectrum = tfs->make<TH1F>("htotHitSpectrum","Total clusters hits spectrum; Hit integral (ADC); Entries/1 ADC", 1000, 0, 1000); 
  // remove? it is fheHitSpectrum + fhgHitSpectrum

  fhEventHitMult = tfs->make<TH1F>("hEventHitMult","Event hit mult.; Hit mult.; Entries/200 hits", 100, 0, 2e4);

  fhEventHitSpectrum = tfs->make<TH1F>("hEventHitSpectrum","Event hits spectrum; Hit integral (ADC); Entries/1 ADC", 1000, 0, 1000);

  // Guess plane from Michel producer
  if( fMichelProducer.find("0") != std::string::npos) fPlane = (geo::View_t)0;
  else if( fMichelProducer.find("1") != std::string::npos) fPlane = (geo::View_t)1;
  else if( fMichelProducer.find("2") != std::string::npos) fPlane = (geo::View_t)2;
  else fPlane = geo::kUnknown; 
  std::cout << "MichelProducer is " << fMichelProducer << ". Guessed plane is " << fPlane << std::endl;
}

void SNMichelAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  ///// Event-wide section /////
  art::InputTag hit_tag { fHitProducer };

  auto const& hit_handle = e.getValidHandle< std::vector<recob::Hit> >(hit_tag);

  size_t evtHits = 0;
  for( auto const& hit : *hit_handle ){
    if( hit.View() != fPlane ) continue;
    // Avoid double counting
    if( hit.PeakTime() < fSamplesOverlapPre || 
	hit.PeakTime() >= fTotalSamplesPerRecord - fSamplesOverlapPost ) continue;
    fhEventHitSpectrum->Fill( hit.Integral() );
    evtHits++;
  }
  fhEventHitMult->Fill( (float)evtHits );

  ///// Michel section /////

  //art::InputTag ecluster_tag { fMichelProducer, "electron", fMichelProcess };
  art::InputTag ecluster_tag { fMichelProducer, "electron" };
  //art::InputTag gcluster_tag { fMichelProducer, "photon", fMichelProcess };
  art::InputTag gcluster_tag { fMichelProducer, "photon" };

  auto const& ecluster_handle = e.getValidHandle< std::vector<recob::Cluster> >(ecluster_tag);
  auto const& gcluster_handle = e.getValidHandle< std::vector<recob::Cluster> >(gcluster_tag);
  //auto const& ecluster_vec( *ecluster_handle );
  //auto const& gcluster_vec( *gcluster_handle );

  if(ecluster_handle->size() < 1) return;

  std::cout << "*** MICHEL ELECTRON EVENT ***" << std::endl;

  std::cout << "\tThere are " << ecluster_handle->size() << " electron clusters in this event." << std::endl;;
  // Reject events with more than one Michel electron since we cannot assign their photon clusters
  if( ecluster_handle->size() > 1 ){
    std::cout << "WARNING: POTENTIAL AMBIGUITY. SKIPPING EVENT." << std::endl;
    return;
    ///////////////////// TO DO: photon clusters can be matched to their corresponding electron clusters by the common start_wire, start_tick
    // Probability for > 1 Michel electron in the event is 0.0325% based on Michel Paper
  }
  std::cout << "\tThere are " << gcluster_handle->size() << " gamma clusters in this event." << std::endl;;

  // Grab hits associated with input clusters
  art::FindManyP<recob::Hit> ecluster_hit_assn_v( ecluster_handle, e, ecluster_tag );
  art::FindManyP<recob::Hit> gcluster_hit_assn_v( gcluster_handle, e, gcluster_tag );

  // Loop over clusters

  // Loop over photon clusters
  float radCharge = 0;
  size_t radHits = 0;
  // To do: if there are several Michel e in the frame, try to assign their photon clusters?
  //  for( auto const& cluster : *gcluster_handle ){
  for( size_t igc = 0; igc < gcluster_handle->size(); igc++ ){

    auto const& cluster = (*gcluster_handle)[igc];

    //auto const& cluster = gcluster_vec[igc];
    // Just count the Michel e in the central frame
    if( cluster.StartTick() < fSamplesOverlapPre || 
	cluster.StartTick() >= fTotalSamplesPerRecord - fSamplesOverlapPost ) continue;

    radCharge += cluster.Integral();

    const std::vector<art::Ptr<recob::Hit> >& ghit_v = gcluster_hit_assn_v.at(igc);
    radHits += ghit_v.size();
    // Loop over photon hits
    for( auto const& hit : ghit_v ){
      fhgHitSpectrum->Fill( hit->Integral() );
      fhtotHitSpectrum->Fill( hit->Integral() );
    }
  }
  std::cout << "Total gamma energy: " << radCharge << " ADC = " << fADC2MeV*radCharge << " MeV" << std::endl;
  // Fill histos only if there were in-time clusters
  if( radCharge != 0.0 ) fhgSpectrum->Fill( fADC2MeV*radCharge );
  if( radHits != 0 ) fhgHitMult->Fill( (float)radHits );

  // Loop over electron clusters
  // This loop is trivial since we are requiring only one electron cluster
  //for( auto const& cluster : *ecluster_handle ){ 
  for( size_t iec = 0; iec < ecluster_handle->size(); iec++ ){

    auto const& cluster = (*ecluster_handle)[iec];

    // Just count the Michel e in the central frame
    if( cluster.StartTick() < fSamplesOverlapPre || 
	cluster.StartTick() >= fTotalSamplesPerRecord - fSamplesOverlapPost ) continue;

    fheSpectrum->Fill( fADC2MeV*cluster.Integral() );
    std::cout << "Electron energy: " << cluster.Integral() << " ADC = " << fADC2MeV*cluster.Integral() << " MeV" << std::endl;
    fhtotSpectrum->Fill( fADC2MeV*radCharge + fADC2MeV*cluster.Integral() );
    std::cout << "Total energy: " << cluster.Integral() + radCharge << " ADC = " << fADC2MeV*(radCharge + cluster.Integral()) << " MeV" << std::endl;

    const std::vector<art::Ptr<recob::Hit> >& ehit_v = ecluster_hit_assn_v.at(iec);
    fheHitMult->Fill( (float)ehit_v.size() );
    fhtotHitMult->Fill( (float)(ehit_v.size() + radHits) );
    // Loop over electron hits
    for( auto const& hit : ehit_v ){
      fheHitSpectrum->Fill( hit->Integral() );
      fhtotHitSpectrum->Fill( hit->Integral() );
    }
  }

}

void SNMichelAna::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SNMichelAna)
