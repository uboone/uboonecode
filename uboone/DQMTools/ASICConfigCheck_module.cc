////////////////////////////////////////////////////////////////////////
// Class:       ASICConfigCheck
// Plugin Type: analyzer (art v2_05_00)
// File:        ASICConfigCheck_module.cc
//
// Generated at Sat Mar  4 14:14:08 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

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

#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TVirtualFFT.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

namespace dqm {
  class ASICConfigCheck;
}


class dqm::ASICConfigCheck : public art::EDAnalyzer {
public:
  explicit ASICConfigCheck(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ASICConfigCheck(ASICConfigCheck const &) = delete;
  ASICConfigCheck(ASICConfigCheck &&) = delete;
  ASICConfigCheck & operator = (ASICConfigCheck const &) = delete;
  ASICConfigCheck & operator = (ASICConfigCheck &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void beginRun(art::Run const & r) override;
  void endRun(art::Run const & r) override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  art::InputTag digit_tag;
  art::InputTag hit_tag;

  int inEventsFFTPerRun;
  size_t nEventsFFTPerRun;
    
  int hitMultiplicityCut;
  int hitRMSCut;

  float fftSignalRemovalCut;
  
  bool verbose;

  std::unordered_map<raw::ChannelID_t,size_t> nhits_map;  
  
  TVirtualFFT *fftr2c;
  std::vector<double> doubleVecInput;
  bool fft_Initialize;
  size_t n_evts_thisrun;

  TNtuple *nt_hits;
  TTree *nt_ch;

  typedef struct chTreeObj{
    unsigned int run;
    unsigned int subrun;
    unsigned int ev;
    unsigned int timestamp;
    unsigned int ch;
    float pedestal;
    float rms;
    float median;
    float rms_trunc;
    int n_hits;
    float fft_m_1q;
    float fft_m_2q;
    float fft_m_3q;
    float fft_m_4q;
    float fft_m_1q_trunc;
    float fft_m_2q_trunc;
    float fft_m_3q_trunc;
    float fft_m_4q_trunc; 
  } chTreeObj_t;
  
  chTreeObj_t chTreeData;

};


dqm::ASICConfigCheck::ASICConfigCheck(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fft_Initialize(true)
{ this->reconfigure(p); }

void dqm::ASICConfigCheck::reconfigure(fhicl::ParameterSet const & pset)
{
  digit_tag = pset.get<art::InputTag>("DigitModuleLabel");
  hit_tag = pset.get<art::InputTag>("HitModuleLabel");

  inEventsFFTPerRun = pset.get<int>("EventsFFTPerRun");
  nEventsFFTPerRun = (inEventsFFTPerRun>0) ? inEventsFFTPerRun : 999999;
    
  hitMultiplicityCut = pset.get<int>("HitMultiplicityCut");
  hitRMSCut = pset.get<int>("HitRMSCut");

  fftSignalRemovalCut = pset.get<float>("FFTSignalRemovalCut");
  if(fftSignalRemovalCut<0)
    fftSignalRemovalCut = 999999.;

  verbose = pset.get<bool>("Verbose",false);
}

void dqm::ASICConfigCheck::analyze(art::Event const & ev)
{

  nhits_map.clear();
  
  auto const& hit_handle = ev.getValidHandle< std::vector<recob::Hit> >(hit_tag);
  auto const& hitVec = *hit_handle;
  
  for(auto const& hit : hitVec){
    
    if(hit.Multiplicity() <= hitMultiplicityCut &&
       hit.RMS() < hitRMSCut ){
      nt_hits->Fill(ev.run(),
		    ev.subRun(),
		    ev.event(),
		    ev.time().timeHigh(),
		    hit.Channel(),
		    hit.PeakAmplitude(),
		    hit.Integral(),
		    hit.SummedADC(),
		    hit.RMS(),
		    hit.Multiplicity());
      nhits_map[hit.Channel()] += 1;
    }      
  }//end loop over hits

  chTreeData.run = ev.run();
  chTreeData.subrun = ev.subRun();
  chTreeData.ev = ev.event();
  chTreeData.timestamp = ev.time().timeHigh();
  chTreeData.fft_m_1q=-1;
  chTreeData.fft_m_2q=-1;
  chTreeData.fft_m_3q=-1;
  chTreeData.fft_m_4q=-1;
  chTreeData.fft_m_1q_trunc=-1;
  chTreeData.fft_m_2q_trunc=-1;
  chTreeData.fft_m_3q_trunc=-1;
  chTreeData.fft_m_4q_trunc=-1;

  auto const& digit_handle = ev.getValidHandle< std::vector<raw::RawDigit> >(digit_tag);
  auto const& digitVec = *digit_handle;
  for(auto const& digit : digitVec){
    
    chTreeData.ch = digit.Channel();
    doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());
    
    //sort to get median
    std::nth_element(doubleVecInput.begin(),doubleVecInput.begin()+(doubleVecInput.size()/2),doubleVecInput.end());
    chTreeData.median = *(doubleVecInput.begin()+(doubleVecInput.size()/2));
    
    //sort to get first and third quartile
    std::nth_element(doubleVecInput.begin(),doubleVecInput.begin()+(doubleVecInput.size()/4),doubleVecInput.begin()+(doubleVecInput.size()/2));
    std::nth_element(doubleVecInput.begin()+(doubleVecInput.size()/2),doubleVecInput.begin()+(3*doubleVecInput.size()/4),doubleVecInput.end());
    
    float q1 = *(doubleVecInput.begin()+(doubleVecInput.size()/4));
    float q3 = *(doubleVecInput.begin()+(3*doubleVecInput.size()/4));
    
    chTreeData.rms_trunc = std::sqrt(0.5*((q1-chTreeData.median)*(q1-chTreeData.median) +
					  (q3-chTreeData.median)*(q3-chTreeData.median))) / 0.6745;
    
    chTreeData.pedestal = (float)std::accumulate(digit.ADCs().begin(),digit.ADCs().end(),0) / (float)(digit.ADCs().size());
    chTreeData.rms = 0.0;
    std::for_each (digit.ADCs().begin(), digit.ADCs().end(), [&](const float d) {
	chTreeData.rms += (d - chTreeData.pedestal) * (d - chTreeData.pedestal);
      });
    chTreeData.rms = std::sqrt(chTreeData.rms / digit.ADCs().size());
    
    
    
    if(n_evts_thisrun < nEventsFFTPerRun){
      
      if(fft_Initialize){
	int ndim = digit.ADCs().size();
	fftr2c = TVirtualFFT::FFT(1, &ndim, "R2C EX K");
	fft_Initialize = false;
      }
      
      doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());
      fftr2c->SetPoints(doubleVecInput.data());
      fftr2c->Transform();
      
      TH1 *hfft_m = 0;
      hfft_m = TH1::TransformHisto(fftr2c,hfft_m,"MAG");
      
      chTreeData.fft_m_1q = hfft_m->Integral(2,hfft_m->GetNbinsX()/8);
      chTreeData.fft_m_2q = hfft_m->Integral(hfft_m->GetNbinsX()/8 + 1,2*hfft_m->GetNbinsX()/8);
      chTreeData.fft_m_3q = hfft_m->Integral(2*hfft_m->GetNbinsX()/8 + 1,3*hfft_m->GetNbinsX()/8);
      chTreeData.fft_m_4q = hfft_m->Integral(3*hfft_m->GetNbinsX()/8 + 1,4*hfft_m->GetNbinsX()/8);
      
      delete hfft_m;
      
      doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());
      for(auto & val : doubleVecInput)
	if( std::abs(val-chTreeData.median)>(fftSignalRemovalCut*chTreeData.rms_trunc) ) val = chTreeData.median;
      
      fftr2c->SetPoints(doubleVecInput.data());
      fftr2c->Transform();
      
      TH1 *hfft_m_trunc = 0;
      hfft_m_trunc = TH1::TransformHisto(fftr2c,hfft_m_trunc,"MAG");
      
      chTreeData.fft_m_1q_trunc = hfft_m_trunc->Integral(2,hfft_m_trunc->GetNbinsX()/8);
      chTreeData.fft_m_2q_trunc = hfft_m_trunc->Integral(hfft_m_trunc->GetNbinsX()/8 + 1,2*hfft_m_trunc->GetNbinsX()/8);
      chTreeData.fft_m_3q_trunc = hfft_m_trunc->Integral(2*hfft_m_trunc->GetNbinsX()/8 + 1,3*hfft_m_trunc->GetNbinsX()/8);
      chTreeData.fft_m_4q_trunc = hfft_m_trunc->Integral(3*hfft_m_trunc->GetNbinsX()/8 + 1,4*hfft_m_trunc->GetNbinsX()/8);
      
      delete hfft_m_trunc;
      
    }
    
    nt_ch->Fill();
  }//end loop over digit vector
  
  ++n_evts_thisrun;
  
}//end analyze function

void dqm::ASICConfigCheck::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  nt_hits = tfs->make<TNtuple>("nt_hits","Hits Ntuple","run:subrun:ev:timestamp:ch:hit_amp:hit_integral:hit_sumadc:hit_rms:hit_mult");
  nt_ch   = tfs->make<TTree>("nt_ch","Channel Ntuple");
  nt_ch->Branch("",&chTreeData,"run/i:subrun/i:ev/i:timestamp/i:ch/i:pedestal/F:rms/F:median/F:rms_trunc/F:n_hits/I:fft_m_1q/F:fft_m_2q/F:fft_m_3q/F:fft_m_4q/F:fft_m_1q_trunc/F:fft_m_2q_trunc/F:fft_m_3q_trunc/F:fft_m_4q_trunc/F");

}

void dqm::ASICConfigCheck::beginRun(art::Run const & r)
{
  n_evts_thisrun = 0;
}

void dqm::ASICConfigCheck::endRun(art::Run const & r){}
void dqm::ASICConfigCheck::endJob(){}

DEFINE_ART_MODULE(dqm::ASICConfigCheck)
