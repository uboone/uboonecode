//
// Name:  SignalShapingMicroBooNETest.h
//
// Purpose: SignalShapingMicroBooNETest module.  Test convolution/deconvolution.
//
// Created:  28-Nov-2011  H. Greenlee

#include <vector>
#include <iostream>
#include <cmath>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardata/Utilities/LArFFT.h"

#include "TComplex.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// Local functions.

namespace {

  // Fill histogram from vector (set underflow/overflow bins to zero).

  void vector_to_hist(const std::vector<float>& v, TH1D* h)
  {
    assert(h != 0);
    int nvec = v.size();
    int nbins = h->GetNbinsX();
    int nfill = std::min(nvec, nbins);
    h->SetBinContent(0, 0.);
    int zerobin = h->GetXaxis()->FindBin(0.);
    for(int i=1; i<=nfill; ++i) {
      if(i >= zerobin)
	h->SetBinContent(i, v[i - zerobin]);
      else
	h->SetBinContent(i, v[i - zerobin + nvec]);
    }
    for(int i=nfill+1; i<=nbins+1; ++i)
      h->SetBinContent(i, 0.);
  }

  // Fill vector with initial delta-function at bin d.

  void fill_delta(std::vector<float>& v, int d)
  {
    int n = v.size();
    assert(d >= 0 && d < n);
    for(int i=0; i<n; ++i)
      v[i] = 0.;
    v[d] = 1.;
  }
}

namespace util
{
  class SignalShapingMicroBooNETest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit SignalShapingMicroBooNETest(fhicl::ParameterSet const& pset);
    virtual ~SignalShapingMicroBooNETest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  };

  DEFINE_ART_MODULE(SignalShapingMicroBooNETest)

  SignalShapingMicroBooNETest::SignalShapingMicroBooNETest(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset)
  {}

  void SignalShapingMicroBooNETest::beginJob()
  {
    // Get services.

    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<util::LArFFT> fft;
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;

    int nticks = fft->FFTSize();
    std::cout << "Number of ticks = " << nticks << std::endl;

    // Make collection plane histograms.

    art::TFileDirectory dirc = tfs->mkdir("Collection", "Collection");
    int nhist = std::min(300, nticks);
    int nkern = nticks/2;

    // Make input pulse.

    std::vector<float> tinc(nticks, 0.);
    fill_delta(tinc, nhist/2);
    TH1D* hinc = dirc.make<TH1D>("input", "Collection Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tinc, hinc);

    // Convoluted pulse.

    std::vector<float> tconvc(tinc);
    sss->Convolute(6000, tconvc, 0.0, 0.0);
    TH1D* hconvc = dirc.make<TH1D>("conv", "Collection Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tconvc, hconvc);

    // Deconvoluted pulse.

    std::vector<float> tdeconvc(tconvc);
    sss->Deconvolute(6000, tdeconvc, 0.0, 0.0);
    TH1D* hdeconvc = dirc.make<TH1D>("deconv", "Collection Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tdeconvc, hdeconvc);

    // Get collection response function and fill histogram.

    const std::vector<float>&  respc = sss->SignalShaping(6000,"nominal").Response();
    TH1D* hrfc = dirc.make<TH1D>("resp", "Collection Response", nhist+1, -nhist/2-0.5, nhist/2+0.5);
    vector_to_hist(respc, hrfc);

    // Get collection convolution kernel and fill histogram.

    const std::vector<ComplexF>&  kernc = sss->SignalShaping(6000,"nominal").ConvKernel();
    std::vector<float> kernrc(kernc.size());
    for(unsigned int i=0; i<kernrc.size(); ++i)
      kernrc[i] = kernc[i].Rho();
    TH1D* hkernc = dirc.make<TH1D>("kern", "Collection Convolution Kernel", nkern+1, -0.5, nkern+0.5);
    hkernc->SetMinimum(0.);
    vector_to_hist(kernrc, hkernc);


    // Make induction plane histograms.

    art::TFileDirectory udiri = tfs->mkdir("InductionU", "InductionU");

    // Make input pulse.

    std::vector<float> utini(nticks, 0.);
    fill_delta(utini, nhist/2);
    TH1D* uhini = udiri.make<TH1D>("input", "InductionU Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(utini, uhini);

    // Convoluted pulse.

    std::vector<float> utconvi(utini);
    sss->Convolute(0, utconvi, 0.0, 0.0);
    TH1D* uhconvi = udiri.make<TH1D>("conv", "InductionU Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(utconvi, uhconvi);

    // Deconvoluted pulse.

    std::vector<float> utdeconvi(utconvi);
    sss->Deconvolute(0, utdeconvi, 0.0, 0.0);
    TH1D* uhdeconvi = udiri.make<TH1D>("deconv", "InductionU Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(utdeconvi, uhdeconvi);

    // Get induction response function and fill histogram.

    const std::vector<float>&  urespi = sss->SignalShaping(0,"nominal").Response();
    TH1D* uhrfi = udiri.make<TH1D>("resp", "InductionU Response", nhist+1, -nhist/2-0.5, nhist/2+0.5);
    vector_to_hist(urespi, uhrfi);

    // Get induction convolution kernel and fill histogram.

    const std::vector<ComplexF>&  ukerni = sss->SignalShaping(0,"nominal").ConvKernel();
    std::vector<float> ukernri(ukerni.size());
    for(unsigned int i=0; i<ukernri.size(); ++i)
      ukernri[i] = ukerni[i].Rho();
    TH1D* hukerni = udiri.make<TH1D>("kern", "InductionU Convolution Kernel", nkern+1, -0.5, nkern+0.5);
    hukerni->SetMinimum(0.);
    vector_to_hist(ukernri, hukerni);



    art::TFileDirectory diri = tfs->mkdir("InductionV", "InductionV");

    // Make input pulse.

    std::vector<float> tini(nticks, 0.);
    fill_delta(tini, nhist/2);
    TH1D* hini = diri.make<TH1D>("input", "InductionV Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tini, hini);

    // Convoluted pulse.

    std::vector<float> tconvi(tini);
    sss->Convolute(3000, tconvi, 0.0, 0.0);
    TH1D* hconvi = diri.make<TH1D>("conv", "InductionV Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tconvi, hconvi);

    // Deconvoluted pulse.

    std::vector<float> tdeconvi(tconvi);
    sss->Deconvolute(3000, tdeconvi, 0.0, 0.0);
    TH1D* hdeconvi = diri.make<TH1D>("deconv", "InductionV Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tdeconvi, hdeconvi);

    // Get induction response function and fill histogram.

    const std::vector<float>&  respi = sss->SignalShaping(3000,"nominal").Response();
    TH1D* hrfi = diri.make<TH1D>("resp", "InductionV Response", nhist+1, -nhist/2-0.5, nhist/2+0.5);
    vector_to_hist(respi, hrfi);

    // Get induction convolution kernel and fill histogram.

    const std::vector<ComplexF>&  kerni = sss->SignalShaping(3000,"nominal").ConvKernel();
    std::vector<float> kernri(kerni.size());
    for(unsigned int i=0; i<kernri.size(); ++i)
      kernri[i] = kerni[i].Rho();
    TH1D* hkerni = diri.make<TH1D>("kern", "InductionV Convolution Kernel", nkern+1, -0.5, nkern+0.5);
    hkerni->SetMinimum(0.);
    vector_to_hist(kernri, hkerni);

  }

  SignalShapingMicroBooNETest::~SignalShapingMicroBooNETest()
  {}

  void SignalShapingMicroBooNETest::analyze(const art::Event& evt)
 {}
}
