////////////////////////////////////////////////////////////////////////
// Class:       BurstNoiseMetrics
// Module Type: producer
// File:        BurstNoiseMetrics_module.cc
//
// Generated at Tue Nov 28 11:16:47 2017 by Ivan Caro Terrazas using artmod
// from cetpkgsupport v1_13_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
//#include "BNMS.h"

// ROOT
#include "TVector3.h"
//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtuple.h"
//

#include <memory>
#include <iostream>
#include <utility>

class BurstNoiseMetrics;

class BurstNoiseMetrics : public art::EDProducer {
public:
  explicit BurstNoiseMetrics(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BurstNoiseMetrics(BurstNoiseMetrics const &) = delete;
  BurstNoiseMetrics(BurstNoiseMetrics &&) = delete;
  BurstNoiseMetrics & operator = (BurstNoiseMetrics const &) = delete;
  BurstNoiseMetrics & operator = (BurstNoiseMetrics &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:
  int median_func(std::vector<double> Chans); 
  // Declare member data here.
  std::string fFlashLabel;
  std::string fRawDigitLabel;

  int fWin_start;
  int fWin_end;
  int fIntWinFFTsum;
  int fIntWinNF;

};


BurstNoiseMetrics::BurstNoiseMetrics(fhicl::ParameterSet const & p)

// Initialize member data here.
{
 
// Data objects being created

// Call appropriate produces<>() functions here.

 produces< double >("CBmetric" ); //Cathode Burst Metric
 produces< int >   ("PMBmetric"); //Purity Monitor Burst Metric

//  produces< std::vector< bnms::MET > >();

  fFlashLabel        = p.get<std::string>("FlashLabel"          );
  fRawDigitLabel     = p.get<std::string>("RawDigitLabel"       );
  fWin_start         = p.get<int>        ("Win_start"           );
  fWin_end           = p.get<int>        ("Win_end"             );
  fIntWinFFTsum      = p.get<int>	 ("IntWinFFTsum"        );
  fIntWinNF	     = p.get<int>        ("IntWinNF"            );
}



void BurstNoiseMetrics::produce(art::Event & e)
{

//Initial parameters. Some have to be parametrized in the .fcl file 
  int numbins = fWin_end - fWin_start;
  double fftsum;
  int sum;
  int nf;
  std::vector<double> U_UberADCvals;
  double ADCval;
  int chanNum;
  double U_median;
 
  std::unique_ptr< double > fs2;
  std::unique_ptr< int > nf5;


//////////////////////////////////HISTOGRAMS////////////////////////////////////
/////////////////Setup histograms used for calculations/////////////////////////

///////////////Histogram storing the event flashes//////////////////////////////
  TH1F *h_time = new TH1F("h_time","Histogram;Counts;Time", numbins, fWin_start, fWin_end);

///////////////Histogram storing the U plane Uberwaveform/////////////////////// 
  TH1F *U_uberwf = new TH1F("U_uberwf", "Event UberWF for U plane;Time Tick; ADC Value",9594, -0.5, 9593.5);

//////////////Histogram used to store the FFT of U plane uber waveform//////////
  TH1  *FFT_U_uberwf = new TH1F("FFT_U_uberwf", "FFT of U_UberWF; ADCvals; Frequencies",9594, -0.5, 9593.5);


///////////////////////////////Locate object////////////////////////////////////
//Locate the flash objects within artroot file/////////////////////////////////
  auto const& opflash_handle = e.getValidHandle< std::vector < recob::OpFlash > > (fFlashLabel);
  auto const& opflash_vec(*opflash_handle);
 
//Locate the raw digit objects within the artroot file///////////////////////// 
  auto const& rawdigit_handle = e.getValidHandle< std::vector < raw::RawDigit > > (fRawDigitLabel);
  auto const& allrawdigits_vec(*rawdigit_handle);

//////////////////////////////nf500 Calculations///////////////////////////////
//Store any flash within the time window -1600 to 3200 ms
  for( auto const& flash : opflash_vec){
    if(flash.Time() < fWin_start || flash.Time() > fWin_end) continue;
    h_time -> Fill(flash.Time());      
  }
  std::vector<int> vec_getmax;
  for (int k = fIntWinNF - 1; k < numbins; k++ ){
     sum = h_time -> Integral(k - fIntWinNF - 1, k);
     vec_getmax.insert(vec_getmax.end(), sum);
  }
  nf = *std::max_element(vec_getmax.begin(), vec_getmax.end());
///////////////////////////////FFTSUM2 Calculations///////////////////////////
/////////////////////////////////////////////////
  int tot_tt = allrawdigits_vec.at(1).Samples();
  double *UChanADCval = new double[tot_tt];
  U_median = 0.;
  U_UberADCvals.clear();
  for(int k = 0; k < tot_tt; k++){	//UberWaveform Calculations
    for (size_t i_ar = 0, size_allrawdigits = rawdigit_handle->size(); i_ar != size_allrawdigits; ++i_ar){
        ADCval = allrawdigits_vec.at(i_ar).ADC(k);
        chanNum = allrawdigits_vec.at(i_ar).Channel();
        if(chanNum < 2400){   
          UChanADCval[k] += ADCval;
        }    
    } 
    U_UberADCvals.insert(U_UberADCvals.end(), UChanADCval[k]);              
  }
  U_median = median_func(U_UberADCvals);
  for(unsigned j = 0; j < U_UberADCvals.size(); j++){
    U_uberwf -> SetBinContent(j+1,(UChanADCval[j] - U_median)/pow(10.,3.));
  }
/////////////////////////////////UWF Calculations////////////////////////////////////////////
  FFT_U_uberwf = U_uberwf -> FFT(FFT_U_uberwf, "MAG");
  fftsum = FFT_U_uberwf -> Integral(0, fIntWinFFTsum); 

  delete h_time; 
  delete U_uberwf;
  delete FFT_U_uberwf;

  nf5 = std::unique_ptr< int >(new int (nf));
  fs2 = std::unique_ptr< double >(new double (fftsum));
   
  e.put(std::move(nf5),"PMBmetric"  );
  e.put(std::move(fs2),"CBmetric"   );

  // Implementation of required member function here.
}

///////////////////////////Functions////////////////////////////
int BurstNoiseMetrics::median_func(std::vector<double> Chans){
  size_t size = Chans.size();
  int median;
  sort(Chans.begin(), Chans.end());
  if (size  % 2 == 0){
     median = (Chans[size / 2 - 1] + Chans[size / 2]) / 2;
  }else{
    median = Chans[size / 2];
  }
  return median;
}
        

DEFINE_ART_MODULE(BurstNoiseMetrics)
