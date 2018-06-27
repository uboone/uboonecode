////////////////////////////////////////////////////////////////////////
// Class:       CRTCalibration
// Module Type: analyzer
// File:        CRTCalibration_module.cc
// Description: Module to obtain pedestals and gains from SiPMs (only pedestals so far)
// Generated at Fri Nov 25 16:03:09 2016 by dalorga using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CRTBernFEBDAQCore//Overlays/BernZMQFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "art/Framework/Services/Optional/TFileService.h"

#include "uboone/CRT/CRTAuxFunctions.hh"

#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TDatime.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>

namespace bernfebdaq {
  class CRTCalibration;
}

class bernfebdaq::CRTCalibration : public art::EDAnalyzer {
public:
  explicit CRTCalibration(fhicl::ParameterSet const & pset);
  
  // Required functions.
  void analyze(art::Event const & evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  //void reconfigure(fhicl::ParameterSet const & p) ;


private:

  // Declare member data here.

  //this is the magic stuff that sets up the tree.
  art::ServiceHandle<art::TFileService> tfs;
  
  std::string  raw_data_label_;
  int print_ped_;
  std::string  SiPMpositions_;
  
  TTree*       my_tree_;
  
  Double_t Hist_min;
  Double_t Hist_max;
  Int_t Hist_bin;

  std::map <int, TH1F*> histos; //key = FEB*100+ch
  TH1F* h0;
  TH1F* hpedpos;
  TH1F* hpedpos_err;
  TH1F* hpedsig;
  int key;

};


bernfebdaq::CRTCalibration::CRTCalibration(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset),
    raw_data_label_(pset.get<std::string>("raw_data_label")),
    print_ped_(pset.get<int>("print_ped")),
    SiPMpositions_(pset.get<std::string>("CRTpositions_file"))//,

{

}

void bernfebdaq::CRTCalibration::analyze(art::Event const & evt)
{

  //time_t evtTime = evt.time().value(); //for file name
  //std::cout<<"time: "<<evtTime<<std::endl;

  // look for raw BernFEB data  
  art::Handle< std::vector<artdaq::Fragment> > rawHandle;
  //evt.getByLabel(raw_data_label_, "BernZMQ", rawHandle); //artdaq //daq                    
  evt.getByLabel(raw_data_label_, rawHandle); //Converted files //crtdaq

  //check to make sure the data we asked for is valid
  if(!rawHandle.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " BernFEB fragments " << " in module " << raw_data_label_ << std::endl;
    std::cout << std::endl;
    return;
  }
  
  //get better access to the data
  std::vector<artdaq::Fragment> const& rawFragments(*rawHandle);

  //loop over the raw data fragments
  //There should be one fragment per FEB in each event.
  for(auto const& frag : rawFragments){
    
    //overlay this so it's in the "BernFragment" format. Same data!
    BernZMQFragment bfrag(frag);
    
    //Grab the metadata.
    //See bernfebdaq-core/bernfebdaq-core/Overlays/BernFEBFragment.hh
    auto bfrag_metadata = bfrag.metadata();
 
    size_t   nevents    = bfrag_metadata->n_events();   //number of BernFEBEvents in this packet

    auto FEB_ID =  bfrag_metadata->feb_id();     //mac addresss of this packet
    //auto IDN = crt::auxfunctions::getFEBN(FEB_ID); //FEB ID   
    auto IDN = FEB_ID;// for converted files  

    for(size_t i_e=0; i_e<nevents; ++i_e){
      BernZMQEvent const* this_event = bfrag.eventdata(i_e); //get the single hit/event
      
      for(size_t i_chan=0; i_chan<32; ++i_chan){
	
	if((this_event->adc[i_chan])>32 && (this_event->adc[i_chan])<4080){//remove fake data(just in case)
	  
	  key = IDN*100 + i_chan;  
	  
	  auto it = histos.find(key);
	  TH1F* h = (*it).second;
	  h->Fill(this_event->adc[i_chan]);
	  
	  h0->Fill(this_event->adc[i_chan]);
	  	  
	}
      }
    }
  }
}//end analyzer

void bernfebdaq::CRTCalibration::beginJob()
{
  // Implementation of optional member function here.
  
  my_tree_ = tfs->make<TTree>("my_tree","CRT Analysis Tree");
  
  Hist_min = 0;
  Hist_max = 4096;
  Hist_bin = 4096;
  

    
  char *namehist = new char[128];

  for (int FEB = 11; FEB<131; FEB++){
    
    if(FEB==13){FEB++;}
    if(FEB==25){FEB++;}
    if(FEB==62){FEB=105;}
    if(FEB==110){FEB++;}
    if(FEB==122){FEB++;}
    if(FEB==130){FEB=195;}

    for(int ch=0; ch<32; ch++){
      
      key = FEB*100 + ch;
      
      sprintf (namehist, "FEB_%d_SiPM_%d", FEB, ch);
      
      TH1F* h = tfs->make<TH1F>(namehist,namehist, Hist_bin, Hist_min, Hist_max);
      h->GetXaxis()->SetTitle("ADC Channel");
      h->GetYaxis()->SetTitle("Entries/bin");
      histos.emplace(key,h);
    }
  }
  
  h0 = tfs->make<TH1F>("h0","h0", Hist_bin, Hist_min, Hist_max);
  
}

void bernfebdaq::CRTCalibration::endJob(){
  
  int bin;
  hpedpos = tfs->make<TH1F>("hpedpos","CRT_pedestal_positions", 18600, 1000, 19600);
  hpedpos->GetXaxis()->SetTitle("Sensor ID (FEB*100+SiPMch)");
  hpedpos->GetYaxis()->SetTitle("ADC Counts");
  
  hpedpos_err = tfs->make<TH1F>("hpedpos_err","CRT_pedestal_position_error", 18600, 1000, 19600);
  hpedpos_err->GetXaxis()->SetTitle("Sensor ID (FEB*100+SiPM)");
  hpedpos_err->GetYaxis()->SetTitle("ADC Counts");

  hpedsig = tfs->make<TH1F>("hpedsig","CRT_pedestal_sigma", 18600, 1000, 19600);
  hpedsig->GetXaxis()->SetTitle("Sensor ID (FEB*100+SiPMch)");
  hpedsig->GetYaxis()->SetTitle("ADC Counts");

  std::map<int,std::pair<double,double> > values;  
  
  for (int FEB = 11; FEB<131; FEB++){

    if(FEB==13){FEB++;}
    if(FEB==25){FEB++;}
    if(FEB==62){FEB=105;}
    if(FEB==110){FEB++;}
    if(FEB==122){FEB++;}
    if(FEB==130){FEB=195;}

    for(int ch=0; ch<32; ch++){

      key = FEB*100 + ch;
      auto it = histos.find(key);
      TH1F* h = (*it).second;

      auto entries = h->GetEntries();
      
      if(entries>0){
	
	double maxat= h->GetMaximumBin();
	TF1 *fitped = new TF1("fitped","gaus", maxat-20, maxat+20);
	h->Fit(fitped, "R+");
	double ped = fitped->GetParameter(1);
	double ped_err= fitped->GetParError(1);
	double sig = fitped->GetParameter(2);

	std::pair<double,double> chped;
	chped = std::make_pair(ped, ped_err);
	values[key]=chped;
	
	bin = key - 1000;
	hpedpos->SetBinContent(bin,ped);
	hpedpos_err->SetBinContent(bin,ped_err);	
	hpedsig->SetBinError(bin, sig);

	//gain stuff
	//TH1F *hnew = (TH1F*)h->Clone("hnew");
	//TF1 *fitmean = new TF1("fitmean","gaus", maxat-20, maxat+20);
	//h->Fit(fitmean, "R+");
	//double mean = fitmean->GetParameter(1);
	//double mean_err= fitmean->GetParError(1);
	
      }
	

    }
  }

  hpedpos->SetMarkerStyle(2);
  hpedpos_err->SetMarkerStyle(2);
  
  if(print_ped_==1){  
    //create file   
    std::ofstream fs("CRTpedestals-V8.txt");
    char *pedname = new char[256];

    for(auto itA = begin(values); itA != end(values); ++itA){
      int key = (*itA).first;
      std::pair<double,double> val = (*itA).second;
      std::cout<<key<<"   "<<val.first<<"   "<<val.second<<std::endl;
      sprintf(pedname,"%d \t %g \t %g \n", key, val.first, val.second);  
      fs<<pedname;
    }
    fs.close();
  }
  
  
}


DEFINE_ART_MODULE(bernfebdaq::CRTCalibration)
