////////////////////////////////////////////////////////////////////////
//
// module to create a TTree for TPCnoise analysis
//
//
////////////////////////////////////////////////////////////////////////
#ifndef ANADATA_H
#define ANADATA_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "SimpleTypesAndConstants/RawTypes.h"                                                                                             
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RawData/BeamInfo.h"
#include "RawData/DAQHeader.h"
#include "RawData/OpDetWaveform.h"
#include "RawData/OpDetPulse.h"
#include "Geometry/Geometry.h"

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>

#include "TTree.h"
#include "TTimeStamp.h"
#include "TH2.h"
#include "TFile.h"


const int kMaxPrimaries = 20000; //maximum number of primary particles
namespace microboone {
  
  class Anadata : public art::EDAnalyzer {
    
  public:
    
    explicit Anadata(fhicl::ParameterSet const& pset); 
    virtual ~Anadata();
    
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);
    
  private:
    
    
    void   ResetVars();
    
    TTree* fTree;
    //run information
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;                
    Double_t evttime;              
    Double_t beamtime;             
    Char_t   isdata;//flag, 0=MC 1=data
    uint32_t channel[kMaxPrimaries];
    Double_t adc_rms[kMaxPrimaries];
    Double_t adc_rms_wire_plane_0[kMaxPrimaries];
    Double_t adc_rms_wire_plane_1[kMaxPrimaries];
    Double_t adc_rms_wire_plane_2[kMaxPrimaries];
    Double_t pedstal_wire_plane_0[kMaxPrimaries];
    Double_t pedstal_wire_plane_1[kMaxPrimaries];
    Double_t pedstal_wire_plane_2[kMaxPrimaries];
    Int_t    tot_channel;
    Int_t    tot_waveform;
    Double_t pedstal_forchannel[kMaxPrimaries];
    Int_t fDataSize;
    std::vector<double>holder;
    std::vector<double>diff_val_holder;
    Int_t event_n;
    std::vector<short> rawadc;
    std::vector<double> channels_pedestal; 
    
    TH2D *ped_histo_nevt;
    TH2D *rms_histo_nevt;
    Int_t N_evt_processed = 10;
    
    
    
    
    std::string processname[kMaxPrimaries];
    std::string fRawDigitModuleLabel;
    std::string fRawOpDetWaveModuleLabel;
    
  };
}//name space

microboone::Anadata::Anadata(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset),
  fRawDigitModuleLabel (pset.get<std::string>("RawDigitModuleLabel"))
{
}

//-------------------------------------------------
microboone::Anadata::~Anadata()
{
}

void microboone::Anadata::beginJob(){
  
  art::ServiceHandle<art::TFileService> tfs;
  fTree= tfs->make<TTree>("noiseonrawtree","noiseonraw");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("beamtime",&beamtime,"beamtime/D");
  fTree->Branch("isdata",&isdata,"isdata/B"); 
  fTree->Branch("tot_channel",&tot_channel,"tot_channel/I");
  fTree->Branch("pedstal_forchannel", pedstal_forchannel, "pedstal_forchannel[tot_channel]/D");
  fTree->Branch("channel", channel, "channel[tot_channel]/I");
  fTree->Branch("adc_rms", adc_rms, "adc_rms[tot_channel]/D");
  fTree->Branch("event_n", &event_n,  "event_n/I");
  fTree->Branch("adc_rms_wire_plane_0", adc_rms_wire_plane_0, "adc_rms_wire_plane_0[tot_channel]/D");
  fTree->Branch("adc_rms_wire_plane_1", adc_rms_wire_plane_1, "adc_rms_wire_plane_1[tot_channel]/D");
  fTree->Branch("adc_rms_wire_plane_2", adc_rms_wire_plane_2, "adc_rms_wire_plane_2[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_0", pedstal_wire_plane_0, "pedstal_wire_plane_0[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_1", pedstal_wire_plane_1, "pedstal_wire_plane_1[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_2", pedstal_wire_plane_2, "pedstal_wire_plane_2[tot_channel]/D");

  ped_histo_nevt = tfs->make<TH2D> ("ped_histo_nevt", "pedestal_forchannel_nevents", 8256, 0 ,8256, 100,0,2500);
  rms_histo_nevt = tfs->make<TH2D> ("rms_histo_nevt", "rms_forchannel_nevents", 8256, 0 ,8256, 100,0,50);
  
}

void microboone::Anadata::beginSubRun(const art::SubRun& sr)
{
}

void microboone::Anadata::analyze(const art::Event& evt)
{  
  ResetVars();
  //services
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();
  
  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beam",beam)){
    beamtime = (double)beam->get_t_ms(); /*millisecond*/
    beamtime/=1000.; //in second
  }

  art::ServiceHandle<geo::Geometry> fGeometry;
  unsigned int maxChannels  = fGeometry->Nchannels();
  unsigned int wireMaxNum[] = {fGeometry->Nwires(0),fGeometry->Nwires(1),fGeometry->Nwires(2)};
  /* 3456 Y wires arrayed vertically and the 2400 U and 2400 V*/
  std::cout<< " /maxChannels:"<< maxChannels << " /wireMaxNum[0]:"<< wireMaxNum[0]<< " /wireMaxNum[1]:"<< wireMaxNum[1]<< " /wireMaxNum[2]:"<< wireMaxNum[2]<<std::endl;
  
  
  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;
  
  
  
  art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
  evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
  tot_channel = rawDigitHandle->size();
  
  
  event_n = 0;
  for(size_t it = 0 ; it< rawDigitHandle->size(); ++it)
    {
      
      holder.clear();
      diff_val_holder.clear();
      rawadc.clear();
      
      art::Ptr<raw::RawDigit> digitVec(rawDigitHandle, it);
      fDataSize= digitVec->Samples();
      holder.resize(fDataSize);
      diff_val_holder.resize(fDataSize);
      rawadc.resize(fDataSize);
      
      double adc_pedestal = 0.;
      double adc_stdev    = 0.;
      double adc_diff_sq  = 0.;
      
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      
      for(int bin = 0; bin < fDataSize; ++bin){
	holder[bin] = rawadc[bin];   // Rawhitfinder_module.
	adc_pedestal += holder[bin]/fDataSize;
	
      }
      
      pedstal_forchannel[it] = adc_pedestal;
      channel[it]            = digitVec->Channel();
      
      ped_histo_nevt->Fill(digitVec->Channel(),adc_pedestal,1/(float)N_evt_processed);
      for (int bin = 0 ; bin < fDataSize; ++bin)
	{ 
	  diff_val_holder[bin] = (rawadc[bin]-adc_pedestal);
	  adc_diff_sq += (diff_val_holder[bin]*diff_val_holder[bin])/fDataSize;
	}
      adc_stdev = sqrt(adc_diff_sq);
      adc_rms[it] = adc_stdev;
      rms_histo_nevt->Fill(digitVec->Channel(),adc_stdev,1/(float)N_evt_processed);
      
      //get the wire plane... stuffs...
      std::vector<geo::WireID> wire_ids;
      wire_ids = fGeometry->ChannelToWire(digitVec->Channel());
      
      /*three wire planes*/ 
      /*std::cout<< "/wire_ID: "<< wire_ids[0].Plane<< " /wire:"<<wire_ids[0].Wire<<  "/pedstal_forchannel[it]:"<<pedstal_forchannel[it]<<std::endl;*/
      
      if(wire_ids[0].Plane==0) {adc_rms_wire_plane_0[it]=adc_stdev; pedstal_wire_plane_0[it]=adc_pedestal;}
      
      if(wire_ids[0].Plane==1) {adc_rms_wire_plane_1[it]=adc_stdev;  pedstal_wire_plane_1[it]=adc_pedestal;}
      
      if(wire_ids[0].Plane==2) {adc_rms_wire_plane_2[it]=adc_stdev;  pedstal_wire_plane_2[it]=adc_pedestal;}    
      
    }
  
  
  event_n++;
  fTree->Fill();
  
}

void microboone::Anadata::ResetVars(){
  
  run      = -99999;
  subrun   = -99999;
  event    = -99999;
  evttime  = -99999;
  beamtime = -99999;
  isdata   = -99;
  event_n = -999999;
  
  for (int i = 0; i<kMaxPrimaries; ++i)
    {
      pedstal_forchannel[i] = -999999;
      channel[i]            = -999999;
      adc_rms[i]            = -999999.;
      adc_rms_wire_plane_0[i]=-100.;
      adc_rms_wire_plane_1[i]=-100.;
      adc_rms_wire_plane_2[i]=-100.;
      pedstal_wire_plane_0[i]=-100.;
      pedstal_wire_plane_1[i]=-100.;
      pedstal_wire_plane_2[i]=-100.;
    }
}
namespace microboone{
  
  DEFINE_ART_MODULE(Anadata)
  
} 

#endif



