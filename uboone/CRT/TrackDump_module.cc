////////////////////////////////////////////////////////////////////////
// Class:       TrackDump
// Module Type: analyzer
// File:        TrackDump_module.cc
//
// Generated at Mon Jul  3 03:51:03 2017 by David Lorca Galindo using artmod
// from cetpkgsupport v1_11_00.
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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <artdaq-core/Data/Fragment.hh>

#include "art/Framework/Services/Optional/TFileService.h"

#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "uboone/CRT/CRTProducts/CRTTzero.hh"
#include "uboone/CRT/CRTProducts/CRTTrack.hh"
#include "uboone/CRT/CRTAuxFunctions.hh"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3S.h"
#include "TProfile.h"
#include "TF1.h"
#include "TDatime.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <typeinfo>

const int kMaxCRThits = 1000;
const int kMaxCRTtzeros = 1000;
const int kMaxCRTtracks = 1000;
const int kMaxTPCtracks = 100;
const int kMaxPMTflashes = 100;


 // namespace crt {
 //   class TrackDump;
 // }

class TrackDump : public art::EDAnalyzer {
public:
  explicit TrackDump(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackDump(TrackDump const &) = delete;
  TrackDump(TrackDump &&) = delete;
  TrackDump & operator = (TrackDump const &) = delete;
  TrackDump & operator = (TrackDump &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void ResetVars();

  art::ServiceHandle<art::TFileService> tfs;
  // Declare member data here.
  
  uint32_t fEvtNum; //Number of current event                       
  uint32_t frunNum;                //Run Number taken from event  
  uint32_t fsubRunNum;             //Subrun Number taken from event         
  std::string  fTrackModuleLabel;
  bool fSaveTPCTrackInfo;
  bool fSavePMTFlashInfo;
  std::string data_labeltrack_;
  std::string data_labeltzero_;
  std::string data_labelhit_;
  std::string data_label_flash_;
  std::string data_label_DAQHeader_;
  int fHardDelay_;
  int verbose_;
  

  //quality plots

  TH2F* hplavspla;
  TH1F* hTlength;
  TH1F* hTtime;
  TH2F* hTlengthvsTime;
  TH2F* hTlengthvsTimeAbs;
  TProfile* hTlengthvsTimeAbs_prof;
  TH1F* htheta;
  TH1F* hphi;
  TH1F* hts0_ns;
  TH2F* hTvsH;

  TH2F* HitDistBot;
  TH2F* HitDistFT;
  TH2F* HitDistPipe;
  TH2F* HitDistTop;
 //quality plots                                                                                                                                          

  TTree*       fTree;
  // run information
  int run;
  int subrun;
  int event;
  double evttime_sec;
  double evttime_nsec;
  double evttime_GPS_sec;
  double evttime_GPS_nsec;
  // CRT hits
  int nCRThits;
  int hit_plane[kMaxCRThits];
  double hit_time_s[kMaxCRThits];
  double hit_time0[kMaxCRThits];
  double hit_time1[kMaxCRThits];
  double hit_charge[kMaxCRThits];
  double hit_posx[kMaxCRThits];
  double hit_posy[kMaxCRThits];
  double hit_posz[kMaxCRThits]; 
  // CRT tzeros
  int nCRTtzeros;
  double tz_time_s[kMaxCRTtzeros];
  double tz_time0[kMaxCRTtzeros];
  double tz_time1[kMaxCRTtzeros];
  // CRT tracks
  int nCRTtracks;
  double ct_theta[kMaxCRTtracks];
  double ct_phi[kMaxCRTtracks];
  double ct_length[kMaxCRTtracks];
  double ct_time_sec[kMaxCRTtracks];
  double ct_time0[kMaxCRTtracks];
  double ct_time1[kMaxCRTtracks];
  double ct_x1[kMaxCRTtracks];
  double ct_y1[kMaxCRTtracks];
  double ct_z1[kMaxCRTtracks];
  double ct_x2[kMaxCRTtracks];
  double ct_y2[kMaxCRTtracks];
  double ct_z2[kMaxCRTtracks];
  // TPC tracks
  int nTPCtracks;
  double trkstartx[kMaxTPCtracks];
  double trkstarty[kMaxTPCtracks];
  double trkstartz[kMaxTPCtracks];
  double trkendx[kMaxTPCtracks];
  double trkendy[kMaxTPCtracks];
  double trkendz[kMaxTPCtracks];
  double trkstartdcosx[kMaxTPCtracks];
  double trkstartdcosy[kMaxTPCtracks];
  double trkstartdcosz[kMaxTPCtracks];
  double trkenddcosx[kMaxTPCtracks];
  double trkenddcosy[kMaxTPCtracks];
  double trkenddcosz[kMaxTPCtracks];
  double trktheta[kMaxTPCtracks];
  double trkphi[kMaxTPCtracks];
  double trklen[kMaxTPCtracks];
  //  Flash information 
  int nPMTflash;
  double Yflash[kMaxPMTflashes];
  double Zflash[kMaxPMTflashes];
  double PEflash[kMaxPMTflashes];
  double Timeflash[kMaxPMTflashes];
  double fbeam[kMaxPMTflashes];

};


TrackDump::TrackDump(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
    fSaveTPCTrackInfo(p.get< bool >("SaveTPCTrackInfo", false)), 
    fSavePMTFlashInfo(p.get< bool >("SavePMTFlashInfo", false)), 
    data_labeltrack_(p.get<std::string>("data_labeltrack")),
    data_labeltzero_(p.get<std::string>("data_labeltzero")),
    data_labelhit_(p.get<std::string>("data_labelhit")),
    data_label_flash_(p.get<std::string>("data_label_flash_")),
    data_label_DAQHeader_(p.get<std::string>("data_label_DAQHeader_")),
    fHardDelay_(p.get<int>("fHardDelay",40000)), // 40 us
    verbose_(p.get<int>("verbose"))
    // More initializers here.    
{
}

void TrackDump::analyze(art::Event const & evt)
{
  
  ResetVars();

  frunNum    = evt.run();
  fsubRunNum = evt.subRun();
  fEvtNum = evt.event();
  
  art::Timestamp evtTime = evt.time();
  auto evt_time_sec = evtTime.timeHigh();
  auto evt_time_nsec = evtTime.timeLow();

   //get DAQ Header                                                                  
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;  
  evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
  

  //check to make sure the data we asked for is valid 
  if(!rawHandle_DAQHeader.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
	      << ", event " << evt.event() << " has zero"
	      << " DAQHeaderTimeUBooNE  " << " in with label " << data_label_DAQHeader_ << std::endl;    
    return;
  }
  
  
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
  double evt_timeGPS_sec = evtTimeGPS.timeHigh();
  double evt_timeGPS_nsec = (double)evtTimeGPS.timeLow();
  art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
  double evt_timeNTP_sec = evtTimeNTP.timeHigh();
  double evt_timeNTP_nsec = (double)evtTimeNTP.timeLow();
  double timstp_diff = std::abs(evt_timeGPS_nsec - evt_timeNTP_nsec);
  //fill tree variables
  evttime_sec=evt_time_sec;
  evttime_nsec=evt_time_nsec;
  evttime_GPS_sec=evt_timeGPS_sec;
  evttime_GPS_nsec=evt_timeGPS_nsec;

  if(verbose_==1){
    std::cout<< "Run:  "<<frunNum << "   subRun: " <<fsubRunNum<<std::endl;
    std::cout<<"event: "<<fEvtNum <<std::endl;
    std::cout.precision(19);
    std::cout<<"  GPS time second:  "<<evt_timeGPS_sec<<std::endl;
    std::cout<<"  GPS time nano_second:  "<<evt_timeGPS_nsec<<std::endl;
    std::cout<<"  NTP time second:  "<<evt_timeNTP_sec<<std::endl;    
    std::cout<<"  NTP time nano_second:  "<<evt_timeNTP_nsec<<std::endl;
    std::cout<<"  event time second:  "<<evt_time_sec<<std::endl;
    std::cout<<"  event time nano_second:  "<<evt_time_nsec<<std::endl;
    std::cout<<"  difference between GPS and NTP:  "<<evt_timeGPS_nsec - evt_timeNTP_nsec<<" ns"<<std::endl;
    std::cout<<"  ABS difference between GPS and NTP:  "<<timstp_diff<<" ns"<<std::endl;
    
    if( (evt_time_sec==evt_timeGPS_sec) && (evt_time_nsec==evt_timeGPS_nsec))  std::cout<<" Event time type is: GPS  "<<std::endl;
    if( (evt_time_sec==evt_timeNTP_sec) && (evt_time_nsec==evt_timeNTP_nsec))  std::cout<<" Event time type is: NTP  "<<std::endl;
    //getchar();
  }  
    

  if (fSavePMTFlashInfo) {

    //get Optical Flash
    art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
    evt.getByLabel(data_label_flash_, rawHandle_OpFlash);
    
    std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
    
    if(verbose_==1){ 
      std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl; 
    }
    
    nPMTflash=OpFlashCollection.size();
    if (nPMTflash>kMaxPMTflashes) nPMTflash=kMaxPMTflashes;
    for(int i = 0; i< nPMTflash; i++) {
      
      recob::OpFlash my_OpFlash = OpFlashCollection[i];      
      Yflash[i]=my_OpFlash.YCenter();
      Zflash[i]=my_OpFlash.ZCenter();
      PEflash[i]=my_OpFlash.TotalPE();
      Timeflash[i]=my_OpFlash.Time(); //in us from trigger time
      fbeam[i]=my_OpFlash.OnBeamTime();
    }
    
  }

  if (fSaveTPCTrackInfo) {

    // get TPC Track List 
    art::Handle< std::vector<recob::Track>  > trackListHandle; 
    std::vector<art::Ptr<recob::Track> >  tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);
    //check to make sure the data we asked for is valid
    if(!trackListHandle.isValid()){
      std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
		<< ", event " << evt.event() << " has zero"
		<< " tracks " << " in module " << fTrackModuleLabel << std::endl;
      std::cout << std::endl;
      return;
    }
    
    
    nTPCtracks = tracklist.size();
    if (nTPCtracks>kMaxTPCtracks) nTPCtracks=kMaxTPCtracks;
    for(int j = 0; j < nTPCtracks; j++) {
      
      TVector3 pos, dir_start, dir_end, end;              
      art::Ptr<recob::Track> ptrack(trackListHandle, j);
      const recob::Track& track = *ptrack;
      
      pos       = track.Vertex();
      dir_start = track.VertexDirection();
      dir_end   = track.EndDirection();
      end       = track.End();
      //
      trklen[j]= track.Length(); //length(track);
      trkstartx[j]=pos.X();
      trkstarty[j]=pos.Y();
      trkstartz[j]=pos.Z();
      trkendx[j]=end.X();
      trkendy[j]=end.Y();
      trkendz[j]=end.Z();
      trkstartdcosx[j]=dir_start.X();
      trkstartdcosy[j]=dir_start.Y();
      trkstartdcosz[j]=dir_start.Z();
      trkenddcosx[j]=dir_end.X();
      trkenddcosy[j]=dir_end.Y();
      trkenddcosz[j]=dir_end.Z();
      trktheta[j]=dir_start.Theta();
      trkphi[j]=dir_start.Phi();
      
    }
  }   //  if (saveTPCtrackinfo)

  /*
  
  //get Optical Flash
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlash;
  evt.getByLabel(data_label_flash_, rawHandle_OpFlash);
  
  std::vector<recob::OpFlash> const& OpFlashCollection(*rawHandle_OpFlash);
  
  if(verbose_==1){ 
    std::cout<<"  OpFlashCollection.size()  "<<OpFlashCollection.size()<<std::endl; 
  }  //get Optical Flash
  
  */


  //  fill tree
  run=frunNum;
  event=fEvtNum;
  subrun=fsubRunNum;
  
  //get CRTHits
  art::Handle< std::vector<crt::CRTHit> > rawHandle_hit;
  evt.getByLabel(data_labelhit_, rawHandle_hit); //
  
  //check to make sure the data we asked for is valid
  if(!rawHandle_hit.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_labelhit_ << std::endl;
    std::cout << std::endl;
    return;
  }
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hit);
  if(verbose_==1){ 
    std::cout<<"  CRTHitCollection.size()  "<<CRTHitCollection.size()<<std::endl; 
    //  getchar();   
  }    //end get CRTHits


  nCRThits = CRTHitCollection.size();
  if (nCRThits>kMaxCRThits) nCRThits=kMaxCRThits;
  for(int j = 0; j < nCRThits; j++) {
    
    //fill tree
    crt::CRTHit my_CRTHit = CRTHitCollection[j];
    hit_time_s[j]=(double)my_CRTHit.ts0_s;
    hit_time0[j]=(double)my_CRTHit.ts0_ns - (double)evt_timeGPS_nsec;
    hit_time1[j]=(double)my_CRTHit.ts1_ns + (double)fHardDelay_;  
    hit_charge[j]=my_CRTHit.peshit;
    hit_plane[j]=my_CRTHit.plane;
    hit_posx[j]=my_CRTHit.x_pos;
    hit_posy[j]=my_CRTHit.y_pos;
    hit_posz[j]=my_CRTHit.z_pos;

    //fillhistograms
    if (my_CRTHit.plane==0) HitDistBot->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos);
    else if (my_CRTHit.plane==1) HitDistFT->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    else if (my_CRTHit.plane==2) HitDistPipe->Fill(my_CRTHit.z_pos,my_CRTHit.y_pos);
    else if (my_CRTHit.plane==3) HitDistTop->Fill(my_CRTHit.z_pos,my_CRTHit.x_pos);



  }

  //get CRTTracks
  art::Handle< std::vector<crt::CRTTrack> > rawHandle_track;
  evt.getByLabel(data_labeltrack_, rawHandle_track); 
  
  //check to make sure the data we asked for is valid
  if(!rawHandle_track.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTTracks " << " in module " << data_labeltrack_ << std::endl;
    std::cout << std::endl;
    return;
  }

  std::vector<crt::CRTTrack> const& CRTTrackCollection(*rawHandle_track);
  if(verbose_==1){ 
    std::cout<<"  CRTTrackCollection.size()  "<<CRTTrackCollection.size()<<std::endl; 
    getchar();   
  }


  nCRTtracks = CRTTrackCollection.size();
  if (nCRTtracks>kMaxCRTtracks) nCRTtracks=kMaxCRTtracks;

  for(int j = 0; j <nCRTtracks; j++) {
    crt::CRTTrack my_CRTTrack = CRTTrackCollection[j];
    double temp = (my_CRTTrack.y1_pos-my_CRTTrack.y2_pos)*(my_CRTTrack.y1_pos-my_CRTTrack.y2_pos)
      +(my_CRTTrack.x1_pos-my_CRTTrack.x2_pos)*(my_CRTTrack.x1_pos-my_CRTTrack.x2_pos);
    double thetatemp =  atan2(sqrt(temp),my_CRTTrack.z1_pos-my_CRTTrack.z2_pos);
    double phitemp = atan2(my_CRTTrack.y1_pos-my_CRTTrack.y2_pos,my_CRTTrack.x1_pos-my_CRTTrack.x2_pos);
    // flip track to point downwards if needed
    if (phitemp<0) { ct_phi[j]=phitemp+3.14159; ct_theta[j]=3.14159-thetatemp;}
    else { ct_theta[j]=thetatemp;     ct_phi[j]=phitemp;}
    ct_length[j]=my_CRTTrack.length;
    ct_time_sec[j]=(double)my_CRTTrack.ts0_s;
    // ct_time0[j]=(double)my_CRTTrack.ts0_ns;
    // ct_time1[j]=(double)my_CRTTrack.ts1_ns;
    ct_time0[j]=(double)my_CRTTrack.ts0_ns - (double)evt_timeGPS_nsec;
    ct_time1[j]=(double)my_CRTTrack.ts1_ns + (double)fHardDelay_;
    // std::cout << my_CRTTrack.ts0_ns << " " << ct_time0[j] << std::endl;
    // std::cout << my_CRTTrack.ts1_ns << " " << ct_time1[j] << std::endl;
    ct_x1[j]=my_CRTTrack.x1_pos;
    ct_y1[j]=my_CRTTrack.y1_pos;
    ct_z1[j]=my_CRTTrack.z1_pos;
    ct_x2[j]=my_CRTTrack.x2_pos;
    ct_y2[j]=my_CRTTrack.y2_pos;
    ct_z2[j]=my_CRTTrack.z2_pos;

    //fill histograms
    hplavspla->Fill(my_CRTTrack.plane1,my_CRTTrack.plane2);
    hTlength->Fill(my_CRTTrack.length);
    double time_diff = my_CRTTrack.ts0_ns_h1-my_CRTTrack.ts0_ns_h2;
    double time_diffABS = fabs(time_diff);
    hTtime->Fill(time_diffABS);
    hTlengthvsTimeAbs->Fill(my_CRTTrack.length,time_diffABS);
    hTlengthvsTimeAbs_prof->Fill(my_CRTTrack.length,time_diffABS);
    hTlengthvsTime->Fill(my_CRTTrack.length,time_diff);
    htheta->Fill(57.30*my_CRTTrack.thetaxy);
    if (my_CRTTrack.phizy>3.14159) 
      hphi->Fill(57.30*(my_CRTTrack.phizy-3.14159));
    else    hphi->Fill(57.30*my_CRTTrack.phizy);
    hts0_ns->Fill(my_CRTTrack.ts0_ns);


  }


  fTree->Fill();

  
}

void TrackDump::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("trackdump","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime_sec",&evttime_sec,"evttime_sec/D");
  fTree->Branch("evttime_nsec",&evttime_nsec,"evttime_nsec/D");
  fTree->Branch("evttime_GPS_sec",&evttime_GPS_sec,"evttime_GPS_sec/D");
  fTree->Branch("evttime_GPS_nsec",&evttime_GPS_nsec,"evttime_GPS_nsec/D");
  //
  fTree->Branch("nCRTtzeros",&nCRTtzeros,"nCRTtzeros/I");
  fTree->Branch("tz_time_s",tz_time_s,"tz_time_s[nCRTtzeros]/D");
  fTree->Branch("tz_time0",tz_time0,"tz_time0[nCRTtzeros]/D");
  fTree->Branch("tz_time1",tz_time1,"tz_time1[nCRTtzeros]/D");
  //
  fTree->Branch("nCRThits",&nCRThits,"nCRThits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nCRThits]/I");
  fTree->Branch("hit_time_s",hit_time_s,"hit_time_s[nCRThits]/D");
  fTree->Branch("hit_time0",hit_time0,"hit_time0[nCRThits]/D");
  fTree->Branch("hit_time1",hit_time1,"hit_time1[nCRThits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nCRThits]/D");
  fTree->Branch("hit_posx",hit_posx,"hit_posx[nCRThits]/D");
  fTree->Branch("hit_posy",hit_posy,"hit_posy[nCRThits]/D");
  fTree->Branch("hit_posz",hit_posz,"hit_posz[nCRThits]/D");
  // CRT tracks
  fTree->Branch("nCRTtracks",&nCRTtracks,"nCRTtracks/I");
  fTree->Branch("ct_theta",ct_theta,"ct_theta[nCRTtracks]/D");
  fTree->Branch("ct_phi",ct_phi,"ct_phi[nCRTtracks]/D");
  fTree->Branch("ct_length",ct_length,"ct_length[nCRTtracks]/D");
  fTree->Branch("ct_time_sec",ct_time_sec,"ct_time_sec[nCRTtracks]/D");
  fTree->Branch("ct_time0",ct_time0,"ct_time0[nCRTtracks]/D");
  fTree->Branch("ct_time1",ct_time1,"ct_time1[nCRTtracks]/D");
  fTree->Branch("ct_x1",ct_x1,"ct_x1[nCRTtracks]/D");
  fTree->Branch("ct_y1",ct_y1,"ct_y1[nCRTtracks]/D");
  fTree->Branch("ct_z1",ct_z1,"ct_z1[nCRTtracks]/D");
  fTree->Branch("ct_x2",ct_x2,"ct_x2[nCRTtracks]/D");
  fTree->Branch("ct_y2",ct_y2,"ct_y2[nCRTtracks]/D");
  fTree->Branch("ct_z2",ct_z2,"ct_z2[nCRTtracks]/D");
  //TPC tracks
  if (fSaveTPCTrackInfo) {
  fTree->Branch("nTPCtracks",&nTPCtracks,"nTPCtracks/I");
  fTree->Branch("trkstartx",trkstartx,"trkstartx[nTPCtracks]/D");
  fTree->Branch("trkstarty",trkstarty,"trkstarty[nTPCtracks]/D");
  fTree->Branch("trkstartz",trkstartz,"trkstartz[nTPCtracks]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[nTPCtracks]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[nTPCtracks]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[nTPCtracks]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[nTPCtracks]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[nTPCtracks]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[nTPCtracks]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[nTPCtracks]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[nTPCtracks]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[nTPCtracks]/D");
  fTree->Branch("trktheta",trktheta,"trktheta[nTPCtracks]/D");
  fTree->Branch("trkphi",trkphi,"trkphi[nTPCtracks]/D");
  fTree->Branch("trklen",trklen,"trklen[nTPCtracks]/D");
  }
  //PMT flashes
  if (fSavePMTFlashInfo) {
  fTree->Branch("nPMTflash",&nPMTflash,"nPMTflash/I");
  fTree->Branch("Yflash",Yflash,"Yflash[nPMTflash]/D");
  fTree->Branch("Zflash",Zflash,"Zflash[nPMTflash]/D");
  fTree->Branch("PEflash",PEflash,"PEflash[nPMTflash]/D");
  fTree->Branch("Timeflash",Timeflash,"Timeflash[nPMTflash]/D");
  fTree->Branch("fbeam",fbeam,"fbeam[nPMTflash]/D");
  }


  hplavspla = tfs->make<TH2F>("hplavspla","PlanevsPlane",4,0,4,4,0,4);
  hplavspla->GetXaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetYaxis()->SetTitle("Plane (0=Bottom, 1=FT, 2=Pipe, 3=Top)");
  hplavspla->GetZaxis()->SetTitle("Entries/bin");
  hplavspla->SetOption("COLZ");

  hTvsH = tfs->make<TH2F>("hTvsH","Track_vs_Hits",500,0,500,500,0,500);
  hTvsH->GetXaxis()->SetTitle("Number of CRTHits per event");
  hTvsH->GetYaxis()->SetTitle("Number of CRTTracks per event");
  hTvsH->GetZaxis()->SetTitle("Entries/bin");
  hTvsH->SetOption("COLZ");

  hTlength = tfs->make<TH1F>("hTlength","Track_Length",1500,0,1500);
  hTlength->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlength->GetYaxis()->SetTitle("Entries/bin");

  hTtime = tfs->make<TH1F>("hTtime","Track_time",120,-10,110);
  hTtime->GetXaxis()->SetTitle("Track_time (ns)");
  hTtime->GetYaxis()->SetTitle("Entries/bin");

  hTlengthvsTime = tfs->make<TH2F>("hTlengthvsTime","Track_LengthvsTime",1500,0,1500,200,-100,100);
  hTlengthvsTime->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTime->GetYaxis()->SetTitle("Track_time (ns)");
  hTlengthvsTime->GetZaxis()->SetTitle("Entries/bin");
  hTlengthvsTime->SetOption("COLZ");

  hTlengthvsTimeAbs = tfs->make<TH2F>("hTlengthvsTimeAbs","Track_LengthvsTimeAbs",1500,0,1500,110,-10,100);
  hTlengthvsTimeAbs->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTimeAbs->GetYaxis()->SetTitle("Track_time (ns)");
  hTlengthvsTimeAbs->GetZaxis()->SetTitle("Entries/bin");
  hTlengthvsTimeAbs->SetOption("COLZ");

  hTlengthvsTimeAbs_prof = tfs->make<TProfile>("hTlengthvsTimeAbs_prof","Track_LengthvsTimeAbs_prof",1500,0,1500,"s");
  hTlengthvsTimeAbs_prof->GetXaxis()->SetTitle("Track_Length (cm)");
  hTlengthvsTimeAbs_prof->GetYaxis()->SetTitle("Track_time (ns)");

  htheta = tfs->make<TH1F>("htheta","Track_theta",900,0,180);
  htheta->GetXaxis()->SetTitle("Theta_xy (º)");
  htheta->GetYaxis()->SetTitle("Entries/bin");
 
  hphi = tfs->make<TH1F>("hphi","Track_phi",900,0,180);
  hphi->GetXaxis()->SetTitle("Phi_zy (º)");
  hphi->GetYaxis()->SetTitle("Entries/bin");

  hts0_ns = tfs->make<TH1F>("hts0_ns","Track_time_ns",100000,0,1e9);
  hts0_ns->GetXaxis()->SetTitle("Track time (ns)");
  hts0_ns->GetYaxis()->SetTitle("Entries/bin");


  double inch =2.54; //inch in cm                                                                                                                          
  HitDistBot = tfs->make<TH2F>("hBottom","Bottom",125,-700+205*inch,-700+205*inch+125*10.89,60,-300+50.4*inch,-300+50.4*inch+60*10.89);
  HitDistBot->GetXaxis()->SetTitle("Length long the beam (cm)");
  HitDistBot->GetYaxis()->SetTitle("Length along the drift (cm)");
  HitDistBot->GetZaxis()->SetTitle("Entries/bin");
  HitDistBot->SetOption("COLZ");

  HitDistFT = tfs->make<TH2F>("hFeedthroughSide","Feedthrough Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-308-19.1*inch,-308-19.1*inch+60*10.89);
  HitDistFT->GetXaxis()->SetTitle("Length along the beam (cm)");
  HitDistFT->GetYaxis()->SetTitle("Height (cm)");
  HitDistFT->GetZaxis()->SetTitle("Entries/bin");
  HitDistFT->SetOption("COLZ");

  HitDistPipe = tfs->make<TH2F>("hPipeSide","Pipe Side",125,-704+205*inch,-704+205*inch+125*10.89,60,-294-19.1*inch,-294-19.1*inch+60*10.89);
  HitDistPipe->GetXaxis()->SetTitle("Length along the beam (cm)");
  HitDistPipe->GetYaxis()->SetTitle("Height (cm)");
  HitDistPipe->GetZaxis()->SetTitle("Entries/bin");
  HitDistPipe->SetOption("COLZ");

  HitDistTop = tfs->make<TH2F>("hTop","Top",125,-701+205*inch,-701+205*inch+125*11.38,80,2-170-300+50.4*inch,2-170-300+50.4*inch+80*11.38);
  HitDistTop->GetXaxis()->SetTitle("Lenght along the beam (cm)");
  HitDistTop->GetYaxis()->SetTitle("Lenght along the drift (cm)");
  HitDistTop->GetZaxis()->SetTitle("Entries/bin");
  HitDistTop->SetOption("COLZ");



}

void TrackDump::endJob()
{
  // Implementation of optional member function here.
  
  
  //fTree->Write();
  
}


void TrackDump::ResetVars()
{
  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime_sec = -99999;
  evttime_nsec = -99999;
  evttime_GPS_sec = -99999;
  evttime_GPS_nsec = -99999;
  nCRThits = 0;
  for (int i = 0; i<kMaxCRThits; ++i){
    hit_plane[i] = -999;
    hit_time_s[i] = -99999.;
    hit_time0[i] = -99999.;
    hit_time1[i] = -99999.;
    hit_charge[i] = -99999.;
    hit_posx[i] = -99999.;
    hit_posy[i] = -99999.;
    hit_posz[i] = -99999.;
  }

  nCRTtzeros = 0;
  for (int i = 0; i<kMaxCRTtzeros; ++i){
    tz_time_s[i] = -99999.;
    tz_time0[i] = -99999.;
    tz_time1[i] = -99999.;
  }


  nCRTtracks=0;
  for (int j = 0; j<kMaxCRTtracks; ++j){
    ct_theta[j]=-99999.;
    ct_phi[j]=-99999.;
    ct_length[j]=-99999.;
    ct_time_sec[j]=-99999.;
    ct_time0[j]=-99999.;
    ct_time1[j]=-99999.;
    ct_x1[j]=-99999.;
    ct_y1[j]=-99999.;
    ct_z1[j]=-99999.;
    ct_x2[j]=-99999.;
    ct_y2[j]=-99999.;
    ct_z2[j]=-99999.;
  }

  if (fSaveTPCTrackInfo) {
    nTPCtracks=0;
    for (int i = 0; i<kMaxTPCtracks; ++i){
      trkstartx[i]=-9999.;
      trkstarty[i]=-9999.;
      trkstartz[i]=-9999.;
      trkendx[i]=-9999.;
      trkendy[i]=-9999.;
      trkendz[i]=-9999.;
      trkstartdcosx[i]=-9999.;
      trkstartdcosy[i]=-9999.;
      trkstartdcosz[i]=-9999.;
      trkenddcosx[i]=-9999.;
      trkenddcosy[i]=-9999.;
      trkenddcosz[i]=-9999.;
      trktheta[i]=-9999.;
      trkphi[i]=-9999.;
      trklen[i]=-9999.;
    }
  }
  
  if (fSavePMTFlashInfo) {
    
    nPMTflash=0;
    for (int i = 0; i<kMaxPMTflashes; ++i){
      Yflash[i]=-9999.;
      Zflash[i]=-9999.;
      PEflash[i]=-9999.;
      Timeflash[i]=-9999.;
      fbeam[i]=-9999.;
    }
  }


}
DEFINE_ART_MODULE(TrackDump)


