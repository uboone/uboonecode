////////////////////////////////////////////////////////////////////////
// Class:       CRTTiming
// Module Type: analyzer
// File:        CRTTiming_module.cc
//
// Generated at Tue Oct 30 11:45:06 2018 by David Lorca galindo using artmod
// from cetpkgsupport v1_14_01.
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

#include "uboone/CRT/CRTProducts/CRTHit.hh"
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

namespace crt {
  class CRTTiming;
}

class crt::CRTTiming : public art::EDAnalyzer {
public:
  explicit CRTTiming(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTiming(CRTTiming const &) = delete;
  CRTTiming(CRTTiming &&) = delete;
  CRTTiming & operator = (CRTTiming const &) = delete;
  CRTTiming & operator = (CRTTiming &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  art::ServiceHandle<art::TFileService> tfs;

  std::string  data_label_crthit_;
  int fHardDelay_;

  TH1F* hCRTHits_T1dis;

  TH1F* hCRTHits_T1dis_B;//bottom
  TH1F* hCRTHits_T1dis_F;//Feedthrough
  TH1F* hCRTHits_T1dis_P;//Pipe
  TH1F* hCRTHits_T1dis_T;//Top

};


crt::CRTTiming::CRTTiming(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  data_label_crthit_(p.get<std::string>("data_label_crthit")),
  fHardDelay_(p.get<int>("fHardDelay",40000))
 // More initializers here.
{}

void crt::CRTTiming::analyze(art::Event const & evt)
{

  //get CRTHits                                                                 
  art::Handle< std::vector<crt::CRTHit> > rawHandle_hit;
  evt.getByLabel(data_label_crthit_, rawHandle_hit); //                                                                                        

  //check to make sure the data we asked for is valid                                                                                      
  if(!rawHandle_hit.isValid()){
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << evt.event() << " has zero"
              << " CRTHits " << " in module " << data_label_crthit_ << std::endl;
    std::cout << std::endl;
    return;
  }

  //get better access to the data    //CRTHit collection on this event                                                                     
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_hit);

  //get CRTHits 

  if(CRTHitCollection.size()>0){//A

    for(std::vector<int>::size_type j = 0; j != CRTHitCollection.size(); j++) {//B
      
      crt::CRTHit my_CRTHit = CRTHitCollection[j];
      
      int Hit_T1_nsec = my_CRTHit.ts1_ns + fHardDelay_;
      //      int plane = my_CRTHit.plane;
      
      hCRTHits_T1dis->Fill(Hit_T1_nsec);
      
      if(my_CRTHit.plane == 0){
	hCRTHits_T1dis_B->Fill(Hit_T1_nsec);
      }

      if(my_CRTHit.plane == 1){
	hCRTHits_T1dis_F->Fill(Hit_T1_nsec);
      }

      if(my_CRTHit.plane == 2){
	hCRTHits_T1dis_P->Fill(Hit_T1_nsec);
      }

      if(my_CRTHit.plane == 3){
	hCRTHits_T1dis_T->Fill(Hit_T1_nsec);
      }


    }//B

  }//A


}

void crt::CRTTiming::beginJob()
{
  // Implementation of optional member function here.

  hCRTHits_T1dis = tfs->make<TH1F>("hCRTHits_T1dis","hCRTHits_T1dis",2000,0,20000);
  hCRTHits_T1dis->GetXaxis()->SetTitle("CRTHit time w.r.t. Beam ref. sig. (ns)");
  hCRTHits_T1dis->GetYaxis()->SetTitle("Entries/bin");

  hCRTHits_T1dis_B = tfs->make<TH1F>("hCRTHits_T1dis_B","hCRTHits_T1dis_B",2000,0,20000);
  hCRTHits_T1dis_B->GetXaxis()->SetTitle("CRTHit (Bottom) time w.r.t. Beam ref. sig. (ns)");
  hCRTHits_T1dis_B->GetYaxis()->SetTitle("Entries/bin");

  hCRTHits_T1dis_F = tfs->make<TH1F>("hCRTHits_T1dis_F","hCRTHits_T1dis_F",2000,0,20000);
  hCRTHits_T1dis_F->GetXaxis()->SetTitle("CRTHit (Feedthrough) time w.r.t. Beam ref. sig. (ns)");
  hCRTHits_T1dis_F->GetYaxis()->SetTitle("Entries/bin");

  hCRTHits_T1dis_P = tfs->make<TH1F>("hCRTHits_T1dis_P","hCRTHits_T1dis_P",2000,0,20000);
  hCRTHits_T1dis_P->GetXaxis()->SetTitle("CRTHit (Pipe) time w.r.t. Beam ref. sig. (ns)");
  hCRTHits_T1dis_P->GetYaxis()->SetTitle("Entries/bin");

  hCRTHits_T1dis_T = tfs->make<TH1F>("hCRTHits_T1dis_T","hCRTHits_T1dis_T",2000,0,20000);
  hCRTHits_T1dis_T->GetXaxis()->SetTitle("CRTHit (Top) time w.r.t. Beam ref. sig. (ns)");
  hCRTHits_T1dis_T->GetYaxis()->SetTitle("Entries/bin");


}

void crt::CRTTiming::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(crt::CRTTiming)
