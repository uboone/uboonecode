
////////////////////////////////////////////////////////////////////////
// Class:       MuCSMerger
// Module Type: producer
// File:        MuCSMerger_module.cc
//
// Generated at Wed May 20 14:08:15 2015 by Leonidas N. Kalousis using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//    trivia : The main routine developed to merge data from the TPC 
//             and the Muon Counter System (MuCS), September 2015
//    author : Leonidas N. Kalousis
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MUCSMERGER_H
#define MUCSMERGER_H 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// #include "EventDisplay/HeaderDrawer.h"
// #include "EventDisplayBase/View2D.h"
// #include "EventDisplayBase/EventHolder.h"
#include "TText.h"
#include "TTimeStamp.h"

#include <memory>
#include <iostream>
#include "MuCSData.h"
#include "RawData/TriggerData.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

class MuCSMerger;

class MuCSMerger : public art::EDProducer {
public:
  explicit MuCSMerger( fhicl::ParameterSet const &pset );
  virtual ~MuCSMerger();
  
  void reconfigure( fhicl::ParameterSet const &pset ); // override;
  void produce( art::Event &evt ) override;
      
  private:
  
  std::string fSwizzlerProducerLabel; 
  
  Int_t group = 0;
  Int_t run = 0;
  
  Int_t trigID = 0;
  TH1F *hDT;
  Int_t run0;
  Int_t srun0;
  TFile *f1;
  TFile *f2;
  
  Int_t seq;
  Float_t time_sec_low;
  Float_t time_sec_high;
  Float_t time_16ns_low;
  Float_t time_16ns_high;
  Double_t t0;
  
  TTree *my_tree;
  Int_t my_entries;
  
  Double_t previous_trigtime;
  Double_t t_start; 
  
  Double_t TOLER = 20.0;
  Double_t offset = -666.0;
  Double_t TOLER2 = 1.6*0.001;

};

void MuCSMerger::reconfigure( fhicl::ParameterSet const &p )
  {
    fSwizzlerProducerLabel = p.get< std::string >( "SwizzlerProducerLabel" );
    
    group = p.get< int >( "group" );
    
    run = p.get< int >( "run" );
    
    return;
    
  }

MuCSMerger::MuCSMerger( fhicl::ParameterSet const &pset )
// :
// Initialize member data here.
{
  this->reconfigure( pset );
  
  produces< std::vector<MuCS::MuCSData> >();  
  
}

MuCSMerger::~MuCSMerger()
{}

void MuCSMerger::produce( art::Event &evt )
{
  if ( trigID==0 )
    {
      cout << "" << endl;
      cout << " starting ... " << endl;
      cout << "" << endl;
      
      cout << " - group : " << group << endl;
      cout << "" << endl;
      
      f1 = new TFile( Form( "/uboone/data/users/kalousis/MuCS/offsets/MuCSDT_%d_%d.root", run, group ), "read" ); 
      TDirectory *dir = f1->GetDirectory( "MuCSDT" );
      hDT = (TH1F*)dir->Get( "hDT" );
      
      offset = hDT->GetXaxis()->GetBinCenter( hDT->GetMaximumBin() );
      // offset = hDT->GetXaxis()->GetBinUpEdge( hDT->GetMaximumBin() );
      cout << " - DT : " << Form( "%.6f", offset ) << endl;
      cout << "" << endl;
            
      f2 = new TFile( Form( "/uboone/data/users/kalousis/MuCS/muons/mega_micro_ana_%d_0.333_0.root", group ), "read" );  
    
      if ( f2->IsZombie() ) 
	{
	  cout << " - mucs file not existing ! " << endl;
	  return;
	  
	}
      
      else
	{
	  my_tree = (TTree*)f2->Get( "preselected" );
	  
	  my_tree->SetBranchStatus( "*", 0 ); 
	  my_tree->SetBranchStatus( "seq", 1 );
	  my_tree->SetBranchStatus( "time_sec_low", 1 );
	  my_tree->SetBranchStatus( "time_sec_high", 1 );
	  my_tree->SetBranchStatus( "time_16ns_low", 1 );
	  my_tree->SetBranchStatus( "time_16ns_high", 1 );
	  my_tree->SetBranchStatus( "t0", 1 );
	  
	  my_tree->SetBranchAddress( "seq", &seq );
	  my_tree->SetBranchAddress( "time_sec_low", &time_sec_low );
	  my_tree->SetBranchAddress( "time_sec_high", &time_sec_high );
	  my_tree->SetBranchAddress( "time_16ns_low", &time_16ns_low );
	  my_tree->SetBranchAddress( "time_16ns_high", &time_16ns_high );
	  my_tree->SetBranchAddress( "t0", &t0 );
	  
	  my_entries = my_tree->GetEntries();
	  cout << " - events in mucs : " << my_entries << endl;
	  cout << "" << endl;
	  cout << "" << endl;
	  
	}
      
      run0 = evt.run();
      srun0 = evt.subRun();
      
      previous_trigtime = 0.0;
            
    }
  
  Int_t event = evt.id().event();
  cout << "" << endl;
  cout << " - run : " << run0 << ", sub. run : " << srun0 << ", event : " << event << endl;
  cout << "" << endl;
  
  unsigned long long int tsval = evt.time().value();
  const unsigned long int mask32 = 0xFFFFFFFFUL;
  unsigned long int unix_time_stamp = ( tsval >> 32 ) & mask32;
  // unsigned long int llo = tsval & mask32;
  // TTimeStamp ts(unix_time_stamp, (int)llo);
  cout << " - unix timestamp : " << unix_time_stamp << endl;
  cout << "" << endl; 
  
  art::Handle< std::vector<raw::Trigger> > trigHandle;
  evt.getByLabel( fSwizzlerProducerLabel, trigHandle );
  std::vector< art::Ptr<raw::Trigger> > trigs;
  art::fill_ptr_vector( trigs, trigHandle );
  Double_t trigtime = ( trigs.at(0)->TriggerTime() )*1.0e-6;
  cout << " - trig. time : " << trigtime << ", diff : " << trigtime-previous_trigtime <<  ", " << trigs.size() << endl;
  cout << "" << endl; 
  
  if ( trigID==0 ) t_start = trigtime;
  Double_t t_rel = trigtime-t_start;
  cout << " - relative trig. time : " << t_rel << endl;
  cout << "" << endl; 
  
  Int_t ntimes=0;
  
  for ( Int_t i=0; i<my_entries; i++ )
      {
	my_tree->GetEntry( i );
	
	Float_t DTunix = TMath::Abs( time_sec_high*65536.0+time_sec_low-unix_time_stamp );
		
	if ( DTunix<=TOLER ) 
	  {
	    Double_t tmucs = t0*1.0e-9;
	    Double_t dt0 =  tmucs-t_rel;
	    Double_t dt = dt0-offset; 
	    
	    if ( TMath::Abs(dt)<TOLER2 )
	      {
		cout << " Gotcha !!! " << endl;
		cout << "" << endl;
		cout << " i : " << i << ", mucs unix timestamp : " << Form( "%.1f", time_sec_high*65536.0+time_sec_low ) << ", diff : " << DTunix << endl; 
		cout << "" << endl;
		cout << " - mucs t0 : " << tmucs << ", " << "diff : " << dt << endl;
		cout << "" << endl;
		
		ntimes++;
		if ( ntimes>=2 ) { cout << " - MULTIPLE PAIRS !!! " << endl; getchar(); }
		
	      }
	    
	  }
		
      }
  
  std::unique_ptr< std::vector<MuCS::MuCSData> > mucsdatacol(new std::vector<MuCS::MuCSData>);
  
  Float_t t0 = 0;
  
  Float_t adc1[24], adc2[24], adc3[24], adc7[24];
  
  std::vector<Int_t> hits1, hits2, hits3, hits7;
  
  t0 = -1.0;  
  for ( Int_t i=0; i<24; i++ ) 
    { 
      adc1[i]=-1.0; adc2[i]=-1.0;
      adc3[i]=-1.0; adc7[i]=-1.0;
      
    }
  hits1.push_back(-1.0);
  hits2.push_back(-1.0);hits2.push_back(-1.0);
  hits3.push_back(-1.0);hits3.push_back(-1.0);hits3.push_back(-1.0);
  hits7.push_back(-1.0);hits7.push_back(-1.0);hits7.push_back(-1.0);hits7.push_back(-1.0);
    
  MuCS::MuCSData mucsevt( t0, adc1, adc2, adc3, adc7, hits1, hits2, hits3, hits7 ); 
    
  
  mucsdatacol->push_back( mucsevt );

  evt.put( std::move( mucsdatacol ) );
  
  previous_trigtime = trigtime;
    
  trigID++;
  return;
    
}

DEFINE_ART_MODULE( MuCSMerger )

#endif

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
