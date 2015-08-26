
#ifndef MuCSDT_Module

#define MuCSDT_Module


#include "Simulation/SimChannel.h"

#include "Simulation/LArG4Parameters.h"


#include "Utilities/LArProperties.h"

#include "Utilities/DetectorProperties.h"


#include "Geometry/Geometry.h"

#include "Geometry/OpDetGeo.h"

#include "SimulationBase/MCParticle.h"

#include "SimulationBase/MCTruth.h"

#include "SimpleTypesAndConstants/geo_types.h"


#include "RecoBase/Hit.h"

#include "RecoAlg/SpacePointAlg.h"

#include "RecoBase/Cluster.h"

#include "RecoBase/Track.h"

#include "RecoBase/PFParticle.h"  

#include "RecoBase/SpacePoint.h"  

#include "RecoBase/OpHit.h"  

#include "RecoBase/OpFlash.h"  

#include "RecoObjects/BezierTrack.h"

#include "RecoAlg/TrackMomentumCalculator.h"


#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
 
#include "art/Framework/Services/Optional/TFileService.h"

#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Core/FindManyP.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"

#include "EventDisplay/HeaderDrawer.h"


#include "TH1.h"

#include "TH2.h"

#include "TTree.h"

#include "TLorentzVector.h"

#include "TVector3.h"


#include <map>

#include <vector>

#include <algorithm>

#include <iostream>

#include <string>

#include <cmath>


using namespace std;


namespace MuCSDT 
{
  class MuCSDT : public art::EDAnalyzer 
  {
  public:
    
    explicit MuCSDT( fhicl::ParameterSet const& pset );
    
    virtual ~MuCSDT();
    
    void beginJob();
        
    void beginRun( const art::Run& run );
        
    void reconfigure( fhicl::ParameterSet const& pset );
        
    void analyze ( const art::Event& evt ); 
    
    void endJob();
    
  private:
    
    Int_t trigID = 0;
    TH1F *hDT;
    Int_t run0;
    Int_t srun0;
    TFile *f1;
        
    Int_t seq;
    Float_t time_sec_low;
    Float_t time_sec_high;
    Float_t time_16ns_low;
    Float_t time_16ns_high;
    Double_t t0;
    
    Int_t my_entries;
    
  }; 
  
  MuCSDT::MuCSDT( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  
  }
  
  MuCSDT::~MuCSDT() 
  {}
    
  void MuCSDT::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    hDT = tfs->make<TH1F>( "hDT", "", 100, -100, 100 );
        
  }
    
  void MuCSDT::beginRun( const art::Run& run )
  {}
  
  void MuCSDT::reconfigure( fhicl::ParameterSet const& p )
  {
    return;
    
  }
    
  void MuCSDT::analyze( const art::Event& evt ) 
  {
    if ( trigID==0 )
      {
	cout << "" << endl;
	cout << " starting ... " << endl;
	cout << "" << endl;
	
	run0 = evt.run();
	srun0 = evt.subRun();
	
	f1 = new TFile( Form( "/uboone/data/users/kalousis/MuCS/muons/mega_micro_ana_669_0.333_0.root" ), "read" );  
	if ( f1->IsZombie() ) 
	  {
	    cout << " - mucs file not existing ! " << endl;
	    return;
	    
	  }
	
	else
	  {
	    TTree *my_tree = (TTree*)f1->Get( "preselected" );
	    
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
	    getchar();
	    
	  }
	
	
      }
    
    Int_t event = evt.id().event();
    cout << "" << endl;
    cout << " - run : " << run0 << ", sub. run : " << srun0 << ", event : " << event << endl;
    cout << "" << endl;
    
    unsigned long long int tsval = evt.time().value();
    const unsigned long int mask32 = 0xFFFFFFFFUL;
    unsigned long int lup = ( tsval >> 32 ) & mask32;
    // unsigned long int llo = tsval & mask32;
    // TTimeStamp ts(lup, (int)llo);
    
    cout << " - unix timestamp : " << lup << endl;
    cout << "" << endl;
    /*
    art::ServiceHandle< util::TimeService > ts;
    Double_t trigT = ts->TriggerClock().Time();
    cout << " - " << trigT-10.0 << endl;
    cout << "" << endl;
    */
        
    getchar();
    trigID++;
    return;
    
  }
  
  void MuCSDT::endJob()
  {
    cout << "" << endl; 
    cout << " - events processed : " << trigID << endl;
    cout << "" << endl;
    cout << " ... ending ! " << endl;
    cout << "" << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSDT )
  
} 

#endif 
