
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
	
	getchar();
	
      }
    
    cout << " - event no. : " << trigID << endl;
    
    cout << "" << endl;
    
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
