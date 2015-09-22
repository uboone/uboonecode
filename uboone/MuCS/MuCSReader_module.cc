
////////////////////////////////////////////////////////////////////////
//
//    trivia : A simple analyser to read merged data, September 2015
//    author : Leonidas N. Kalousis
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MuCSReader_Module

#define MuCSReader_Module


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

#include "RawData/TriggerData.h"


#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
 
#include "art/Framework/Services/Optional/TFileService.h"

#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Core/FindManyP.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"

// #include "EventDisplay/HeaderDrawer.h"


#include "TMath.h"

#include "TH1.h"

#include "TAxis.h"

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


#include "MuCSData.h"

using namespace std;


namespace MuCSReader 
{
  class MuCSReader : public art::EDAnalyzer 
  {
  public:
    
    explicit MuCSReader( fhicl::ParameterSet const& pset );
    
    virtual ~MuCSReader();
    
    void beginJob();
        
    void beginRun( const art::Run& run );
        
    void reconfigure( fhicl::ParameterSet const& pset );
        
    void analyze ( const art::Event& evt ); 
    
    void endJob();
    
  private:
    
    std::string fMergerProducerLabel; 
    Int_t group = 0;
    Int_t trigID = 0;
    Int_t run0;
    Int_t srun0;
    
  }; 
  
  MuCSReader::MuCSReader( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  
  }
  
  MuCSReader::~MuCSReader() 
  {}
    
  void MuCSReader::beginJob()
  {}
  
  void MuCSReader::beginRun( const art::Run& run )
  {}
  
  void MuCSReader::reconfigure( fhicl::ParameterSet const& p )
  {
    fMergerProducerLabel = p.get< std::string >( "MergerProducerLabel" );
        
    return;
    
  }
    
  void MuCSReader::analyze( const art::Event& evt ) 
  {
    if ( trigID==0 )
      {
	cout << "" << endl;
	cout << " starting ... " << endl;
	cout << "" << endl;
		
	run0 = evt.run();
	srun0 = evt.subRun();
		
      }
    
    Int_t event = evt.id().event();
    cout << "" << endl;
    cout << " - run : " << run0 << ", sub. run : " << srun0 << ", event : " << event << endl;
    cout << "" << endl;
    /*
    art::Handle< std::vector<MuCS::MuCSData> > mucsHandle;
    evt.getByLabel( fMergerProducerLabel, mucsHandle );
    std::vector< art::Ptr<MuCS::MuCSData> > mucs;
    art::fill_ptr_vector( mucs, mucsHandle );
    */
    std::vector< art::Handle< std::vector<MuCS::MuCSData> > > mucslist;
    evt.getManyByType( mucslist );
    art::Handle< std::vector<MuCS::MuCSData> > mucs = mucslist[0]; 
    cout << mucs->size() << endl;
    
    // Float_t time0 = mucs->at(0).T0();
    // cout << time0 << endl;
        
    trigID++;
    return;
    
  }
  
  void MuCSReader::endJob()
  {
    cout << "" << endl; 
    cout << " - events processed : " << trigID << endl;
    cout << "" << endl;
    cout << " ... ending ! " << endl;
    cout << "" << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSReader )
  
} 

#endif 

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
