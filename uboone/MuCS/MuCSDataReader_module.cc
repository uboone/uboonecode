
////////////////////////////////////////////////////////////////////////
//
//    trivia : A simple example to read MuCS merged data
//    author : Odysseas Kanenas
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MuCSDataReader_Module

#define MuCSDataReader_Module


#include "fstream"

#include "sstream"

#include "iomanip"


#include "Simulation/SimChannel.h"

#include "Simulation/LArG4Parameters.h"


#include "Utilities/LArProperties.h"

#include "Utilities/DetectorProperties.h"


#include "Geometry/Geometry.h"

#include "Geometry/TPCGeo.h"

#include "Geometry/PlaneGeo.h"

#include "Geometry/WireGeo.h"

#include "Geometry/OpDetGeo.h"

#include "SimulationBase/MCParticle.h"

#include "SimulationBase/MCTruth.h"

#include "SimpleTypesAndConstants/geo_types.h"


#include "RecoBase/Hit.h"

#include "RecoBase/Cluster.h"

#include "RecoBase/Track.h"

#include "RecoBase/PFParticle.h"  

#include "RecoBase/SpacePoint.h"  

#include "RecoBase/OpHit.h"  

#include "RecoBase/OpFlash.h"  

#include "RecoObjects/BezierTrack.h"

#include "RecoAlg/TrackMomentumCalculator.h"

#include "AnalysisBase/ParticleID.h"

#include "uboonecode/uboone/MuCS/MuCSData.h"

/* ... */


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

class MuCSDataReader; 

  class MuCSDataReader : public art::EDAnalyzer 
  {
  public:
    
    explicit MuCSDataReader( fhicl::ParameterSet const& pset );
    
    virtual ~MuCSDataReader();
    
    void beginJob();
        
    void beginRun( const art::Run& run );
        
    void reconfigure( fhicl::ParameterSet const& pset );
        
    void analyze ( const art::Event& evt ); 
    
    void endJob();
    
  private:
        
    std::string fMergerProducerLabel; 
        
  }; 
  
  MuCSDataReader::MuCSDataReader( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  
  }
  
  MuCSDataReader::~MuCSDataReader() 
  {}
    
  void MuCSDataReader::beginJob()
  {
    cout << " Begin ...  " << endl;
    
  }
  
  void MuCSDataReader::beginRun( const art::Run& run )
  {}
  
  void MuCSDataReader::reconfigure( fhicl::ParameterSet const& p )
  {
    fMergerProducerLabel = p.get< std::string >( "MuCSMergerLabel" );
    
    return;
    
  }
    
  void MuCSDataReader::analyze( const art::Event& evt ) 
  {
    art::Handle< std::vector<MuCS::MuCSData> > muHandle;
    
    evt.getByLabel( fMergerProducerLabel, muHandle );
    
    std::vector< art::Ptr<MuCS::MuCSData> > mus;
    
    art::fill_ptr_vector( mus, muHandle );
    
    UInt_t muno = mus.size();
    
    for( size_t t=0; t<muno; t++ )
      {	
	cout << "!" << endl;
	
      }
    
    // ..
    
    return;
    
  }
  
  void MuCSDataReader::endJob()
  {
    cout << " ... ending ! " << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSDataReader )

#endif          
