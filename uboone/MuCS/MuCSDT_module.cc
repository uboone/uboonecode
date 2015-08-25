
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
    
    Int_t ntot = 0;
    
    Int_t no_planes;
    
    Int_t no_wires[3];
    
    Float_t wy[3456];
    Float_t wt[3456];
    
    // Float_t  wu[2399];
    
    // Float_t  wv[2399];
        
    std::string fSimulationProducerLabel; 
    
    std::string fBezierRecoProducerLabel;
    
    std::string fHitRecoProducerLabel;
        
    std::vector<Float_t> *muE = new std::vector<Float_t>; 
    
    std::vector<Float_t> *muPx = new std::vector<Float_t>; std::vector<Float_t> *muPy = new std::vector<Float_t>; std::vector<Float_t> *muPz = new std::vector<Float_t>;
    
    Double_t p = -1;
        
    std::vector<Float_t> *muX = new std::vector<Float_t>; std::vector<Float_t> *muY = new std::vector<Float_t>; std::vector<Float_t> *muZ = new std::vector<Float_t>;
    
    std::vector<Float_t> *muT = new std::vector<Float_t>;
    
    std::vector<Float_t> *recoX = new std::vector<Float_t>; std::vector<Float_t> *recoY = new std::vector<Float_t>; std::vector<Float_t> *recoZ = new std::vector<Float_t>;
    
    Int_t flag2 = 0;
    
    Double_t L0 = -1;
    
    Double_t p0 = -1;
    
    Double_t p1 = -1;
        
    Int_t sw = 0;
    
    Double_t LLHDp = -1;
    
    std::vector<Float_t> *hit_yz = new std::vector<Float_t>;
    std::vector<Float_t> *hit_yx = new std::vector<Float_t>; std::vector<Float_t> *hit_yx_e = new std::vector<Float_t>;
    std::vector<Float_t> *hit_yq = new std::vector<Float_t>;
    /*
    std::vector<Float_t> *hit_uz = new std::vector<Float_t>; 
    std::vector<Float_t> *hit_ux = new std::vector<Float_t>; std::vector<Float_t> *hit_ux_e = new std::vector<Float_t>;
    
    std::vector<Float_t> *hit_vz = new std::vector<Float_t>; 
    std::vector<Float_t> *hit_vx = new std::vector<Float_t>; std::vector<Float_t> *hit_vx_e = new std::vector<Float_t>;
    */      
    
    art::ServiceHandle<geo::Geometry> geo;
    
    Double_t x1 = 0.0; Double_t x2 = 256.35;
    
    Double_t y1 = -116.5; Double_t y2 = 116.5;
    
    Double_t z1 = 0.0; Double_t z2 = 1036.8;
    
    TTree *my_tree;
    
    Double_t tickToDist;
    
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
    muE->clear();
    
    muPx->clear(); muPy->clear(); muPz->clear();
    
    muX->clear(); muY->clear(); muZ->clear();
    
    muT->clear();
    
    p = -1;
    
    recoX->clear(); recoY->clear(); recoZ->clear();
    
    L0 = -1;
    
    flag2 = 0;
    
    p0 = -1;
    
    p1 = -1;
        
    sw = 0;
    
    hit_yz->clear(); hit_yx->clear(); hit_yx_e->clear();
    
    hit_yq->clear();
    
    LLHDp = -1;

    
    // hit_uz->clear(); hit_ux->clear(); hit_ux_e->clear();
    
    // hit_vz->clear(); hit_vx->clear(); hit_vx_e->clear();
    
    art::ServiceHandle<art::TFileService> tfs;
    
    my_tree = tfs->make<TTree>( "my_tree", "my_tree" );
    
    my_tree->Branch( "trgId", &trigID, "trgId/I" );
    
    my_tree->Branch( "muE", &muE ); 
    
    my_tree->Branch( "muPx", &muPx ); my_tree->Branch( "muPy", &muPy ); my_tree->Branch( "muPz", &muPz ); 
    
    my_tree->Branch( "muX", &muX ); my_tree->Branch( "muY", &muY ); my_tree->Branch( "muZ", &muZ ); 
    
    my_tree->Branch( "muT", &muT ); 
    
    my_tree->Branch( "p", &p ); 
    
    my_tree->Branch( "recoX", &recoX ); my_tree->Branch( "recoY", &recoY ); my_tree->Branch( "recoZ", &recoZ ); 
    
    my_tree->Branch( "flag2", &flag2 ); 
    
    my_tree->Branch( "L", &L0 ); 
    
    my_tree->Branch( "pmcs0", &p0 ); 
    
    my_tree->Branch( "pmcs1", &p1 ); 
    
    my_tree->Branch( "sw", &sw ); 
    
    my_tree->Branch( "hit_yz", &hit_yz ); my_tree->Branch( "hit_yx", &hit_yx ); my_tree->Branch( "hit_yx_e", &hit_yx_e );
    
    my_tree->Branch( "hit_yq", &hit_yq );
    
    my_tree->Branch( "LLHDp", &LLHDp );
    
  }
    
  void MuCSDT::beginRun( const art::Run& run )
  {}
  
  void MuCSDT::reconfigure( fhicl::ParameterSet const& p )
  {
    fSimulationProducerLabel = p.get< std::string >( "SimulationLabel" );
    
    fBezierRecoProducerLabel = p.get< std::string >( "BezierRecoLabel" );
    
    fHitRecoProducerLabel = p.get< std::string >( "HitLabel" );
                
    return;
  
  }
    
  void MuCSDT::analyze( const art::Event& evt ) 
  {
    if ( trigID==0 )
      {
	cout << "" << endl;
	
	cout << " starting ... " << endl;
	
	cout << "" << endl;
      
	no_planes = geo->Nplanes( 0, 0 );
      
	cout << " * planes        : " << no_planes << endl;
	
	cout << "" << endl; 
	
	no_wires[0] = geo->Nwires( 0, 0, 0 ); 
	
	no_wires[1] = geo->Nwires( 1, 0, 0 ); 
	
	no_wires[2] = geo->Nwires( 2, 0, 0 ); 
	
	cout << " * U-plane ( 0 ) : " << no_wires[0] << " wires " << endl; 
	
	cout << "" << endl;
	
	cout << " * V-plane ( 1 ) : " << no_wires[1] << " wires " << endl; 
	
	cout << "" << endl;
	
	cout << " * Y-plane ( 2 ) : " << no_wires[2] << " wires " << endl; 
	
	cout << "" << endl;
	
	cout << "" << endl; 
	
	
	Double_t e_1[3]; 
	
	Double_t e_2[3];
	
	for ( int jjj=0; jjj<no_wires[2]; jjj++ ) 
	  {
	    geo->WireEndPoints( 0, 0, 2, jjj, e_1, e_2 );
	    
	    // cout << e_1[0] << ", " << e_1[1] << ", " << e_1[2] << endl;
	    // cout << e_2[0] << ", " << e_2[1] << ", " << e_2[2] << endl;
	    	    
	    // output << setw(15) << 4798+jjj << setw(15) << jjj << setw(15) << 2 << setw(15) << e_1[0] 
	    // << setw(15) << e_1[1] << setw(15) << e_1[2] << setw(15) << e_2[0] << setw(15) << e_2[1] << setw(15) << e_2[2] << endl;
	    
	    wy[jjj] = e_1[2];
	    wt[jjj] = e_1[0];
	    
	    // getchar();
	    
	  }
	
			
	Int_t nop = geo->Cryostat( 0 ).NOpDet();
	
	cout << " Number of op. detectors : " << nop << endl;  
	
	cout << "" << endl;  
		
	cout << "" << endl;
	
	getchar();
	
      }
    
    muE->clear();
    
    muPx->clear(); muPy->clear(); muPz->clear();
    
    muX->clear(); muY->clear(); muZ->clear();
    
    muT->clear();
    
    recoX->clear(); recoY->clear(); recoZ->clear();
    
    flag2 = 0;
        
    hit_yz->clear(); hit_yx->clear(); hit_yx_e->clear();
    
    hit_yq->clear();
    
    LLHDp = -1;
    
    // hit_uz->clear(); hit_ux->clear(); hit_ux_e->clear();
    
    // hit_vz->clear(); hit_vx->clear(); hit_vx_e->clear();
    
    Double_t xTRUE[3];
    
    xTRUE[0]=100000.0; xTRUE[1]=100000.0; xTRUE[2]=100000.0; 
    
    xTRUE[0]*=1.0; xTRUE[1]*=1.0; xTRUE[2]*=1.0; 
            
    std::vector< art::Handle< std::vector<simb::MCParticle> > > mclist;
    
    evt.getManyByType( mclist );
    
    art::Handle< std::vector<simb::MCParticle> > mc = mclist[0]; 
    
    for ( UInt_t i=0; i<mc->size(); i++ )
      {  
	string is_pr = mc->at(i).Process(); 
	
	Int_t pr_t = mc->at(i).PdgCode();
		
	if ( is_pr.compare( "primary" )==0 && TMath::Abs( pr_t )==13 )
	  {
	    Int_t last = mc->at(i).NumberTrajectoryPoints();
	    
	    for ( int j=0; j<last; j++ )
	      {
		const TLorentzVector &pos = mc->at(i).Position( j );
		
		if ( j==0 ) 
		  { 
		    xTRUE[0]=pos.X(); 
		    xTRUE[1]=pos.Y(); 
		    xTRUE[2]=pos.Z(); 
		    		    
		    const TLorentzVector &pp = mc->at(i).Momentum( j ); 		    
		    p = sqrt( pow( pp.Px(), 2.0 )+pow( pp.Py(), 2.0 )+pow( pp.Pz(), 2.0 ) );
		    
		  }
						
		if ( pos.X()>=x1 && pos.X()<=x2 && pos.Y()>=y1 && pos.Y()<=y2 && pos.Z()>=z1 && pos.Z()<=z2  )
		  {
		    const TLorentzVector &ene = mc->at(i).Momentum( j ); 
		    
		    Double_t eee = ene.E(); 
		    
		    muE->push_back( eee ); 
		    
		    muPx->push_back( ene.Px() ); muPy->push_back( ene.Py() ); muPz->push_back( ene.Pz() );
		    
		    muX->push_back( pos.X() ); muY->push_back( pos.Y() ); muZ->push_back( pos.Z() ); 
		    
		    muT->push_back( pos.T() );
		    
		  }
		
		else { flag2=2; }
		
	      }
	    
	  }
	
      }
    
    // cout << "t0: " << muT->at( 0 ) << endl; 
        
    art::Handle< std::vector<recob::Track> > trackHandle;
    
    evt.getByLabel( fBezierRecoProducerLabel, trackHandle );
    
    std::vector< art::Ptr<recob::Track> > tracks;
    
    art::fill_ptr_vector( tracks, trackHandle );
    
    UInt_t trno = tracks.size();
    
    Double_t xRECO[3]; 
    
    Int_t doit = -1; doit*=1.0;
    
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<util::DetectorProperties> detp;
    tickToDist = larp->DriftVelocity( larp->Efield(), larp->Temperature() );
    tickToDist *= ( 1.e-3*detp->SamplingRate() );
    
    for( size_t t=0; t<trno; t++ )
      {
	art::Ptr<recob::Track> trk = tracks.at( t );
	std::vector< art::Ptr<recob::Hit> > hits = art::FindManyP<recob::Hit>( tracks, evt, fHitRecoProducerLabel ).at( t );
		
	UInt_t nhits = hits.size(); 
	UInt_t npoints = trk->NumberTrajectoryPoints();
	// cout << nhits << ", " << npoints << endl; getchar();
	
	const TVector3 &VP = trk->LocationAtPoint( 0 );
	
	xRECO[0]=VP.X(); xRECO[1]=VP.Y(); xRECO[2]=VP.Z(); 
	
	Double_t xR = sqrt( pow( xRECO[0]-xTRUE[0], 2.0 )+pow( xRECO[1]-xTRUE[1], 2.0 )+pow( xRECO[2]-xTRUE[2], 2.0 ) );
	
	Double_t LL0 = trk->Length();
	
	recoX->clear(); recoY->clear(); recoZ->clear();
	
	hit_yz->clear(); hit_yx->clear(); hit_yx_e->clear();
	hit_yq->clear();
	// hit_uz->clear(); hit_ux->clear(); hit_ux_e->clear();
	// hit_vz->clear(); hit_vx->clear(); hit_vx_e->clear();
    
	if ( xR<=20.0 && LL0>17.5 )
	  {
	    ntot++;
	    	    
	    doit = 1;
	    
	    L0 = trk->Length();
	    
	    Double_t pp0, pp1, pp2;
	    
	    sw = 0;
	    
	    for ( UInt_t u=0; u<npoints; u++ )
	      {
		const TVector3 &muP = trk->LocationAtPoint( u );
		
		recoX->push_back( muP.X() ); recoY->push_back( muP.Y() ); recoZ->push_back( muP.Z() ); 
		
		if ( u!=0 )
		  {
		    Double_t che = sqrt( pow( muP.X()-pp0, 2.0 )+pow( muP.Y()-pp1, 2.0 )+pow( muP.Z()-pp2, 2.0 ) );
		    
		    if ( che>5.0 ) sw++;
		    
		  }
		  
		pp0 = muP.X(); pp1 = muP.Y(); pp2 = muP.Z();
				
	      }
	    	    
	    Int_t NoY=0;Int_t NoU=0;Int_t NoV=0;
	    
	    // cout << " 1 " << endl;
	    
	    for ( UInt_t ii=0; ii<nhits; ii++ )
	      {  
		UInt_t curr_plane = hits.at(ii)->WireID().Plane; 
		
		UInt_t curr_wire = hits.at(ii)->WireID().Wire; 
				
		if ( curr_plane==0 ) NoU++;
		if ( curr_plane==1 ) NoV++;
		
		if ( curr_plane==2 ) 
		  {
		    NoY++;
		    
		    Float_t co_z = wy[ curr_wire ];
		    
		    hit_yz->push_back( co_z );
		    
		    Float_t DrT = hits.at(ii)->PeakTime(); 
		    
		    hit_yx->push_back( ( DrT )*tickToDist-257.36806744 ); // T-muT->at(0)
		    
		    Float_t DrTe = hits.at(ii)->RMS(); 
		    
		    hit_yx_e->push_back( ( DrTe )*tickToDist );
		    
		    hit_yq->push_back( hits.at(ii)->Integral() );
		    		    
		  }
		
	      }
	    
	    // cout << " 2 " << endl;
	    
	    // cout << " ---> " << NoU << ", " << NoV << ", " << NoY << " ---> " <<  NoU+NoV+NoY << endl; getchar();
	    
	    // !!!
	    
	    trkf::TrackMomentumCalculator trkm;
	    Double_t pmcs0 = trkm.GetMomentumMuCSDTChi2( trk );	    
	    Double_t pmcs = trkm.GetMomentumMuCSDTLLHD( trk );
	    p0 = pmcs0; p1 = pmcs; cout << " : " << pmcs0 << " : " << pmcs << endl; 
	    if ( isinf( float(pmcs0) ) || isnan( float(pmcs0) ) || isinf( float(pmcs) ) || isnan( float(pmcs) ) ) { cout << " Problem !" << endl; getchar(); }
	    
	    LLHDp = trkm.GetMuMuCSDTLLHD3( trk, true );
	    
	    // <--- !!!
	    
	  }
	
      }
    
    if ( doit!=-1 ) my_tree->Fill();
    
    // my_tree->Fill();
        
    trigID++;
    
    return;
    
  }
  
  void MuCSDT::endJob()
  {
    cout << "--> !" << endl; 
    
    cout << ntot << ", " << ( 1.0*ntot*1.0 )/( 1.0*trigID*1.0 )*100.0 << endl;
    
    cout << "" << endl;
        
    cout << " ... ending ! " << endl;
    
    cout << "" << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSDT )
  
} 

#endif 
