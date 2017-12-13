////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEAnalyzer
// Plugin Type: analyzer (art v2_05_00)
// File:        ErezCCQEAnalyzer_module.cc
//
// Generated at Wed Jul 12 16:10:16 2017 by Erez Cohen using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"


// ROOT includes
#include "TTree.h"

// C/C++ libraries
#include <memory>
#include <utility>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>

// my pandoraNu track...
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/flash.h"
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"
#include "uboone/ErezCCQEana/MyObjects/pairVertex.h"

// constants
constexpr int debug          = 1;
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTruth      = 100;
constexpr int kMaxNgenie     = 100;
constexpr float kMaxInterTrackDistance = 11; // 11 cm between tracks - maximal distance for clustering
constexpr float EPSILON      = 0.1;   // tollerance for equations

// charge deposition around the vertex in a box of N(wires) x N(time-ticks)
constexpr int N_box_sizes    = 30;
constexpr int MinNwiresBox   = 5;
constexpr int dNwiresBox     = 5;
constexpr int MinNticksBox   = 10;
constexpr int dNticksBox     = 10;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class ErezCCQEAnalyzer; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::ErezCCQEAnalyzer : public art::EDAnalyzer {
public:
    explicit ErezCCQEAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEAnalyzer(ErezCCQEAnalyzer const &) = delete;
    ErezCCQEAnalyzer(ErezCCQEAnalyzer &&) = delete;
    ErezCCQEAnalyzer & operator = (ErezCCQEAnalyzer const &) = delete;
    ErezCCQEAnalyzer & operator = (ErezCCQEAnalyzer &&) = delete;
    
    // Required functions.
    void                   analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                  beginJob () override;
    void               reconfigure (fhicl::ParameterSet const& p) override;
    
    void                 endSubRun (const art::SubRun& sr);

    
    // functionallity
    void         ConstructVertices ();
    void   ClusterTracksToVertices ();
    void           AnalyzeVertices ();
    void          FilterGoodPairVertices ();
    void               TagVertices ();
    void          PrintInformation ();
    bool    TrackAlreadyInVertices (int ftrack_id);
    void       HeaderVerticesInCSV ();
    void       StreamVerticesToCSV ();
    
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=0;
    void Debug (Int_t verobosity_level, std::string text){
        if ( debug > verobosity_level ) cout << text << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    

    
    
private:
    
    // Declare member data here.
    void ResetVars();
    
    // Declare member data here.
    TTree *fTree;
    TTree* fPOTTree;

    Int_t run, subrun, event;
    
    short   isdata;
    
    bool    MCmode;
    
    int     Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    int     Nflashes;
    int     Nvertices;
    int     vertices_ctr;
    
    int     NwiresBox[N_box_sizes], NticksBox[N_box_sizes];
    
    double  pot, pot_total;

    // my objects
    std::vector<PandoraNuTrack>     tracks;
    std::vector<hit>                hits;
    std::vector<GENIEinteraction>   genie_interactions;
    std::vector<pairVertex>         vertices;
    std::vector<flash>              flashes;
    
    
    //Module labels to get data products
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDataSampleLabel;
    std::string fPOTModuleLabel;

    
    //mctruth information
    Int_t    mcevts_truth;    //number of neutrino Int_teractions in the spill
    
    // time stamp
    std::chrono::time_point<std::chrono::system_clock> start_ana_time, end_ana_time;
    
    // output csv file of vertices
    ofstream vertices_file;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::ErezCCQEAnalyzer::ErezCCQEAnalyzer(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::analyze(art::Event const & evt){
    
    ResetVars();
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<cheat::BackTracker> bt;
    //    const sim::ParticleList& plist = bt->ParticleList();
    isdata = evt.isRealData();
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    
    // * hits
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

    // * flashes
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fFlashModuleLabel, flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);

    
    // * tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
    
    // * associations
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

    
    // ----------------------------------------
    // hits information
    // ----------------------------------------
    Nhits = hitlist.size();
    Nhits_stored = std::min(Nhits, kMaxHits);
    for (int i = 0; i < Nhits_stored ; ++i){//loop over hits
        
        hit fhit(
                  hitlist[i]->WireID().Plane    // hit plane
                  ,hitlist[i]->WireID().Wire    // hit wire
                  ,i                            // hit id
                  ,hitlist[i]->PeakTime()       // hit peak time
                  ,hitlist[i]->Integral()       // hit charge
                  );
        
        hits.push_back( fhit );
        
    }
    
    
    // ----------------------------------------
    // flash information
    // ----------------------------------------
    Nflashes = flashlist.size();
    if (debug>0) SHOW(Nflashes);
    for ( int f = 0; f < std::min(Nflashes,kMaxHits); f++ ) {
        
        if (debug>0){ SHOW2(flashlist[f]->Time() , flashlist[f]->TotalPE() ) };
        flash fflash(
                     flashlist[f]->Time(),         // flash time
                     flashlist[f]->TimeWidth(),    // flash time width
                     flashlist[f]->ZCenter(),      // flash z-center
                     flashlist[f]->ZWidth(),       // flash z-width
                     flashlist[f]->YCenter(),      // flash y-center
                     flashlist[f]->YWidth(),       // flash y-width
                     flashlist[f]->TotalPE()       // flash total PE
        );
        
        // keep only in-time flashes for the event
        // following Katherine conditions
        if (debug>0){ SHOW2(fflash.GetTime() , fflash.GetTotalPE()) };
        if( (0.0 < fflash.GetTime()) && (fflash.GetTime() < 10.0) && (6.5 < fflash.GetTotalPE()) ){
            flashes.push_back( fflash );
        }
    }


    // ----------------------------------------
    // tracks information
    // ----------------------------------------
    Ntracks = tracklist.size();
    for(int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i ){
        recob::Track::Point_t start_pos, end_pos;
        std::tie( start_pos, end_pos ) = tracklist[i]->Extent();
        
        PandoraNuTrack track(
                             run , subrun, event        // r/s/e
                             ,tracklist[i]->ID()        // track id
                             ,tracklist[i]->Length()    // length
                             ,TVector3(start_pos.X(),start_pos.Y(),start_pos.Z())   // start position
                             ,TVector3(end_pos.X(),end_pos.Y(),end_pos.Z())         // end position
                            );
        
        double StartLoc[3] = {start_pos.X(), start_pos.Y(), start_pos.Z()};
        
        // U / V / Y coordinates
        for (int plane = 0; plane < 3; plane++){
            
            
            Debug(5,"geo::TPCID tpcID = fGeom->FindTPCAtPosition( StartLoc );");
            geo::TPCID tpcID = geom->FindTPCAtPosition( StartLoc );
            int tpc = 0;
            if (tpcID.isValid) tpc = tpcID.TPC;
            else continue;
            
            // Construct wire ID for this point projected onto the plane
            Debug(5,"Construct wire ID for this point projected onto the plane");
            geo::PlaneID planeID = geo::PlaneID( 0 , tpc , plane ); // cryostat=0
            
            // get start point
            Int_t start_wire = (Int_t) ( geom->WireCoordinate( start_pos.Y() , start_pos.Z() ,  planeID ) );
            Int_t start_time = (Int_t) ( detprop->ConvertXToTicks( start_pos.X() , planeID ) ) ;
            // get end point
            Int_t end_wire = (Int_t) ( geom->WireCoordinate( end_pos.Y() , end_pos.Z() ,  planeID ) );
            Int_t end_time = (Int_t) ( detprop->ConvertXToTicks( end_pos.X() , planeID ) ) ;
            // plug into the track
            track.SetStartEndPlane( plane , start_wire , start_time , end_wire , end_time );
        }

        
        // Hits-Tracks association
        if (fmth.isValid()){
            std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(i);
            for (size_t h = 0; h < vhit.size(); ++h){
                if (vhit[h].key()<kMaxHits){
                    hits.at( vhit[h].key() ).SetTrackKey( tracklist[i].key() );
                }
            }
        }
        
        // PIDa and calorimetric KE
        if (fmcal.isValid()){
            unsigned maxnumhits = 0;
            std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
            for (auto const& calo : calos){
                if (calo->PlaneID().isValid){
                    int plane = calo->PlaneID().Plane;
                    
                    // get the calorimetric kinetic energy of the track
                    track.SetCaloKEPerPlane( plane , calo->KineticEnergy() );
                    
                    // select the best plane as the one with the maximal number of charge deposition points
                    if (calo->dEdx().size() > maxnumhits){
                        maxnumhits = calo->dEdx().size();
                        track.SetBestPlane ( plane );
                        track.SetMaxNHits ( maxnumhits );
                    }
                    // build the PIDa as a fit the reduced Bethe Bloch
                    // dE/dx = A * R^{0.42}
                    double pida = 0;
                    int used_trkres = 0;
                    for (size_t ip = 0; ip < calo->dEdx().size(); ++ip){
                        if (calo->ResidualRange()[ip]<30){
                            pida += calo->dEdx()[ip]*pow( calo->ResidualRange()[ip],0.42);
                            ++used_trkres;
                        }
                    }
                    if (used_trkres) pida /= used_trkres;
                    track.SetPIDaPerPlane( plane , pida );
                }
            }
            track.SetPIDa();
        }
        
        // MC information
        Debug(0,"before MC information (if(!isdata&&fmth.isValid()))");
        SHOW2(isdata,fmth.isValid())
        if (!isdata&&fmth.isValid()){
            Debug(0,"MC information !isdata&&fmth.isValid()");
            
            // Find true track for each reconstructed track
            int TrackID = 0;
            std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);
            
            std::map<int,double> trkide;
            for(size_t h = 0; h < allHits.size(); ++h){
                art::Ptr<recob::Hit> hit = allHits[h];
                std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
                for( size_t e = 0; e < TrackIDs.size(); ++e ){
                    trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
                }
            }
            // Work out which IDE despoited the most charge in the hit if there was more than one.
            double max_Edep = -1;
            for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
                if ((ii->second)>max_Edep){
                    max_Edep = ii->second;
                    TrackID = ii->first;
                }
            }
            // Now have trackID, so get PDG code and T0 etc.
            const simb::MCParticle *particle = bt->TrackIDToParticle( TrackID );
            Debug(0,"particle = bt->TrackIDToParticle( TrackID )");
            SHOW(particle);
            if (particle){
                Debug(0, Form("particle: pdg=%d",particle->PdgCode()));
                track.SetMCpdgCode( particle->PdgCode() );
                track.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                track.SetTruthDirection();
                
                track.SetTruthLength();
                track.SetTruthMomentum( particle -> Momentum() );
                track.SetTruthMother( particle -> Mother() );
                track.SetTruthProcess( particle -> Process() );
                
                const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(particle->TrackId());
                if (mc_truth->Origin() == simb::kBeamNeutrino)      track.SetOrigin( "beam neutrino" );
                else if (mc_truth->Origin() == simb::kCosmicRay)    track.SetOrigin( "cosmic ray" );
            }//if (particle)
            
        }//MC
        
        // flash - matching
        // find the closest flash to the track
        if (flashes.size()){
            float YZdistance = 1000, ClosestFlashYZdistance = 1000;
            for (auto f : flashes){
                YZdistance = track.GetDis2Flash(f);
                if ( YZdistance < ClosestFlashYZdistance ){
                    ClosestFlashYZdistance = YZdistance;
                    track.SetClosestFlash( f );
                }
            }
        }
        tracks.push_back( track );
    }
    
    
    
    // ----------------------------------------
    // MC truth information
    // ----------------------------------------
    if (!isdata){
        
        // * MC truth information
        art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
        std::vector<art::Ptr<simb::MCTruth> > mclist;
        if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
        
        art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
        std::vector<art::Ptr<simb::MCFlux> > fluxlist;
        if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
        art::fill_ptr_vector(fluxlist, mcfluxListHandle);
        
        
        mcevts_truth = mclist.size();
        if (mcevts_truth){
            MCmode = true;
            
            for( int mc_evend_id = 0; (mc_evend_id < mcevts_truth) && (mc_evend_id < kMaxTruth) ; mc_evend_id++ ){
                art::Ptr<simb::MCTruth> mctruth = mclist[mc_evend_id];
                
                if (mctruth->Origin() == simb::kBeamNeutrino){
                    
                    GENIEinteraction genie_interaction( run , subrun , event , mc_evend_id );
                    genie_interaction.SetNuPDG( mctruth->GetNeutrino().Nu().PdgCode() );
                    genie_interaction.SetKinematics(
                                                    mctruth->GetNeutrino().QSqr() // Q2
                                                    ,mctruth->GetNeutrino().W()     // W
                                                    ,mctruth->GetNeutrino().X()     // Bjorken x
                                                    ,mctruth->GetNeutrino().Y()     // y
                                                    ,mctruth->GetNeutrino().CCNC()  // CC=0/NC=1
                                                    ,mctruth->GetNeutrino().Mode()  // QE=0
                                                    );
                    genie_interaction.SetVertexPosition( TVector3(mctruth->GetNeutrino().Nu().Vx()
                                                                  ,mctruth->GetNeutrino().Nu().Vy()
                                                                  ,mctruth->GetNeutrino().Nu().Vz()));
                    
                    genie_interaction.SetNuMomentum( mctruth->GetNeutrino().Nu().Momentum() );
                    genie_interaction.SetLeptonMomentum( mctruth->GetNeutrino().Lepton().Momentum() );
                    genie_interaction.SetMomentumTransfer();
                    
                    // all other particles in the interaction
                    Int_t NgenieParticles = (Int_t)mctruth->NParticles();
                    if ( NgenieParticles ){
                        for( int iPart = 0; iPart < std::min( NgenieParticles , kMaxNgenie ); iPart++ ){
                            const simb::MCParticle& part( mctruth->GetParticle(iPart) );
                            // add a primary to genie-interaction
                            genie_interaction.AddPrimary( part.PdgCode()    // pdg code
                                                         ,part.Momentum()   // 4-momentum
                                                         ,part.StatusCode() // status code
                                                         ,part.Mother()     // mother
                                                         ,part.Process()    // process
                                                         ) ;
                            
                            // match the primary particle with a track
                            for (auto & track : tracks){
                                if (
                                    ( part.PdgCode() == track.GetMCpdgCode() ) // the same particle type (pdg code)
                                    &&
                                    ( (part.Momentum() - track.GetTruthMomentum()).Mag() < EPSILON ) // the same truth momentum
                                    &&
                                    ( (genie_interaction.GetVertexPosition() - track.GetTruthStartPos()).Mag() < EPSILON ) // the same start position
                                    &&
                                    ( genie_interaction.IncludesTrack( track.GetTrackID() ) == false ) // does not already include this track
                                    ) {
                                    // Printf("found a primary-track match! plugging mc_evend_id=%d into track %d",mc_evend_id,track.GetTrackID());
                                    genie_interaction.AddTrack ( track );
                                    track.SetMCeventID( mc_evend_id );
                                }
                            }

                        } // for particle
                    }
                    genie_interaction.SortNucleons();
                    genie_interaction.ComputePmissPrec();
                    genie_interaction.SetTruthTopology();
                    genie_interaction.SetReconstructedTopology();

                    genie_interactions.push_back( genie_interaction );
                    
                }//mctruth->Origin()
            }
        }
    }//is neutrino


    // ----------------------------------------
    // event topology (my-vertex....)
    // ----------------------------------------
    ConstructVertices();
    
    
    PrintInformation();
    fTree -> Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ConstructVertices(){
    
    // cluster all tracks at close proximity to vertices
    ClusterTracksToVertices();
    // analyze these vertices: inter-tracks distances, angles...
    AnalyzeVertices();
    // retain only vertices with pairs of 2-tracks at close proximity
    FilterGoodPairVertices();
    // if its a MC event, tag the vertex by their MC information
    if (MCmode) TagVertices();
    // output to csv file
    StreamVerticesToCSV();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ClusterTracksToVertices(){
    // July-25, 2017
    // cluster all tracks at close proximity to vertices
    // and fix the position of each vertex
    bool    FoundCloseTracks , AlreadySetPosition;
    float   closest_distance_ij;
    TVector3 vertex_position;
    
    for (int i=0; i < Ntracks; i++){
        
        // if (!tracks[i].IsFullyContained) continue;
        if (!tracks[i].IsTrackContainedSoft()) continue;
        
        //        // skip if track was clustered to a vertex by in one of the previous loop steps
        //        if ( TrackAlreadyInVertices( tracks[i].GetTrackID() )) continue;
        
        pairVertex vertex( run, subrun, event , vertices.size() );
        vertex.AddTrack( tracks[i] );
        
        FoundCloseTracks = AlreadySetPosition = false;
        
        for ( int j=i+1 ; j < Ntracks ; j++ ){ // i+1?
            
            // if (!tracks[j].IsFullyContained) continue;
            if (tracks[j].IsTrackContainedSoft() && j!=i){
                
                // if this is the first time we go over these two tracks
                // and they are close enough to define a vertex,
                // we also define the position of their mutual vertex
                if (!AlreadySetPosition){
                    
                    // two close tracks (at a separation distance smaller that max_mu_p_distance)
                    std::string StartOrEnd = "None";
                    closest_distance_ij = tracks[i].ClosestDistanceToOtherTrack(tracks[j],&StartOrEnd);
                    
                    if ( closest_distance_ij < kMaxInterTrackDistance ){
                        
                        vertex.AddTrack( tracks[j] );
                        FoundCloseTracks = true;
                        
                        if (StartOrEnd.compare("Start")==0)     vertex_position = tracks[i].GetStartPos();
                        else if (StartOrEnd.compare("End")==0)  vertex_position = tracks[i].GetEndPos() ;
                        else                                    vertex_position = TVector3(-1000,-1000,-1000) ;
                        
                        vertex.SetPosition( vertex_position );
                        AlreadySetPosition = true;
                    }
                }
                
                // else, namely if we have already clustered a vertex
                // and positioned it in space,
                // we only need to check wether the new track (j) is close enough to this vertex
                // to be also associated with it
                else {
                    // Debug(3, Form("cheking track %d distance from vertex %d ",tracks[j].track_id,c_vertex.vertex_id));
                    if ( tracks[j].DistanceFromPoint(vertex.GetPosition()) < kMaxInterTrackDistance ){
                        
                        // SHOWTVector3(c_vertex.position);
                        // Printf("track %d close enough...",tracks[j].track_id);
                        // SHOW(tracks[j].DistanceFromPoint(c_vertex.position));
                        
                        vertex.AddTrack( tracks[j] );
                    }
                }
            }
        }
        
        Debug( 4 , "if (FoundCloseTracks) {");
        if (FoundCloseTracks) {
            
            // if this is an MC event,
            // match a GENIE interaction to the vertex
            Debug( 5 , "if (MCmode)");
            if (MCmode){

                // 1st (and best) method: match the GENIE interaction by the mc event-id
                // use the mc event-id of the first/second track. This is the mc event-id of the proper GENIE interaction
                Debug( 6 , "(vertex.GetTracks().size()>1){");
                if (vertex.GetTracks().size()>1){
                    int mc_id_t0 = vertex.GetTracks().at(0).GetMCeventID();
                    int mc_id_t1 = vertex.GetTracks().at(1).GetMCeventID();
                    Debug( 7 , "if ( mc_id_t0 >= 0 && mc_id_t0 == mc_id_t1 )");
                    if ( mc_id_t0 >= 0
                        &&
                        mc_id_t0 == mc_id_t1 ){
                        Debug( 8 , "vertex.SetGENIEinfo( genie_interactions.at( mc_id_t0 ) );");
                        vertex.SetGENIEinfo( genie_interactions.at( mc_id_t0 ) );
                        Printf("matched GENIE information");
                        vertex.GetGENIEinfo().Print(); // PRINTOUT
                    }
                }
                Debug( 6 , "} (vertex.GetTracks().size()>1)");
                // the problem with this method is that not allways the two tracks came from
                // the same GENIE interaction
                // for example, happenings in which one track came from a GENIE int. and the other from cosmic...
                // so,
                // 2nd method: look for the closest GENIE interaction in the event, spatially
                GENIEinteraction closest_genie;
                float closest_genie_interaction_distance = 10000; // [cm]
                for (auto genie_interaction : genie_interactions){
                    float genie_distance = (genie_interaction.GetVertexPosition() - vertex.GetPosition()).Mag();
                    if ( genie_distance < closest_genie_interaction_distance ){
                        closest_genie_interaction_distance = genie_distance;
                        closest_genie = genie_interaction;
                        break;
                    }
                }
                Debug( 6 , "} for (auto genie_interaction : genie_interactions)");
                if ( closest_genie_interaction_distance < 1000 ){
                    vertex.SetClosestGENIE( closest_genie );
                }
            }
            Debug( 5 , "} if (MCmode)");
            // plug into vertices list
            vertices.push_back( vertex );
        }
        Debug( 4 , "} if (FoundCloseTracks)");
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::ErezCCQEAnalyzer::TrackAlreadyInVertices(int ftrack_id){
    for (auto v:vertices){
        if ( v.IncludesTrack( ftrack_id ) ) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::AnalyzeVertices(){
    if (vertices.size()>0){
        for (auto & v:vertices){
            // after fixing the vertext position, remove far tracks
            Debug( 4 , "v.RemoveFarTracks( kMaxInterTrackDistance );");
            if (v.GetTracks().size()<2) continue;
            
            v.RemoveFarTracks( kMaxInterTrackDistance );
            
            // now sort the tracks
            Debug( 4 ,"SortTracksByPIDA;" );
            v.SortTracksByPIDA ();
            
            Debug( 4 , "SortTracksByLength;" );
            v.SortTracksByLength ();
            
            // and the relations between the tracks
            // inter-track distances, delta-theta, delta-phi...
            Debug( 4 , "SetTracksRelations;" );
            v.SetTracksRelations ();
            
            // flash - matching
            // find the closest flash to the vertex
            Debug( 4 , "vertices flash matching" );
            
            if (flashes.size()){
                float YZdistance = 1000, ClosestFlashYZdistance = 1000;
                for (auto f : flashes){
                    YZdistance = v.GetDis2Flash(f);
                    if ( YZdistance < ClosestFlashYZdistance ){
                        ClosestFlashYZdistance = YZdistance;
                        v.SetClosestFlash( f );
                    }
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::FilterGoodPairVertices(){

    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    
    // for momentum reconstuction from range / otherwise
    trkf::TrackMomentumCalculator trkm;

    // this funciton also kiils all vertices which are not up to standard.
    // to this end, we create a temporary vertices array and clear all vertices
    std::vector<pairVertex> tmp_vertices = vertices;
    vertices.clear();
    
    Debug(3 , "ub::ErezCCQEAnalyzer::FilterGoodPairVertices()");
    for (auto & v:tmp_vertices) {
        // vertices with only two tracks at close proximity and nothing else
        if (
            // we are looking for clusters of only two tracks that are fully contained
            v.GetNtracks() == 2
            // if there is a semi-contained track, with start/end point too close to the vertex, we don't want the vertex...
            &&  v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size() == 0
            ){
            
            // assign muon and proton tracks by PID-A
            auto AssignedMuonTrack = v.GetSmallPIDaTrack();
            auto PmuFromRange = trkm.GetTrackMomentum( AssignedMuonTrack.GetLength()  , 13  );
            v.AssignMuonTrack( AssignedMuonTrack  );
            
            
            auto AssignedProtonTrack = v.GetLargePIDaTrack();
            auto PpFromRange = trkm.GetTrackMomentum( AssignedProtonTrack.GetLength() , 2212);
            v.AssignProtonTrack( AssignedProtonTrack );
            
            v.FixTracksDirections ();
            v.SetReconstructedFeatures ( PmuFromRange , PpFromRange );
           
            // set vertex position in the three wire planes
            geo::TPCID tpcID = geom->FindTPCAtPosition( v.GetPosition() );
            if (tpcID.isValid) {
                int tpc = tpcID.TPC;
                for (int plane = 0 ; plane < 3 ; plane ++){
                    geo::PlaneID planeID = geo::PlaneID( 0 , tpc , plane ); // cryostat=0
                    float wire = geom->WireCoordinate( v.GetPosition().Y() , v.GetPosition().Z() ,  planeID );
                    float time = detprop->ConvertXToTicks( v.GetPosition().X() , planeID ) ;
                    // plug into the track
                    v.SetPlaneProjection( plane , wire , time );
                }
            }
            // plug hits associated with the proton and the muon tracks
            v.AssociateHitsToTracks( hits );
            
            vertices.push_back( v );
        }
        // if the vertex is not such, kill it
        else {
            Debug( 0 , Form("erasing vertex %d, since its not up to standards: N(tracks)=%d, N(close semi-contained tracks):%d"
                            , v.GetVertexID(), v.GetNtracks(), (int)(v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size()) ) );
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::TagVertices(){
    
    // tag vertices
    // ------------
    // µp           : event with a muon and a proton detected, and nothing else at close proximity
    // CC 1p 0π     : a subset of µp, in which only one proton with momentum > 200 MeV/c was produced, and no pions
    // non µp       : other pairs
    // cosmic       : pairs from cosmic (at least one of the tracks is cosmic, i.e. without MC information)
    //
    // Note that this functionallity is only relevant for MC events
    //
    for (auto & v:vertices){
        
        PandoraNuTrack t1 = v.GetAssignedMuonTrack();
        PandoraNuTrack t2 = v.GetAssignedProtonTrack();
        
        // for MC-BNB and cosmic-data overlay, we determine that a vertex is cosmic if one of its tracks is unrecognized (-9999)
        if ( (t1.GetMCpdgCode() * t2.GetMCpdgCode())==(13*2212) ){
            
            v.SetAs1mu1p();
            v.SetTrueMuonProton( t1 , t2 );
            
            // tag the vertex as true CC1p0π if it has two MC tracks, µ and p, which are
            // (1) close enough,
            // (2) come from a (close-enough) CC1p0π genie interaction
            if (
                (t1.GetTruthStartPos() - t2.GetTruthStartPos()).Mag() < 1.                  // distance between the true position of the two tracks is small
                && (v.GetClosestGENIE().GetVertexPosition() - v.GetPosition()).Mag() < 10   // distance from the closest genie vertex
                && (v.GetClosestGENIE().GetIsCC_1p_200MeVc_0pi()==true)                     // the closest GENIE is a CC1p0π
                ){
                v.SetAsCC1p0pi();
            }
        }
        else if ( t1.GetMCpdgCode()!=-9999 && t2.GetMCpdgCode()!=-9999 ){
            v.SetAsNon1mu1p();
        }
        else {
            v.SetAsCosmic();
        }
        
        // for MC-BNB and cosmic-MC overlay, we determine that a vertex is cosmic if one of its tracks is known to come from a cosmic ray
        if (t1.GetOrigin()=="cosmic ray" || t2.GetOrigin()=="cosmic ray"){
            v.SetAsCosmic();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::HeaderVerticesInCSV(){
    
    vertices_ctr = 0;
    
    vertices_file
    << "run" << "," << "subrun" << "," << "event" << "," << "vertex_id" << ","
    << "x" << "," << "y" << "," << "z" << ","
    << "track_id" << ","
    // tracks sorted by long / short
    << "PIDa_long"<< "," << "PIDa_short" << ","
    << "l_long"<< "," << "l_short" << ","
    // tracks sorted by small / large PIDa
    << "PIDa_small_PIDa"<< "," << "PIDa_large_PIDa" << ","
    << "l_small_PIDa"<< "," << "l_large_PIDa" << ","
    // µ/p assigned tracks
    << "PIDa_assigned_muon"<< "," << "PIDa_assigned_proton" << ","
    << "l_assigned_muon"<< "," << "l_assigned_proton" << ","
    // flash matching of tracks
    << "ClosestFlash_YZdistance_assigned_muon"<< "," << "ClosestFlash_YZdistance_assigned_proton" << ","
    << "ClosestFlash_TotalPE_assigned_muon"<< "," << "ClosestFlash_TotalPE_assigned_proton" << ","
    
    // start/end points, for FV cuts
    << "startx_assigned_muon" << ","
    << "starty_assigned_muon" << ","
    << "startz_assigned_muon" << ","
    << "startx_assigned_proton" << ","
    << "starty_assigned_proton" << ","
    << "startz_assigned_proton" << ","
    << "endx_assigned_muon" << ","
    << "endy_assigned_muon" << ","
    << "endz_assigned_muon" << ","
    << "endx_assigned_proton" << ","
    << "endy_assigned_proton" << ","
    << "endz_assigned_proton" << ","
    
    
    // CC1p0π reconstructed featues
    << "distance" << "," << "delta_phi" << "," << "delta_theta" << ","
    << "theta_12" << ","
    
    // reconstructed kinematics
    << "reco_Ev" << "," << "reco_Q2" << "," << "reco_Xb" << "," << "reco_y" << "," << "reco_W2" << ","
    << "reco_Pt" << "," << "reco_theta_pq" << ","
    << "reco_Pmu" << "," << "reco_Pmu_x" << "," << "reco_Pmu_y" << "," << "reco_Pmu_z" << "," << "reco_Pmu_theta" << "," << "reco_Pmu_phi" << ","
    << "reco_Pp" << "," << "reco_Pp_x" << "," << "reco_Pp_y" << "," << "reco_Pp_z" << "," << "reco_Pp_theta" << "," << "reco_Pp_phi" << ","

    // missing momentum (reconstructed struck neutron)
    << "reco_Pmiss" << ","  << "reco_Pmiss_x" << ","  << "reco_Pmiss_y" << ","  << "reco_Pmiss_z" << "," << "reco_Pmiss_t" << ","
    
    
    // truth MC information
    << "truth_l_assigned_muon"<< "," << "truth_l_assigned_proton" << ","
    
    << "truth_Pmu" << "," << "truth_Pmu_x" << "," << "truth_Pmu_y" << "," << "truth_Pmu_z" << "," << "truth_Pmu_theta" << ","<< "truth_Pmu_phi" << ","
    << "truth_Pp" << "," << "truth_Pp_x" << "," << "truth_Pp_y" << "," << "truth_Pp_z" << "," << "truth_Pp_theta" << ","<< "truth_Pp_phi" << ","
    
    // mathing genie interaction
    << "genie_distance" << ","
    << "truth_Ev" << "," << "truth_Q2" << "," << "truth_Xb" << "," << "truth_y" << "," << "truth_W2" << ","
    << "truth_Pt" << "," << "truth_theta_pq" << ","

    // closest genie (e.g. a proton was detected in a µp event, which is not the original proton in a CC interaction, since the real proton rescattered)
    << "closest_genie_distance" << ","
    << "closest_genie_Ev" << "," << "closest_genie_Q2" << "," << "closest_genie_Xb" << "," << "closest_genie_y" << "," << "closest_genie_W2" << ","
    << "closest_genie_Pt" << "," << "closest_genie_theta_pq" << ","
    
    // truth delta-phi
    << "truth_delta_phi" << ","
    // pdg code of long / short tracks and small / large PIDa
    << "pdg_long"<< "," << "pdg_short" << ","
    << "pdg_small_PIDa"<< "," << "pdg_large_PIDa" << ",";

    
    // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
    // see description of the observable at docdb-10958
    for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
        for (int plane = 0; plane < 3; plane++) {
            vertices_file << Form( "RdQaroundVertex[plane %d][%d wires x %d ticks]"
                                  , plane , NwiresBox[i_box_size] , NticksBox[i_box_size] ) << "," ;
        }
    }

    // flash matching of vertex
    vertices_file << "ClosestFlash_YZdistance" << "," << "ClosestFlash_TotalPE" << ",";

    // vertex truth-topology in MC
    vertices_file
    << "1mu-1p" << "," << "CC 1p 0pi" << "," << "other pairs" << "," << "cosmic"
    
    << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::StreamVerticesToCSV(){
    // July-25, 2017
    // whatever you add here - must add also in header - ub::ErezCCQEAnalyzer::HeaderVerticesInCSV()
    for (auto v:vertices){
        
        
        
        vertices_ctr++;
        
        vertices_file
        << v.GetRun() << "," << v.GetSubrun() << "," << v.GetEvent() << "," << v.GetVertexID() << ","
        << v.GetPosition().x() << "," << v.GetPosition().y() << "," << v.GetPosition().z() << ",";
        
        // tracks id
        for ( auto track : v.GetTracks()){
            vertices_file << track.GetTrackID() << ";" ;
        }
        vertices_file << ",";
        
        
        // tracks sorted by long / short
        vertices_file
        << v.GetLongestTrack().GetPIDa() << "," << v.GetShortestTrack().GetPIDa() << ","
        << v.GetLongestTrack().GetLength() << "," << v.GetShortestTrack().GetLength() << "," ;
        
        
        // tracks sorted by small / large PIDa
        vertices_file
        << v.GetSmallPIDaTrack().GetPIDa() << "," << v.GetLargePIDaTrack().GetPIDa() << ","
        << v.GetSmallPIDaTrack().GetLength() << "," << v.GetLargePIDaTrack().GetLength() << "," ;
        
        // µ/p assigned tracks
        vertices_file
        << v.GetAssignedMuonTrack().GetPIDa() << "," << v.GetAssignedProtonTrack().GetPIDa() << ","
        << v.GetAssignedMuonTrack().GetLength() << "," << v.GetAssignedProtonTrack().GetLength() << "," ;
        // flash matching of tracks
        vertices_file
        << v.GetAssignedMuonTrack().GetDis2ClosestFlash() << "," << v.GetAssignedProtonTrack().GetDis2ClosestFlash() << ","
        << v.GetAssignedMuonTrack().GetClosestFlash().GetTotalPE() << "," << v.GetAssignedProtonTrack().GetClosestFlash().GetTotalPE() << "," ;

        // start/end points, for FV cuts
        vertices_file
        << v.GetAssignedMuonTrack().GetStartPos().x() << ","
        << v.GetAssignedMuonTrack().GetStartPos().y() << ","
        << v.GetAssignedMuonTrack().GetStartPos().z() << ","
        << v.GetAssignedProtonTrack().GetStartPos().x() << ","
        << v.GetAssignedProtonTrack().GetStartPos().y() << ","
        << v.GetAssignedProtonTrack().GetStartPos().z() << ","
        << v.GetAssignedMuonTrack().GetEndPos().x() << ","
        << v.GetAssignedMuonTrack().GetEndPos().y() << ","
        << v.GetAssignedMuonTrack().GetEndPos().z() << ","
        << v.GetAssignedProtonTrack().GetEndPos().x() << ","
        << v.GetAssignedProtonTrack().GetEndPos().y() << ","
        << v.GetAssignedProtonTrack().GetEndPos().z() << ",";
        
        
        // CC1p0π reconstructed featues
        // distances between the tracks
        for ( auto d : v.Get_distances_ij()){
            vertices_file << d ; //  for more than 2 tracks add << ";" ;
        }
        vertices_file << ",";
        
        // delta_phi
        for ( auto delta_phi : v.Get_delta_phi_ij()){
            // in degrees
            vertices_file << delta_phi ; //  for more than 2 tracks add << ";" ;
        }
        vertices_file << ",";
        
        // delta_theta
        for ( auto delta_theta : v.Get_delta_theta_ij()){
            // in degrees
            vertices_file << delta_theta ; //  for more than 2 tracks add << ";" ;
        }
        vertices_file << ",";
        // theta_12 (the 3D angle between the two tracks)
        vertices_file << v.GetAngleBetween2tracks() << ","; // in degrees
        
        // reconstructed kinematics
        vertices_file << v.GetRecoEv() << "," << v.GetRecoQ2() << "," <<  v.GetRecoXb() << "," << v.GetRecoY() << "," << v.GetRecoW2() << ","  ;
        vertices_file << v.GetRecoPt() << "," << v.GetReco_theta_pq() << ",";
        vertices_file << v.GetRecoPmu().P() << "," << v.GetRecoPmu().Px() << "," << v.GetRecoPmu().Py() << "," << v.GetRecoPmu().Pz() << "," << v.GetRecoPmu().Theta() << "," << v.GetRecoPmu().Phi() << ",";
        vertices_file << v.GetRecoPp().P() << "," << v.GetRecoPp().Px() << "," << v.GetRecoPp().Py() << "," << v.GetRecoPp().Pz() << "," << v.GetRecoPp().Theta() << "," << v.GetRecoPp().Phi() << ",";

        // missing momentum (reconstructed struck neutron)
        vertices_file << v.GetRecoPmiss().P() << ","  << v.GetRecoPmiss().Px() << ","  << v.GetRecoPmiss().Py() << ","  << v.GetRecoPmiss().Pz() << "," << v.GetRecoPmiss().Pt() << ",";
        
        // truth MC information
        vertices_file << v.GetAssignedMuonTrack().GetTruthLength() << "," << v.GetAssignedProtonTrack().GetTruthLength() << ",";
        
        vertices_file
        << v.GetAssignedMuonTrack().GetTruthMomentum().P() << ","
        << v.GetAssignedMuonTrack().GetTruthMomentum().Px() << "," << v.GetAssignedMuonTrack().GetTruthMomentum().Py() << "," << v.GetAssignedMuonTrack().GetTruthMomentum().Pz() << "," << v.GetAssignedMuonTrack().GetTruthMomentum().Theta() << "," << v.GetAssignedMuonTrack().GetTruthMomentum().Phi() << ",";
        
        vertices_file
        << v.GetAssignedProtonTrack().GetTruthMomentum().P() << ","
        << v.GetAssignedProtonTrack().GetTruthMomentum().Px() << "," << v.GetAssignedProtonTrack().GetTruthMomentum().Py() << "," << v.GetAssignedProtonTrack().GetTruthMomentum().Pz() << "," << v.GetAssignedProtonTrack().GetTruthMomentum().Theta() << "," << v.GetAssignedProtonTrack().GetTruthMomentum().Phi() << ",";
        

        // mathing genie interaction
        vertices_file << v.GetDistanceToGENIE() << ",";
        vertices_file << v.GetGENIEinfo().GetEv() << "," << v.GetGENIEinfo().GetQ2() << "," << v.GetGENIEinfo().GetXb() << "," << v.GetGENIEinfo().GetY() << "," << v.GetGENIEinfo().GetW2() << ",";
        vertices_file << v.GetGENIEinfo().GetPt() << "," << v.GetGENIEinfo().Get_theta_pq() << ",";
        
        
        
        // closest genie (e.g. a proton was detected in a µp event, which is not the original proton in a CC interaction, since the real proton rescattered)
        vertices_file << v.GetDistanceToClosestGENIE() << ",";
        vertices_file << v.GetClosestGENIE().GetEv() << "," << v.GetClosestGENIE().GetQ2() << "," << v.GetClosestGENIE().GetXb() << "," << v.GetClosestGENIE().GetY() << "," << v.GetClosestGENIE().GetW2() << ",";
        vertices_file << v.GetClosestGENIE().GetPt() << "," << v.GetClosestGENIE().Get_theta_pq() << ",";
        
       
        // truth delta-phi
        vertices_file << v.GetTruthDeltaPhi() << ",";
        // pdg code of long / short tracks and small / large PIDa
        vertices_file << v.GetLongestTrack().GetMCpdgCode() << "," << v.GetShortestTrack().GetMCpdgCode() << ",";
        vertices_file << v.GetSmallPIDaTrack().GetMCpdgCode() << "," << v.GetLargePIDaTrack().GetMCpdgCode() << ",";

        
        
        // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
        // see description of the observable at docdb-10958
        for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
            for (int plane = 0; plane < 3; plane++) {
                vertices_file << v.GetRdQaroundVertex( plane, NwiresBox[i_box_size] , NticksBox[i_box_size] , hits ) << "," ;
            }
        }
        // flash matching of vertex
        vertices_file << v.GetDis2ClosestFlash() << "," << v.GetClosestFlash().GetTotalPE() << ",";

        
        // vertex truth-topology in MC
        vertices_file << v.GetIs1mu1p() << "," << v.GetIsGENIECC_1p_200MeVc_0pi() << "," << v.GetIsNon1mu1p() << "," << v.GetIsCosmic();
        
        // finish
        vertices_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::PrintInformation(){
    
    PrintXLine();
    Printf( "processed so far %d events", (int)(fTree->GetEntries()) );
    SHOW3( run , subrun , event );
    
    
    if (MCmode){
        if(!genie_interactions.empty()){
            cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << genie_interactions.size() << " genie interactions\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
            for (auto g: genie_interactions) {
                g.Print();
            }
        } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no interactions\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
    }
    
    if(!flashes.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << flashes.size() << " flashes\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto f: flashes) {
            f.Print();
        }
    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no flashes\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
   
    if(!tracks.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << tracks.size() << " pandoraNu tracks\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto t: tracks) {
            t.Print( true );
        }
    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no reco tracks\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
    
    
    if(!vertices.empty()){
        cout << "\033[36m" << "xxxxxxxxxxxxxx\n\n" << vertices.size() << " vertices\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto v: vertices) {
            v.Print( (debug>2) ? true : false       // do print pandoraNu tracks
                    );
        }
    } else {cout << "\033[36m" << "xxxxxxxxxxxxxx\n" << "no vertices\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}

    // time stamp
    PrintLine();
    end_ana_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
    std::cout << "\033[33m"
    << "finished analysis of this event at " << std::ctime(&end_time)
    << "time elapsed: " << elapsed_seconds.count() << "s"
    << "\033[31m" << endl;

    cout << "wrote " << vertices_ctr << " vertices to output file " << endl;
    EndEventBlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::beginJob(){
    
    // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
    for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
        NwiresBox[i_box_size] = MinNwiresBox + i_box_size * dNwiresBox;
        NticksBox[i_box_size] = MinNticksBox + i_box_size * dNticksBox;
    }

    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of events");
    fTree->Branch("run"     ,&run           ,"run/I");
    fTree->Branch("subrun"  ,&subrun        ,"subrun/I");
    fTree->Branch("event"   ,&event         ,"event/I");
    fTree->Branch("Ntracks" ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("Nhits"   ,&Nhits_stored  ,"Nhits/I");
    fTree->Branch("Nvertices",&Nvertices    ,"Nvertices/I");
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("hits"                ,&hits);
    fTree->Branch("genie_interactions"  ,&genie_interactions);
    fTree->Branch("vertices"            ,&vertices);

    
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("pot",&pot,"pot/D");
    fPOTTree->Branch("run",&run,"run/I");
    fPOTTree->Branch("subrun",&subrun,"subrun/I");

    // output csv file
    //    vertices_file.open("/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+fDataSampleLabel+"_vertices.csv");
    vertices_file.open(fDataSampleLabel+"_vertices.csv");
    cout << "opened vertices file: "+fDataSampleLabel+"_vertices.csv" << endl;
    HeaderVerticesInCSV();
    
    pot_total = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::endSubRun(const art::SubRun& sr){
    
    art::Handle< sumdata::POTSummary > potListHandle;
    
    if(
       ((fPOTModuleLabel=="generator") && (sr.getByLabel(fPOTModuleLabel,potListHandle))) // for MC-BNB
        ||
       ((fPOTModuleLabel=="beamdata") && (sr.getByLabel(fPOTModuleLabel,"bnbETOR860",potListHandle))) // for BNB, Marco Del Tutto, Aug-27, 2017
        )
        pot = potListHandle->totpot;
    else
        pot = 0.;
    if (fPOTTree) fPOTTree->Fill();
    pot_total += pot;
    
    if (debug){
        Printf("end subrun %d",subrun);
        SHOW2( pot , pot_total );
    }
    Printf( "POT from this subrun: %16.0lf" ,pot );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel        = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel    = p.get< std::string >("GenieGenModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
    fDataSampleLabel        = p.get< std::string >("DataSampleLabel");
    fPOTModuleLabel         = p.get< std::string >("POTModuleLabel");
    fFlashModuleLabel       = p.get< std::string >("FlashModuleLabel");

    debug = p.get< int >("VerbosityLevel");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ResetVars(){
    
    MCmode = false;
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = Nvertices = 0 ;
    pot = 0;
    tracks.clear();
    hits.clear();
    flashes.clear();
    genie_interactions.clear();
    vertices.clear();
    
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezCCQEAnalyzer)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -