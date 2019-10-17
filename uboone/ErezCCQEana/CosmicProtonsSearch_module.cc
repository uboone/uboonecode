////////////////////////////////////////////////////////////////////////
// Class:       CosmicProtonsSearch
// Plugin Type: analyzer (art v2_05_00)
// File:        CosmicProtonsSearch_module.cc
//
// Generated at 16-March 2018 Erez Cohen
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
#include "lardata/Utilities/PxUtils.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"

// for the new MC truth matching
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
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
#include <cstdarg>

// my pandoraNu track...
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/flash.h"
#include "uboone/ErezCCQEana/MyObjects/tripleVertex.h"

// constants
constexpr int debug          = 1;
constexpr int kMaxTrack      = 4000;  //maximum number of tracks
constexpr float kMaxInterTrackDistance = 11; // 11 cm between tracks - maximal distance for clustering
constexpr float kMin_dQ_inTruthMatchedHits = 0.1; // the minimal fraction of hits-charge for MC-truth matching
constexpr float EPSILON      = 0.1;   // tollerance for equations


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class CosmicProtonsSearch; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::CosmicProtonsSearch : public art::EDAnalyzer {
public:
    explicit CosmicProtonsSearch(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    CosmicProtonsSearch(CosmicProtonsSearch const &) = delete;
    CosmicProtonsSearch(CosmicProtonsSearch &&) = delete;
    CosmicProtonsSearch & operator = (CosmicProtonsSearch const &) = delete;
    CosmicProtonsSearch & operator = (CosmicProtonsSearch &&) = delete;
    
    // Required functions.
    void                   analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                  beginJob () override;
    void                    endJob () override;
    void               reconfigure (fhicl::ParameterSet const& p) override;
    
    void                 endSubRun (const art::SubRun& sr);

    
    // functionallity
    void         ConstructVertices ();
    void   ClusterTracksToVertices ();
    void           AnalyzeVertices ();
    void      FilterWantedVertices ();
    void          PrintInformation (bool DoPrintTracksFull=false);
    bool    TrackAlreadyInVertices (int ftrack_id);
    void       HeaderVerticesInCSV ();
    void       StreamVerticesToCSV ();

    bool ParticleAlreadyMatchedInThisHit ( std::vector<int> ,int );
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=0;
    void Debug(Int_t verobosity_level, const char* format) // base function
    {
        if ( debug < verobosity_level ) return;
        std::cout << format << std::endl;
    }
    template<typename T, typename... Targs>
    void Debug(Int_t verobosity_level, const char* format, T value, Targs... Fargs) // recursive variadic function
    {
        if ( debug < verobosity_level ) return;
        for ( ; *format != '\0'; format++ ) {
            if ( *format == '%' ) {
                std::cout << value << " ";
                Debug(verobosity_level, format+1, Fargs...); // recursive call
                return;
            }
            std::cout << *format;
        }
        std::cout << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    
    
    
private:
    
    // Declare member data here.
    void ResetVars();
    
    // Declare member data here.
    TTree *fTree;

    Int_t run, subrun, event;
    
    short   isdata;
    
    bool    MCmode=false;

    int     Ntracks;                // number of reconstructed tracks
    int     Nvertices;
    int     tracks_ctr, vertices_ctr;
    

    // my objects
    tripleVertex                    cvertex=tripleVertex();
    std::vector<PandoraNuTrack>     tracks;
    std::vector<tripleVertex>       vertices;
    
    
    //Module labels to get data products
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDataSampleLabel;
    std::string fPOTModuleLabel;
    std::string fHitParticleAssnsModuleLabel;
    std::string fG4ModuleLabel;

    //mctruth information
    Int_t    mcevts_truth;    //number of neutrino Int_teractions in the spill
    
    // time stamp
    std::chrono::time_point<std::chrono::system_clock> start_ana_time, end_ana_time;
    
    // output csv file of vertices
    ofstream vertices_file, summary_file;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::CosmicProtonsSearch::CosmicProtonsSearch(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::analyze(art::Event const & evt){
    
    ResetVars();
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    isdata = evt.isRealData();
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    
    // * tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
    
    // * associations
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);     // uncalibrated calorimetry
    
    // * MCTruth information
    auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, fTrackModuleLabel);
    
    art::Handle< std::vector<simb::MCParticle> >    pHandle;
    art::Handle< std::vector<simb::MCTruth> >       mctruthListHandle;
    art::Handle< std::vector<simb::MCFlux> >        mcfluxListHandle;
    std::vector< art::Ptr<simb::MCTruth> >          mclist;
    std::vector< art::Ptr<simb::MCFlux> >           fluxlist;
    
    if (MCmode) {
        evt.getByLabel(fG4ModuleLabel, pHandle);
        art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
    }
    
    // check if this makes data runs crash. If so, change to: fMCmodeLabel != "MC"
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
    
    
    

    // ----------------------------------------
    // tracks information
    // ----------------------------------------
    Ntracks = tracklist.size();
    tracks_ctr += Ntracks;
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
                        if (debug>3) SHOW3(track.GetTrackID(),plane,calo->dEdx().size());
                        
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
                            Debug(3 , "pida: %",pida);
                        }
                        // record the track dE/dx vs. residual range
                        track.FillResRange( plane , calo->ResidualRange()[ip] );
                        track.FilldEdxHit( plane , calo->dEdx()[ip] );
                    }
                    if (used_trkres) pida /= used_trkres;
                    track.SetPIDaPerPlane( plane , pida );
                }
            }
            track.SetPIDa();
            Debug(3 , "track.GetPIDa(): %",track.GetPIDa());
        }
        
        
        
        // MC-truth mathching for the tracks
        Debug(4,"before if (MCmode=%)",MCmode);
        if (MCmode){
            
            Debug(3,"inside (MCmode)");

            evt.getByLabel(fG4ModuleLabel, pHandle);
            art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);

            
            std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i);
            Debug(4,"\tThere are % associated hits to track %." ,(int)trk_hits_ptrs.size(), track.GetTrackID());
            
            std::unordered_map<int,double> trkide;
            double maxe=-1, tote=0, purity=0;
            art::Ptr< simb::MCParticle > maxp_me; //pointer for the particle match we will calculate
            
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle , evt , fHitParticleAssnsModuleLabel);
            std::vector<art::Ptr<simb::MCParticle>> particle_vec;
            std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
            
            //loop only over our hits
            std::unordered_map<int,double> trackid_dQinTruthMatchedHits;
            double max_dQinTruthMatchedHits=-1, dQinAllHits=0;

            for(size_t i_h=0; i_h < trk_hits_ptrs.size(); ++i_h){
                Debug(4,"--------------------------------");
                float dQinHit = trk_hits_ptrs[i_h]->Integral();
                dQinAllHits += dQinHit;
                Debug(4,"in track %, hit %, dQinAllHits: %",track.GetTrackID(),(int)i_h,dQinAllHits);
                // for each hit, ask how many particles match this hit
                particle_vec.clear(); match_vec.clear();
                particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection
                Debug(4,"in track %, There are % particles matched to hit %" ,track.GetTrackID(), (int)particle_vec.size() ,(int)i_h );
                
                // to avoid from matching the same particle more than once
                // we introduce a vector of matched TrackId-s
                // and require that the matched particle has not been mathced already for this hit
                std::vector <int> ParticlesMatchedInThisHit;
                if (particle_vec.size()>0){
                    //loop over particles that match this hit and ask which one deposited the most energy
                    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                        Debug(5,"i_p: %, particle_vec[i_p]->TrackId(): %...",(int)i_p,(int)particle_vec[i_p]->TrackId());
                        float Edep_particle = match_vec[i_p]->energy;  // energy deposited by ionization by this track ID [MeV]
                        float Edep_particle_frac = match_vec[i_p]->ideFraction; // fraction of energy in hit from this particle
                        Debug(5, "Edep_particle:%, Edep_particle_frac:%",Edep_particle,Edep_particle_frac);
                        
                        trkide[ particle_vec[i_p]->TrackId() ] += Edep_particle; //store energy [MeV] deposited by track id

                        tote += Edep_particle; //calculate total energy deposited
                        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                            maxe = trkide[ particle_vec[i_p]->TrackId() ];
                            maxp_me = particle_vec[i_p];
                            Debug(5,"changing maxp_me to particle_vec[%], pdg:%, deposited maxe:%f",(int)i_p,maxp_me->PdgCode(),maxe);
                            
                            if (!ParticleAlreadyMatchedInThisHit( ParticlesMatchedInThisHit , (int)particle_vec[i_p]->TrackId()) ) {
                                ParticlesMatchedInThisHit.push_back( (int)particle_vec[i_p]->TrackId() );
                                trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ] += dQinHit; // store the integral on the hit by the track id
                                Debug(5, "dQinHit:%, trackid_dQinTruthMatchedHits[%]:%",dQinHit,particle_vec[i_p]->TrackId(),trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ]);
                                max_dQinTruthMatchedHits = trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ];
                                Debug(5,"max_dQinTruthMatchedHits:%, dQinAllHits:%",max_dQinTruthMatchedHits,dQinAllHits);
                            } else {
                                Debug(5,"particle of TrackID % was already matched in this hit",(int)particle_vec[i_p]->TrackId());
                            }
                        }
                    }//end loop over particles per hit
                }
                Debug(5,"after if (particle_vec.size()>0); maxe=%",maxe);
                purity = (fabs(tote)>0) ? maxe/tote : 0.0;
            }
            
            const art::Ptr< simb::MCParticle > particle = maxp_me;
            // Now have matched the truth information - plug into the track object
            // first, check completeness of the track MC-truth matching:
            // if a particle is to be matched with a track, require that it deposited more than some fraction kMin_dQ_inTruthMatchedHits in the charge deposited in the track
            // keep track of the ratio max_dQinTruthMatchedHits / dQinAllHits,
            // in order to select the best value to cut on...
            Debug(5,"track.SetMaxdQinTruthMatchedHits(%)",max_dQinTruthMatchedHits);
            track.SetMaxdQinTruthMatchedHits(max_dQinTruthMatchedHits);
            Debug(5,"track.SetdQinAllHits(%)",dQinAllHits);
            track.SetdQinAllHits(dQinAllHits);
            if (debug>5) SHOW3( max_dQinTruthMatchedHits , dQinAllHits, track.GetRatiodQinTruthMatchedHits() );
            if ( (!particle.isNull()) && (track.GetRatiodQinTruthMatchedHits() >= kMin_dQ_inTruthMatchedHits )) {
                track.SetMCpdgCode( particle->PdgCode() );
                track.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                track.SetTruthDirection();
                track.SetTruthLength();
                track.SetTruthMomentum( particle -> Momentum() );
                track.SetTruthMother( particle -> Mother() );
                track.SetTruthProcess( particle -> Process() );
                track.SetTruthPurity( purity );
                // * MC-truth information related to the associated particle
                int particle_key = (int)particle.key();
                if( fo.isValid() ){
                    auto pHandle_at_particle_key = pHandle->at(particle_key);
                    art::Ptr<simb::MCTruth> mc_truth = fo.at(particle_key);
                    track.SetOrigin( (mc_truth->Origin() == simb::kBeamNeutrino)
                                    ? "beam neutrino"
                                    : ((mc_truth->Origin() == simb::kCosmicRay)
                                       ? "cosmic ray"
                                       : "unkwon origin")  );
                }// end if fo.isValid()
            }else {
                Debug(3,"particle.isNull() = true!...");
                if (!maxp_me.isNull()) {
                    Debug(5,  "Un-Matching particle % since it deposited too few of hit-charge in track (max_dQinTruthMatchedHits/dQinAllHits = %)" , maxp_me->PdgCode() , max_dQinTruthMatchedHits / dQinAllHits );
                }
            }
        }//MC
        tracks.push_back( track );
    }
    
    
    ConstructVertices();
    fTree -> Fill();
    // ----------------------------------------
    // print and finish
    // ----------------------------------------
    PrintInformation( (debug>0) ? true : false );
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::CosmicProtonsSearch::ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs
                                         ,int cTrackID ){
    // to avoid from matching the same particle more than once
    // we introduce a vector of matched TrackId-s for each hit
    // and require that the matched particle has not been mathced already for this hit
    if (debug>5) {
        cout << "ub::CosmicProtonsSearch::ParticleAlreadyMatchedInThisHit()" << endl;
        cout << "looking if track " <<  cTrackID << " has already matched in the list: " << endl;
        for (auto trk_id:AlreadyMatched_TrackIDs) {
            cout << trk_id << ",";
        }
        cout << endl;
    }
    for (auto trk_id:AlreadyMatched_TrackIDs) {
        if (trk_id == cTrackID) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::ConstructVertices(){
    
    // cluster all tracks at close proximity to vertices
    ClusterTracksToVertices();
    // analyze these vertices: inter-tracks distances, angles...
    AnalyzeVertices();
    // retain only vertices with pairs of 3-tracks at close proximity
    FilterWantedVertices();
    StreamVerticesToCSV();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::ClusterTracksToVertices(){
    
    // Jan-18, 2018
    // cluster all tracks at close proximity to vertices, and fix the position of each vertex
    bool    FoundCloseTracks=false;
    float   closest_distance_ij;
    TVector3 vertex_position;
    
    for (int i=0; i < Ntracks; i++){
        if (!tracks[i].IsTrackContainedSoft()) continue; // 3<x<257, -115<y<115, 5<z<1045
        tripleVertex vertex( run, subrun, event , vertices.size() );
        vertex.AddTrack( tracks[i] );
        FoundCloseTracks = false;
        for ( int j=i+1 ; j < Ntracks ; j++ ){ // for each track i, loop over all other tracks j!=i
            if (tracks[j].IsTrackContainedSoft() && j!=i){ // 3<x<257, -115<y<115, 5<z<1045
                // if this is the first time we go over these two tracks
                // and they are close enough to define a vertex,
                // we also define the position of their mutual vertex
                if (!FoundCloseTracks){
                    std::string StartOrEnd = "None";
                    closest_distance_ij = tracks[i].ClosestDistanceToOtherTrack(tracks[j],&StartOrEnd);
                    if ( closest_distance_ij < kMaxInterTrackDistance ){ // if track-i and track-j are close enough we want to keep this vertex
                        vertex.AddTrack( tracks[j] );
                        FoundCloseTracks = true;
                        vertex.SetPosition( (StartOrEnd.compare("Start")==0)
                                           ?
                                           tracks[i].GetStartPos()
                                           :
                                           tracks[i].GetEndPos());
                    }// end if keep vertex
                }// end if !FoundCloseTracks
                else { // if we have already clustered a vertex, and positioned it in space, we only need to check wether the new track (j) is close enough to this vertex
                    if ( tracks[j].DistanceFromPoint(vertex.GetPosition()) < kMaxInterTrackDistance ){
                        vertex.AddTrack( tracks[j] );
                    }
                }
            }
        }// end loop on track-j
        // plug into list of vertices
        if (FoundCloseTracks) {
            vertices.push_back( vertex );
        }
        Debug( 4 , "} if (FoundCloseTracks)");
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::CosmicProtonsSearch::TrackAlreadyInVertices(int ftrack_id){
    for (auto v:vertices){
        if ( v.IncludesTrack( ftrack_id ) ) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::AnalyzeVertices(){
    if (vertices.size()>0){
        for (auto & v:vertices){
            if (v.GetTracks().size()<2) continue;
            // here we can do some analysis...
            Debug( 4 , "SortTracksByLength;" );
            v.SortTracksByLength ();
            v.SetTracksRelations ();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::FilterWantedVertices(){
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    
    // this funciton also kiils all vertices which are not up to standard.
    // to this end, we create a temporary vertices array and clear all vertices
    std::vector<tripleVertex> tmp_vertices = vertices;
    vertices.clear();
    
    Debug(3 , "ub::CosmicProtonsSearch::FilterWantedVertices()");
    for (auto & v:tmp_vertices) {
        if (
            // we are looking for clusters of only three tracks,
            // we can later ask that they were contained or at least the proton contained,
            // or that there are no near-by tracks in the surroundings
            v.GetNtracks() == 3
            ){
            
            // set vertex position in the three wire planes
            // version for v06_26_01
            double const v_position[3] = {v.GetPosition().X(),v.GetPosition().Y(),v.GetPosition().Z()};
            geo::TPCID tpcID = geom->FindTPCAtPosition( v_position );
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
            vertices.push_back( v );
        }
        // if the vertex is not like a topology we wanted...
        else { Debug( 0 , "erasing vertex %, since N(tracks)=%", v.GetVertexID(), v.GetNtracks()) ;}
    }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::HeaderVerticesInCSV(){
    // March-16, 2018
    vertices_ctr = 0;
    
    vertices_file
    << "run" << "," << "subrun" << "," << "event" << "," << "vertex_id" << ","
    << "x" << "," << "y" << "," << "z" << ","
    << "TrackIDs" << ","
    
    // truth MC information
    << "truth_pdg_1" << ","
    << "truth_pdg_2" << ","
    << "truth_pdg_3"
    
    // finish
    << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::StreamVerticesToCSV(){
    // March-16, 2018
    // whatever you add here - must add also in header - ub::CosmicProtonsSearch::HeaderVerticesInCSV()
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
        
        
        // truth MC information
        vertices_file
        << v.GetTracks().at(0).GetMCpdgCode() << ","
        << v.GetTracks().at(1).GetMCpdgCode() << ","
        << v.GetTracks().at(2).GetMCpdgCode() ;
        
        // finish
        vertices_file << endl;
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::PrintInformation(bool DoPrintTracksFull){
    
    PrintXLine();
    SHOW3( run , subrun , event );
    
    
    
    if(!tracks.empty()){ PrintHeader(Form("%d pandoraNu-tracks",(int)tracks.size()));
        for (auto t: tracks) {
            if (DoPrintTracksFull) t.Print( true );
            else cout << t.GetTrackID() << "\t";
        }
        cout << endl;
    } else { PrintHeader("no pandoraNu-tracks");    }
    
    if(!vertices.empty()){
        PrintHeader(Form("%d triple-vertices",(int)vertices.size()));
        for (auto v: vertices) {
            v.Print(false,false);
            if (v.GetTracks().at(0).GetMCpdgCode()*v.GetTracks().at(1).GetMCpdgCode()*v.GetTracks().at(2).GetMCpdgCode()
                ==
                13*13*2212) {
                Printf("This is a vertex of two muons and a proton!");
            }
        }
    } else { PrintHeader("no triple-vertices"); }
    
    
    // time stamp
    if (debug>2){
        PrintLine();
        end_ana_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
        std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
        std::cout << "\033[33m"
        << "finished analysis of this event at " << std::ctime(&end_time)
        << "time elapsed: " << elapsed_seconds.count() << "s"
        << "\033[31m" << endl;
    }
    PrintLine();
    std::cout << "\033[33m"
    << "processed so far " << (int)(fTree->GetEntries()) << " events" << endl;
    cout << "wrote " << vertices_ctr << " vertices to output file " << "\033[31m"  << endl;
    EndEventBlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::beginJob(){
    
    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of events");
    fTree->Branch("run"     ,&run           ,"run/I");
    fTree->Branch("subrun"  ,&subrun        ,"subrun/I");
    fTree->Branch("event"   ,&event         ,"event/I");
    fTree->Branch("Ntracks" ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("Nvertices",&Nvertices    ,"Nvertices/I");
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("vertices"            ,&vertices);
    
    // output csv files
    vertices_file.open(fDataSampleLabel+"_vertices.csv");
    cout << "opened vertices file: "+fDataSampleLabel+"_vertices.csv" << endl;
    HeaderVerticesInCSV();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::endJob(){
    Debug(3,"ub::CosmicProtonsSearch::endJob()");
    
    std::time_t now_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    summary_file.open(fDataSampleLabel+"_summary.csv");
    summary_file << "time" << ","
    << "Nevents" << ","
    << "Ntracks" << ","
    << "Nvertices"
    << endl;
    
    std::string sTimeS = std::ctime(&now_time);
    
    summary_file << sTimeS.substr(0,sTimeS.length()-1) << ","
    << fTree->GetEntries() << ","
    << tracks_ctr << ","
    << vertices_ctr
    << endl;
    
    summary_file.close();
    
    cout << "Ended job. Wrote summary into file: "+fDataSampleLabel+"_summary.csv" << endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::endSubRun(const art::SubRun& sr){
    
    Debug(0,"end subrun %",subrun);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel        = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel    = p.get< std::string >("GenieGenModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
    
    fDataSampleLabel        = p.get< std::string >("DataSampleLabel");
    fPOTModuleLabel         = p.get< std::string >("POTModuleLabel");
    fFlashModuleLabel       = p.get< std::string >("FlashModuleLabel");
    debug                   = p.get< int >("VerbosityLevel");
    MCmode                  = p.get< bool >("MCmodeLabel",false);
    fHitParticleAssnsModuleLabel = p.get< std::string >("HitParticleAssnsModuleLabel");
    fG4ModuleLabel          = p.get< std::string >("G4ModuleLabel","largeant");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicProtonsSearch::ResetVars(){
    
    run = subrun = event = -9999;
    Ntracks = Nvertices = 0 ;
    tracks.clear();
    vertices.clear();
    
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::CosmicProtonsSearch)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
