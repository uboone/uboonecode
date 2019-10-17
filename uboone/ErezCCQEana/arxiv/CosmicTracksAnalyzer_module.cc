////////////////////////////////////////////////////////////////////////
// Class:       CosmicTracksAnalyzer
// Plugin Type: analyzer (art v2_05_00)
// File:        CosmicTracksAnalyzer_module.cc
//
// Generated at Jan 14 2018 by Erez Cohen using cetskelgen
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"
// for the new MC truth matching (by Wes)
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
// for MCtrack
#include "lardataobj/MCBase/MCTrack.h"
// lar core algorithm
#include "uboone/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoLineSegment.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoLine.h"

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
constexpr int kMaxTruth      = 100000;
constexpr int kMaxNgenie     = 100;
// 11 cm between tracks - maximal distance for clustering, here we take 30 cm to open the possible clustering window
constexpr float kMaxInterTrackDistance = 30;
constexpr float EPSILON      = 0.1;   // tollerance for equations

// charge deposition around the vertex in a box of N(wires) x N(time-ticks)
constexpr int N_box_sizes    = 30;
constexpr int MinNwiresBox   = 5;
constexpr int dNwiresBox     = 5;
constexpr int MinNticksBox   = 10;
constexpr int dNticksBox     = 10;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class CosmicTracksAnalyzer; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::CosmicTracksAnalyzer : public art::EDAnalyzer {
public:
    explicit CosmicTracksAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    CosmicTracksAnalyzer(CosmicTracksAnalyzer const &) = delete;
    CosmicTracksAnalyzer(CosmicTracksAnalyzer &&) = delete;
    CosmicTracksAnalyzer & operator = (CosmicTracksAnalyzer const &) = delete;
    CosmicTracksAnalyzer & operator = (CosmicTracksAnalyzer &&) = delete;
    
    // Required functions.
    void                   analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                  beginJob () override;
    void                    endJob () override;
    void               reconfigure (fhicl::ParameterSet const& p) override;
    
    void                 endSubRun (const art::SubRun& sr);

    
    // functionallity
    void     FindTwoTracksVertices (std::vector<PandoraNuTrack> & tracks_vector // tracks to be clustered
                                     ,std::vector<pairVertex> & vertices_vector  // vertices vector to hold the pairs
                                    );
    void             FindTrajPairs (std::vector<PandoraNuTrack> & tracks_vector // tracks to be clustered
    );
    void           AnalyzeVertices (std::vector<pairVertex> & vertices_vector  // vertices vector to hold the pairs
                                    );
    void    FilterGoodPairVertices (std::vector<pairVertex> & vertices_vector );
    void        FilterGoodTrajPair ();

    void          PrintInformation ();
    void       HeaderVerticesInCSV ();
    void       StreamVerticesToCSV ();
    void HeaderTrajectoryPairsInCSV();
    void StreamTrajectoryPairsToCSV();
    void         MatchPairVertices ();
    void      CheckIfReconstructed ( std::vector<pairVertex> & truth_pairs_vector );
    void    CheckIfReconstructable ( std::vector<pairVertex> & truth_pairs_vector );
    bool IsTrajecPairReconstructed ( PandoraNuTrack truth_trajectory );
    bool HasMatchingRecoPairVertex ( pairVertex truth_pair );
    pairVertex  FormPairFrom2Trajs (PandoraNuTrack ti,PandoraNuTrack tj, TVector3 position);
    void     SearchOtherCloseTrajs (pairVertex & pair , int i , int j);

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
    TTree* fPOTTree;

    Int_t   run, subrun, event;
    
    short   isdata;
    
    bool    MCmode;
    bool    DoPrintTruthTrajectories = true;
    bool    DoPrintCosmicTracks = false;
    bool    DoPrintPandoraNuTracks = false;
    int     NCosmicTracks_total=0, Ntracks_total=0;
    int     truth_pairs_ctr=0, truth_pairs_ctr_reconstructable=0;
    int     truth_pairs_ctr_reconstructed=0, truth_pairs_ctr_contained=0;
    int     NCosmicTracks, Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    int     Nflashes;
    int     vertices_counter;
    int     NwiresBox[N_box_sizes], NticksBox[N_box_sizes];
    

    long    Nhits_tot=0;
    long    NMCParticles_total=0;

    double  pot, pot_total;

    // my objects
    std::vector<PandoraNuTrack>     tracks, cosmic_tracks;
    std::vector<PandoraNuTrack>     truth_trajectories; // trajectories of truth particles
    std::vector<hit>                hits;
    std::vector<pairVertex>         vertices, cosmic_vertices;
    std::vector<pairVertex>         truth_pairs;    // pairs of trajectories of truth particles
    std::vector<flash>              flashes;
    
    
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
    std::string fMCTrackModuleLabel;
    std::string fCosmicTrackModuleLabel;
    
    //mctruth information
    Int_t    mcevts_truth;    //number of neutrino Int_teractions in the spill
    
    // time stamp
    std::chrono::time_point<std::chrono::system_clock> start_ana_time, end_ana_time;
    
    // output csv file of vertices
    ofstream vertices_file, summary_file, truth_pairs_file;
    geoalgo::GeoAlgo mygeoalgo;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::CosmicTracksAnalyzer::CosmicTracksAnalyzer(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::analyze(art::Event const & evt){
    
    ResetVars();
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    if (debug>0) SHOW3(run , subrun , event);
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // * run data
    isdata = evt.isRealData();
    
    // * hits
    Debug(2,"// * hits");
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
    
    // * flashes
    Debug(2,"// * flashes");
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fFlashModuleLabel, flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);

    // * MCtracks
    Debug(2,"// * MCtracks");
    art::Handle< std::vector<sim::MCTrack> > MCtrackListHandle;
    std::vector<art::Ptr<sim::MCTrack> > MCtracklist;
    if (evt.getByLabel(fMCTrackModuleLabel,MCtrackListHandle))
        art::fill_ptr_vector(MCtracklist, MCtrackListHandle);
    
    // * tracks
    Debug(2,"// * tracks");
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

    // * cosmic-tracks
    Debug(2,"// * cosmic-tracks");
    art::Handle< std::vector<recob::Track> > CosmicTrackListHandle;
    std::vector<art::Ptr<recob::Track> > CosmicTracklist;
    if (evt.getByLabel(fCosmicTrackModuleLabel,CosmicTrackListHandle))
        art::fill_ptr_vector(CosmicTracklist, CosmicTrackListHandle);

    // * associations
    Debug(2,"// * associations");
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

    
    // * new truth matching from /uboone/app/users/wketchum/dev_areas/mcc8_4_drop/gallery_macros/TruthMatchTracks.C
    auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, fTrackModuleLabel);
    
    // * MCTruth information
    // get the particles from the event
    Debug(2,"// * MCTruth information");
    art::Handle<std::vector<simb::MCParticle>> pHandle;
    evt.getByLabel(fG4ModuleLabel, pHandle);
    art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
    
    Debug(2,"// * MCTruth information");
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
    
    // MCParticle from largeant
    Debug(2,"// * MCParticle information");
    art::Handle< std::vector<simb::MCParticle> > MCParticleListHandle;
    std::vector<art::Ptr<simb::MCParticle> > mcparticlelist;
    if (evt.getByLabel("largeant",MCParticleListHandle))
        art::fill_ptr_vector(mcparticlelist, MCParticleListHandle);
    auto Nparticles = (int)mcparticlelist.size();
    Debug(1,"Nparticles: %",Nparticles);

    // MCFlux
    Debug(2,"// * MCFlux information");
    art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
    std::vector<art::Ptr<simb::MCFlux> > fluxlist;
    if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
        art::fill_ptr_vector(fluxlist, mcfluxListHandle);

    
    // ----------------------------------------
    // hits information
    // ----------------------------------------
    Nhits = hitlist.size();
    Nhits_tot += Nhits;
    Nhits_stored = std::min(Nhits, kMaxHits);
    Debug(2,"Nhits: %, Nhits_stored: %",Nhits,Nhits_stored);
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
    Debug(4,"Nflashes: %",Nflashes);
    for ( int f = 0; f < std::min(Nflashes,kMaxHits); f++ ) {
        
        if (debug>5){ SHOW2(flashlist[f]->Time() , flashlist[f]->TotalPE() ) };
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
        if (debug>5){ SHOW2(fflash.GetTime() , fflash.GetTotalPE()) };
        if( (0.0 < fflash.GetTime()) && (fflash.GetTime() < 10.0) && (6.5 < fflash.GetTotalPE()) ){
            flashes.push_back( fflash );
        }
    }


    // ----------------------------------------
    // cosmic-tracks information
    // ----------------------------------------
    NCosmicTracks = CosmicTracklist.size();
    NCosmicTracks_total += NCosmicTracks;
    Debug(2,"NCosmicTracks: %, NCosmicTracks_total: %",NCosmicTracks,NCosmicTracks_total);
    for(int i=0; i < std::min(int(CosmicTracklist.size()),kMaxTrack); ++i ){
        recob::Track::Point_t start_pos, end_pos;
        std::tie( start_pos, end_pos ) = CosmicTracklist[i]->Extent();
        
        Debug(5 , Form("analyzing cosmic track %d in this event",CosmicTracklist[i]->ID()) );
        PandoraNuTrack cosmic_track(
                                    run , subrun, event        // r/s/e
                                    ,CosmicTracklist[i]->ID()        // track id
                                    ,CosmicTracklist[i]->Length()    // length
                                    ,TVector3(start_pos.X(),start_pos.Y(),start_pos.Z())   // start position
                                    ,TVector3(end_pos.X(),end_pos.Y(),end_pos.Z())         // end position
                                    );
        
        cosmic_tracks.push_back( cosmic_track );
    } // end loop of cosmic tracks
    FindTwoTracksVertices( cosmic_tracks , cosmic_vertices );
    AnalyzeVertices( cosmic_vertices );
    FilterGoodPairVertices( cosmic_vertices );
    
    // ----------------------------------------
    // tracks information
    // ----------------------------------------
    Ntracks = tracklist.size();
    Ntracks_total += Ntracks;
    Debug(2,"Ntracks: %, Ntracks_total: %",Ntracks,Ntracks_total);
    Debug(2,"std::min(int(tracklist.size()),kMaxTrack):%",std::min(int(tracklist.size()),kMaxTrack));
    for(int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i ){
        Debug(2, "analyzing track %",i);
        if (debug>2) SHOW(tracklist[i]->ID());
        //        Debug(3 , "analyzing track tracklist[%]->ID()=%",i,tracklist[i]->ID() );
        recob::Track::Point_t start_pos, end_pos;
        std::tie( start_pos, end_pos ) = tracklist[i]->Extent();
        
        Debug(2,"PandoraNuTrack track");
        PandoraNuTrack track(
                             run , subrun, event        // r/s/e
                             ,tracklist[i]->ID()        // track id
                             ,tracklist[i]->Length()    // length
                             ,TVector3(start_pos.X(),start_pos.Y(),start_pos.Z())   // start position
                             ,TVector3(end_pos.X(),end_pos.Y(),end_pos.Z())         // end position
                            );
        
        double StartLoc[3] = {start_pos.X(), start_pos.Y(), start_pos.Z()};
        
        // U / V / Y coordinates
        Debug(3,"// U / V / Y coordinates");
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
        Debug(2,"// PIDa and calorimetric KE");
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
        
        // flash - matching
        // find the closest flash to the track
        Debug(2,"// flash - matching");
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
        
        // MC information
        Debug(2,"// MC information");
        bool DoNewMCtruthMatching = true;
        if (mclist.size() && DoNewMCtruthMatching){

            std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i);
            std::unordered_map<int,double> trkide;
            double maxe=-1, tote=0;
            art::Ptr< simb::MCParticle > maxp_me; //pointer for the particle match we will calculate
            
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle , evt , fHitParticleAssnsModuleLabel);
            Debug(2,"art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle , evt , fHitParticleAssnsModuleLabel);");
            std::vector<art::Ptr<simb::MCParticle>> particle_vec;
            Debug(2,"std::vector<art::Ptr<simb::MCParticle>> particle_vec;");
            std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
            Debug(2,"std::vector<anab::BackTrackerHitMatchingData const*> match_vec;");
            
            //loop only over our hits
            for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h){
                
                // for each hit, ask how many particles match this hit
                Debug(4,Form("i_h: %d, particle_vec.clear(); match_vec.clear();",(int)i_h));
                particle_vec.clear(); match_vec.clear();
                particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection
                Debug(4,Form("There are %d particles matched to hit %d" , (int)particle_vec.size() ,(int)i_h ));
                
                if (particle_vec.size()>0){
                    Debug(5,Form("%d",(int)particle_vec.size()));
                    //loop over particles
                    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                        Debug(5,Form("i_p: %d, particle_vec[i_p]->TrackId(): %d...",(int)i_p,(int)particle_vec[i_p]->TrackId()));
                        trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
                        tote += match_vec[i_p]->energy; //calculate total energy deposited
                        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                            maxe = trkide[ particle_vec[i_p]->TrackId() ];
                            maxp_me = particle_vec[i_p];
                        }
                    }//end loop over particles per hit
                }
                Debug(4,"after if (particle_vec.size()>0)");
            }
            Debug(2,"after for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h) loop");
            // Now have matched the truth information - plug into the track object
            const art::Ptr< simb::MCParticle > particle = maxp_me;
            Debug(2 , "const art::Ptr< simb::MCParticle > particle = maxp_me;" );
            if (!particle.isNull()){
                Debug(2,"particle.isNull() = false!...");

                track.SetMCpdgCode( particle->PdgCode() );
                track.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                track.SetTruthDirection();
                track.SetTruthLength();
                track.SetTruthMomentum( particle -> Momentum() );
                track.SetTruthMother( particle -> Mother() );
                track.SetTruthProcess( particle -> Process() );
                
                // from which particle did it come in the largeant list
                for( int i_mcparticle = 0; (i_mcparticle < std::min(Nparticles,kMaxTruth)) ; i_mcparticle++ ){
                    art::Ptr<simb::MCParticle> mcparticle = mcparticlelist.at(i_mcparticle);
                    if (mcparticle == particle) {
                        track.SetTruthMCParticle( i_mcparticle );
                    }
                }
                // art::Ptr< simb::MCParticle >  track_truth_particle = particle;
                
                // * MC-truth information of this particle
                int particle_key = (int)particle.key();
                Debug(2,"before if( fo.isValid() )");
                if( fo.isValid() ){
                    auto pHandle_at_particle_key = pHandle->at(particle_key);
                    art::Ptr<simb::MCTruth> mc_truth = fo.at(particle_key);
                    if (mc_truth->Origin() == simb::kBeamNeutrino)      track.SetOrigin( "beam neutrino" );
                    else if (mc_truth->Origin() == simb::kCosmicRay)    track.SetOrigin( "cosmic ray" );
                    if (debug>5) {SHOW( mc_truth -> Origin() );}
                }// end if fo.isValid()
                Debug(2,"end if fo.isValid()");
            } else {
                Debug(2,"particle.isNull() = true!...");
            }
        }//MC

        tracks.push_back( track );
        Debug(3 , "done analyzing track tracklist[%]->ID()=%",i,tracklist[i]->ID() );
    }// for loop on tracks: for(int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i )
    FindTwoTracksVertices( tracks , vertices );
    AnalyzeVertices( vertices );
    FilterGoodPairVertices( vertices );

    // check if the pandoraCosmic vertices match the pandoraNu vertices
    MatchPairVertices();
    Debug(2,"passed MatchPairVertices()");
    // this is seemingly irrelevant since the pandoraNu tracks != pandoraCosmic tracks
    
    
    StreamVerticesToCSV();
    Debug(2,"StreamVerticesToCSV()");

    
    
    Debug(2,"before if (mclist.size()=%)",mcparticlelist.size());
    // record all start- and end-points of the truth-particles with path-length longer than the 3-wires requirement (1 cm)
    int Nparticles_passed_cut=0;
    if (Nparticles){
        MCmode = true;
        for( int i_mcparticle = 0; (i_mcparticle < std::min(Nparticles,kMaxTruth)) ; i_mcparticle++ ){
            art::Ptr<simb::MCParticle> particle = mcparticlelist.at(i_mcparticle);

            if (particle->StatusCode()!=1) continue;
            TVector3 particle_start(particle->Vx() , particle->Vy() , particle->Vz());
            TVector3 particle_end(particle->EndX() , particle->EndY() , particle->EndZ());
            PandoraNuTrack truth_trajectory(
                                            run , subrun, event                    // r/s/e
                                            ,i_mcparticle                          // id
                                            ,(particle_start-particle_end).Mag()   // length
                                            ,particle_start                        // start position
                                            ,particle_end                          // end position
                                            );

            Debug(4,"particle %: pdg %: status-code: %, path length: %, process: %"
                  ,i_mcparticle , particle->PdgCode(), particle->StatusCode(), truth_trajectory.GetLength(), particle->Process());
            truth_trajectory.SetTruthProcess( particle->Process() );
            truth_trajectory.SetMCpdgCode( particle->PdgCode() );
            truth_trajectory.SetTruthMomentum( particle -> Momentum() );


            // we would like to keep only the following trajectories
            if (
                // we found that typically,
                // pandoraCosmic and pandoraNu reconstruct trajectories originating from
                // primary particles or Compton scattering or Decay
                (
                 truth_trajectory.GetTruthProcess() == "primary"
//                  ||
//                 truth_trajectory.GetTruthProcess() == "compt"
                 ||
                 truth_trajectory.GetTruthProcess() == "Decay"
                 ||
                 truth_trajectory.GetTruthProcess() == "muMinusCaptureAtRest"
                 ||
                 truth_trajectory.GetTruthProcess() == "muPlusCaptureAtRest"
                 )
                &&
                // kinetic energy greater than 10 MeV, follwing D. Caratelli paper [https://arxiv.org/pdf/1704.02927.pdf]
                (truth_trajectory.GetTruthMomentum().E() - truth_trajectory.GetTruthMomentum().M()) > 0.01
                &&
                // longer than 3 wires-pitch, i.e. 1 cm
                truth_trajectory.GetLength() > 1.0
                &&
                // not a neutron or a neutrino
                truth_trajectory.GetMCpdgCode() != 2112
                &&
                fabs(truth_trajectory.GetMCpdgCode()) != 12
                &&
                fabs(truth_trajectory.GetMCpdgCode()) != 14
                &&
                fabs(truth_trajectory.GetMCpdgCode()) != 16
                ) {
                
                Nparticles_passed_cut++;

                Debug(2,"adding particle % to the truth-traectories list, truth_trajectory.IsTrackContainedSoft(): %",i_mcparticle,truth_trajectory.IsTrackContainedSoft());
                
                
                truth_trajectory.SetTruthMCParticle( i_mcparticle );
                truth_trajectory.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                truth_trajectory.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                truth_trajectory.SetTruthLength();
                truth_trajectory.SetTruthMother( particle -> Mother() );
                
                truth_trajectories.push_back(truth_trajectory);
            }
            
            // loop over all start- and end-points of the truth-particles
            // to find pairs at close proximity in the truth-level
            //            FindTrueCosmicPairs( particles_start , particles_end );
        }
        Debug(1,"Nparticles_passed_cut: %",Nparticles_passed_cut);
        NMCParticles_total += Nparticles_passed_cut;
        
        FindTrajPairs( truth_trajectories );
        // AnalyzeVertices( truth_pairs );
        FilterGoodTrajPair();
        // check if this truth-trajectory pair was reconstructable
        CheckIfReconstructable( truth_pairs );
        // now check if this truth-trajectory pair was reconstructed
        CheckIfReconstructed( truth_pairs );
    }
    StreamTrajectoryPairsToCSV();
    Debug(2,"StreamTrajectoryPairsToCSV()");

    // print out
    if (debug>0) {
        if(!tracks.empty() && !cosmic_tracks.empty()){
            PrintInformation();
        }    else {
            cout << "tracks and cosmic_tracks are empty in this event" << endl;
        }
    }
    fTree -> Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::CheckIfReconstructed(std::vector<pairVertex> & truth_pairs_vector){
    // check if a pair of truth-trajectories at close proximity
    // (1) has been reconstructed
    // (2) if reconsutrcted - also if the two tracks are contained
    
    for (auto & truth_pair:truth_pairs_vector) {
        // to check if the truth-trajectories pair was reconstructed,
        // loop over the two truth-trajectories,
        // and check for each trajectory was reconsutrcted
        bool IsPairReconstructed = true;
        for (auto truth_trajectory: truth_pair.GetTracks()) {
            Debug(2,"lets check if truth_trajectory % is reconstructed",truth_trajectory.GetTrackID());
            bool IsTrajReconstructed = IsTrajecPairReconstructed( truth_trajectory );
            IsPairReconstructed = IsPairReconstructed && IsTrajReconstructed;
        }
        if (IsPairReconstructed) {

            // apart from asking if each of the particles was reconstructed,
            // we want in addition to ask if it was reconsutructed
            // as a pandoraNu pair-vertex
            IsPairReconstructed = IsPairReconstructed && HasMatchingRecoPairVertex(truth_pair);
            truth_pair.SetIsVertexReconstructed( IsPairReconstructed );
            
            if(IsPairReconstructed){
                bool IsPairContained = (truth_pair.GetTracks().at(0).IsTrackContainedSoft()
                                        &&
                                        truth_pair.GetTracks().at(1).IsTrackContainedSoft());
                truth_pair.SetIsVertexContained( IsPairContained );
            }
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::CheckIfReconstructable(std::vector<pairVertex> & truth_pairs_vector){
    // check if a pair of truth-trajectories at close proximity can be reconsutrcted
    // by asking if the "reconstructed" vertex position of the pair is inside the TPC.
    // We use the "reconstructed" vertex position,
    // since it is corrected at the pairVertex::FixTracksDirections() stage
    // while the truth-vertex position is unmodified
    for (auto & truth_pair:truth_pairs_vector) {
        Debug(3,"CheckIfReconstructable(truth_pair %) at position %,%,%"
              ,truth_pair.GetVertexID()
              ,truth_pair.GetPosition().X(),truth_pair.GetPosition().y(),truth_pair.GetPosition().Z());

        truth_pair.CheckIsVertexInTPC();
        truth_pair.SetIsVertexReconstructable( truth_pair.GetIsRecoVertexInTPC() );
    }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::CosmicTracksAnalyzer::HasMatchingRecoPairVertex( pairVertex truth_pair ){
    // March-14, 2018
    // check if the truth-trajectories pair has a matching
    // pandoraNu-tracks pair that was reconstructed at close proximity
    for (auto PandoraNuTrackPair: vertices) {
        Int_t tA_particle = PandoraNuTrackPair.GetTracks().at(0).GetTruthMCParticle();
        Int_t tB_particle = PandoraNuTrackPair.GetTracks().at(1).GetTruthMCParticle();
        Int_t particleA = truth_pair.GetTracks().at(0).GetTruthMCParticle();
        Int_t particleB = truth_pair.GetTracks().at(1).GetTruthMCParticle();
        
        // determine that the pair matches
        // if the particles assigned with the tracks are the same particle in the truth-trajectory pair
        if (
            (tA_particle == particleA && tB_particle == particleB)
            ||
            (tA_particle == particleB && tB_particle == particleA)
            ) {
                return true;
        }
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::CosmicTracksAnalyzer::IsTrajecPairReconstructed(PandoraNuTrack truth_trajectory){
    int fdebug=3;
    // loop over all reconstructed tracks,
    // and check for each track if it matches one of the trajectories
    Debug(fdebug,"ub::CosmicTracksAnalyzer::IsTrajecPairReconstructed(truth_trajectory %)",truth_trajectory.GetTrackID());

    Debug(fdebug,"start: %,%,%",truth_trajectory.GetTruthStartPos().X(),truth_trajectory.GetTruthStartPos().Y(),truth_trajectory.GetTruthStartPos().Z());
    Debug(fdebug,"truth_trajectory.GetMCpdgCode(): %",truth_trajectory.GetMCpdgCode());
    
    for (auto & track : tracks){
        Debug(fdebug,"comparing to see if this can match track %",track.GetTrackID());

        Debug(fdebug,"track.GetMCpdgCode(): %, track.GetTruthMCParticle(): %",track.GetMCpdgCode(),track.GetTruthMCParticle());
        
        if ( track.GetTruthMCParticle() == truth_trajectory.GetTruthMCParticle()
            ) {
            Debug(fdebug,"matched truth_trajectory % with reco. track % ! \n -- - -- - -- - -- -- -- -- -- -- -- --"
                  ,truth_trajectory.GetTrackID(),track.GetTrackID());
            return true;
        }
        Debug(fdebug,"-- - -- - -- - -- -- -- -- -- -- -- -- ");
    }
    return false;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::FindTrajPairs(std::vector<PandoraNuTrack> & traj_vector // tracks to be clustered
){
    Debug(2,"ub::CosmicTracksAnalyzer::FindTrajPairs()");
    // Mar-14, 2018
    // cluster pairs of truth-trajectories at close proximity to pair-vertices
    // since these are truth-trajectories,
    // we want to ask what is the closest distance between their trajectories,
    // and not between their start or end points
    float   closest_SqDist_ij;
    TVector3 QuasiIntersectionPoint;
    geoalgo::Point_t Pt_i, Pt_j;

    // loop over all tracks - to look for partners at close proximity
    for ( size_t i=0; i < traj_vector.size(); i++){
        geoalgo::Line_t line_i(traj_vector[i].GetStartPos().y()
                               ,traj_vector[i].GetStartPos().y()
                               ,traj_vector[i].GetStartPos().z()
                               ,traj_vector[i].GetEndPos().x()
                               ,traj_vector[i].GetEndPos().y()
                               ,traj_vector[i].GetEndPos().z());
        for ( size_t j=i+1 ; j < traj_vector.size() ; j++ ){
            geoalgo::Line_t line_j(traj_vector[j].GetStartPos().x()
                                   ,traj_vector[j].GetStartPos().y()
                                   ,traj_vector[j].GetStartPos().z()
                                   ,traj_vector[j].GetEndPos().x()
                                   ,traj_vector[j].GetEndPos().y()
                                   ,traj_vector[j].GetEndPos().z());
            
            closest_SqDist_ij = mygeoalgo.SqDist(line_i , line_j , Pt_i , Pt_j );
            
            Debug(3,"closest SqDist of i=% & j=%: %",i,j,closest_SqDist_ij);
            if ( closest_SqDist_ij < kMaxInterTrackDistance*kMaxInterTrackDistance ){
                
                TVector3 TVec3_i = Pt_i.ToTLorentzVector().Vect();
                TVector3 TVec3_j = Pt_j.ToTLorentzVector().Vect();
                
                QuasiIntersectionPoint = 0.5*( TVec3_i + TVec3_j );
                
                pairVertex pair = FormPairFrom2Trajs( traj_vector[i] , traj_vector[j]
                                     , QuasiIntersectionPoint );
                pair.SetClosestDistance( (TVec3_i - TVec3_j).Mag() );
                
                Debug(3,"formed pair % from i=% & j=% at %,%,%"
                      ,pair.GetVertexID()
                      ,i,j
                      ,pair.GetPosition().X()
                      ,pair.GetPosition().Y()
                      ,pair.GetPosition().Z());
                
                SearchOtherCloseTrajs( pair , i , j );
                if (pair.GetNtracks() == 2) truth_pairs.push_back( pair );
            }
        }
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
pairVertex ub::CosmicTracksAnalyzer::FormPairFrom2Trajs(PandoraNuTrack ti,PandoraNuTrack tj, TVector3 pos){
    pairVertex pair( run, subrun, event , truth_pairs.size() );
    pair.AddTrack( ti );
    pair.AddTrack( tj );
    pair.SetPosition( pos );
    return pair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::SearchOtherCloseTrajs( pairVertex & pair , int i , int j ){
    Debug(2 , "SearchOtherCloseTrajs( pairVertex % , i=% , j=% )",pair.GetVertexID(),i,j);
    for ( size_t k=0 ; k < truth_trajectories.size() ; k++ ){
        if ((int)k!=i && (int)k!=j){
            geoalgo::Point_t VtxPos( pair.GetPosition() );
            geoalgo::Line_t line_k(truth_trajectories[k].GetStartPos().x()
                                   ,truth_trajectories[k].GetStartPos().y()
                                   ,truth_trajectories[k].GetStartPos().z()
                                   ,truth_trajectories[k].GetEndPos().x()
                                   ,truth_trajectories[k].GetEndPos().y()
                                   ,truth_trajectories[k].GetEndPos().z());

            if(debug>3){
                geoalgo::Point_t Pt_k = mygeoalgo.ClosestPt( line_k , VtxPos );
                TVector3 Pt_k_vect = Pt_k.ToTLorentzVector().Vect();
                TVector3 Pcm_tracks = (1./(Ntracks-1))*pair.GetPosition() + (1./Ntracks)*Pt_k_vect;

                Debug(0,"distance of line k=% from the vertex at %,%,% is % at %,%,%"
                  ,k,pair.GetPosition().X(),pair.GetPosition().Y(),pair.GetPosition().Z()
                  ,mygeoalgo.SqDist( line_k, VtxPos)
                      ,Pt_k_vect.x(),Pt_k_vect.y(),Pt_k_vect.z());
            }
            if( mygeoalgo.SqDist( line_k, VtxPos) < kMaxInterTrackDistance*kMaxInterTrackDistance ){
                geoalgo::Point_t Pt_k = mygeoalgo.ClosestPt( line_k , VtxPos );
                pair.AddTrack( truth_trajectories[k] );
                int Ntracks = (int)pair.GetTracks().size();
                TVector3 Pt_k_vect = Pt_k.ToTLorentzVector().Vect();
                TVector3 Pcm_tracks = (1./(Ntracks-1))*pair.GetPosition() + (1./Ntracks)*Pt_k_vect;
                pair.SetPosition( Pcm_tracks );
                Debug(3,"Adding traj k=% to traj j=% and i=%, that now has % trajectories and is located at %,%,%..."
                      ,k,j,i,pair.GetTracks().size()
                      ,pair.GetPosition().X(),pair.GetPosition().Y(),pair.GetPosition().Z());
            }
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::FilterGoodTrajPair(){
    
    std::vector<pairVertex> tmp_pairs = truth_pairs;
    truth_pairs.clear();
    
    Debug(2 , "ub::CosmicTracksAnalyzer::FilterGoodTrajPair()");
    for (auto & pair:tmp_pairs) {
        
        if (debug>2) {
            pair.Print(false, false);
            cout << "filtering Trajectories-pair "
            << pair.GetVertexID()
            << " with " << pair.GetNtracks() << " trajectories: ";
            
            for (auto t: pair.GetTracks()) {
                cout
                << t.GetTrackID()
                <<"(pdg "<< t.GetMCpdgCode()
                <<", MCParticle " << t.GetTruthMCParticle() <<"),  ";
            }
            cout << endl;
        }

        
        // vertices with only two tracks at close proximity and nothing else
        if (pair.GetNtracks() == 2){
            // we have to 'assign' a muon and a proton track,
            // so lets do it arbitrarily
            pair.AssignMuonTrack (pair.GetTracks().at(0));
            pair.AssignProtonTrack (pair.GetTracks().at(1));

            pair.FixTracksDirections ();
            Debug(2,"Kept Trajectories-pair % at position %,%,% as a good pair!"
                  ,pair.GetVertexID(),pair.GetPosition().X(),pair.GetPosition().y(),pair.GetPosition().Z());
            
            truth_pairs.push_back( pair );
        } else {
            Debug(2,"Filtered out Trajectories-pair % at position %,%,%, being not a good pair..."
                  ,pair.GetVertexID(),pair.GetPosition().X(),pair.GetPosition().y(),pair.GetPosition().Z());
        }
    }
    Debug(2 , "-- -- - -- - -- -- - --- - -- -- -- -- - -- -- - -- -- - - ----");
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::FindTwoTracksVertices(std::vector<PandoraNuTrack> & tracks_vector // tracks to be clustered
                                                     ,std::vector<pairVertex> & vertices_vector  // vertices vector to hold the pairs
){
    Debug(2,"ub::CosmicTracksAnalyzer::FindTwoTracksVertices()");
    
    // Jan-14, 2018
    // cluster pairs of tracks at close proximity to two-tracks vertices
    bool    FoundCloseTracks=false;
    float   closest_distance_ij;
    TVector3 vertex_position;
    
    // loop over all tracks - to look for partners at close proximity
    for ( size_t i=0; i < tracks_vector.size(); i++){
        if (!tracks_vector[i].IsTrackContainedSoft()) continue; // if this track is not contained, don't even bother...
        pairVertex vertex( run, subrun, event , vertices_vector.size() );
        vertex.AddTrack( tracks_vector[i] );
        FoundCloseTracks = false;
        for ( size_t j=i+1 ; j < tracks_vector.size() ; j++ ){
            if (tracks_vector[j].IsTrackContainedSoft() && j!=i){ // if this track is not contained, don't even bother...
                if (!FoundCloseTracks){
                    std::string StartOrEnd = "None";
                    closest_distance_ij = tracks_vector[i].ClosestDistanceToOtherTrack( tracks_vector[j] , &StartOrEnd);
                    if ( closest_distance_ij < kMaxInterTrackDistance ){
                        // if trk-i and trk-j are close enough - we want to keep this vertex
                        vertex.AddTrack( tracks_vector[j] );
                        FoundCloseTracks = true;
                        vertex.SetPosition( (StartOrEnd.compare("Start")==0)
                                           ?
                                           tracks_vector[i].GetStartPos()
                                           :
                                           tracks_vector[i].GetEndPos());
                    } // end if keep vertex
                } // end if !FoundCloseTracks
                else { // (namely if FoundCloseTracks==true)
                    if ( tracks_vector[j].DistanceFromPoint(vertex.GetPosition()) < kMaxInterTrackDistance ){
                        vertex.AddTrack( tracks_vector[j] );
                    }
                }
            }
        } // end loop on track-j
        if (FoundCloseTracks) { // plug into vertices list
            vertex.BuildVertexIDFromTracks();
            vertices_vector.push_back( vertex );
        }
    } // end loop on track-i
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::AnalyzeVertices( std::vector<pairVertex> & vertices_vector  ){
    if (vertices_vector.size()>0){
        for (auto & v:vertices_vector){
            // after fixing the vertext position, remove far tracks
            Debug( 4 , "v.RemoveFarTracks( kMaxInterTrackDistance );");
            if (v.GetTracks().size()<2) continue;
            v.RemoveFarTracks( kMaxInterTrackDistance );
            
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
            
            Debug(4,"analyzed vertex % at position %,%,%"
                  ,v.GetVertexID(),v.GetPosition().X(),v.GetPosition().y(),v.GetPosition().Z());
        }
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::FilterGoodPairVertices(std::vector<pairVertex> & vertices_vector ){
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    std::vector<pairVertex> tmp_vertices = vertices_vector;
    vertices_vector.clear();
    trkf::TrackMomentumCalculator trkm;

    
    Debug(2 , "ub::CosmicTracksAnalyzer::FilterGoodPairVertices()");

    for (auto & v:tmp_vertices) {
        if (debug>2) {
            cout << "filtering vertex "<< v.GetVertexID() << " with tracks: ";
            for (auto t: v.GetTracks()) {
                cout
                << t.GetTrackID()
                <<"(pdg "<< t.GetMCpdgCode()
                <<", MCParticle " << t.GetTruthMCParticle() <<"),  ";
            }
            cout << endl;
        }
        // vertices with only two tracks at close proximity and nothing else
        if (
            // we are looking for clusters of only two tracks that are fully contained
            v.GetNtracks() == 2
            // if there is a semi-contained track, with start/end point too close to the vertex, we don't want the vertex...
            &&
            v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size() == 0
            ){
            v.AssignMuonTrack (v.GetTracks().at(0));
            v.AssignProtonTrack (v.GetTracks().at(1));

            v.FixTracksDirections ();
            
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
            Debug(2,"Filtered vertex % at position %,%,% as a good pair-vertex!"
                  ,v.GetVertexID(),v.GetPosition().X(),v.GetPosition().y(),v.GetPosition().Z());
            vertices_vector.push_back( v );
        }
    }
    Debug(2 , "-- -- - -- - -- -- - --- - -- -- -- -- - -- -- - -- -- - - ----");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::MatchPairVertices(){
    // check if the pandoraCosmic vertices match the pandoraNu vertices
    int fdebug=3;
    Debug(fdebug,"ub::CosmicTracksAnalyzer::MatchPairVertices()");
    for( auto pandoraCosmic_vertex: cosmic_vertices){
        auto pandoraCosmic_vertex_tracks = pandoraCosmic_vertex.GetTracks();
        if (debug>fdebug) {
            cout << "pandoraCosmic_vertex" << endl;
            pandoraCosmic_vertex.Print(false,false);
        }
        
        
        int N_pandoraCosmic_tracks_in_pandoraNu_vertex = 0;
        for( auto pandoraNu_vertex: vertices){
            if (debug>fdebug) {
                cout << "pandoraNu_vertex" << endl;
                pandoraNu_vertex.Print(false,false);
            }
            
            // the vertices are matched if all of the tracks in the pandoraCosmic-vertex exist in the pandoraNu-vertex
            for (auto pandoraCosmic_track:pandoraCosmic_vertex_tracks) {
                Debug(fdebug,Form("looking for match for track %d in pandoraNu_vertex",pandoraCosmic_track.GetTrackID()));
                if ( pandoraNu_vertex.IncludesTrack( pandoraCosmic_track.GetTrackID() ) ) {
                    Debug(fdebug,Form("!found a match for track %d in pandoraNu_vertex!",pandoraCosmic_track.GetTrackID()));
                    N_pandoraCosmic_tracks_in_pandoraNu_vertex += 1;
                } else {
                    Debug(fdebug,Form("!did not find a match for track %d in pandoraNu_vertex!",pandoraCosmic_track.GetTrackID()));
                }
            }
            
        }
        if (N_pandoraCosmic_tracks_in_pandoraNu_vertex == (int)pandoraCosmic_vertex_tracks.size()){
            Debug(fdebug,Form("Reproduction: cosmic vertex %d was reproduced also in pandoraNu stage",pandoraCosmic_vertex.GetVertexID()));
        } else {
            Debug(fdebug,Form("Cosmic rejection: cosmic vertex %d was not reproduced also in pandoraNu stage",pandoraCosmic_vertex.GetVertexID()));
        }
    }
    Debug(fdebug,"ub::CosmicTracksAnalyzer::MatchPairVertices()");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::HeaderVerticesInCSV(){
    
    // we will stream to the output file all the vertices,
    // both those of PandoraCosmic and those from PandoraNu,
    // and all the relevant information
    // ( for example, the reconstructed distacne between the tracks in the pair )
    // Then, in the analysis stage
    // we can separate to the desired vertices and draw the distributions
    vertices_counter = 0;
    
    vertices_file
    << "run" << "," << "subrun" << "," << "event" << ","
    << "vertex_id" << ","
    << "Ntracks" << ","
    << "pdg_1" << ","
    << "pdg_2" << ","
    << "MCParticle_1" << ","
    << "MCParticle_2" << ","
    ;
    
    
    // which tracks are populating this vertex topology
    vertices_file
    << "isPandoraCosmic" << ","
    << "isPandoraNu" << ","
    << "isBrokenTrajectory" << ",";
    
    // start/end points, for FV cuts
    vertices_file
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
    << "endz_assigned_proton" << ",";


    // reconstructed distance between the tracks
    vertices_file
    << "distance" ;
    
    
    // finish
    vertices_file << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::StreamVerticesToCSV(){
    // Oct-22, 2017
    // whatever you add here - must add also in header - ub::CosmicTracksAnalyzer::HeaderVerticesInCSV()
    for (int i_vertices_vector = 0; i_vertices_vector < 2 ; i_vertices_vector++){
        
        std::vector<pairVertex> vertices_vector = (i_vertices_vector==0) ? cosmic_vertices : vertices;
        bool isPandoraCosmic = (i_vertices_vector==0) ? true : false;
        
        for (auto v : vertices_vector) {
            if (!isPandoraCosmic) {vertices_counter++;}
            
            vertices_file
            << v.GetRun() << "," << v.GetSubrun() << "," << v.GetEvent() << ","
            << v.GetVertexID() << ","
            << v.GetTracks().size() << ","
            << v.GetTracks().at(0).GetMCpdgCode() << ","
            << v.GetTracks().at(1).GetMCpdgCode() << ","
            << v.GetTracks().at(0).GetTruthMCParticle() << ","
            << v.GetTracks().at(1).GetTruthMCParticle() << ","
            ;
            
            // which tracks are populating this vertex topology
            vertices_file
            << isPandoraCosmic << ","
            << (!isPandoraCosmic) << ","
            << v.GetIsBrokenTrajectory() << ",";

            
            // start/end points, for FV cuts
            vertices_file
            << v.GetTrack_muCandidate().GetStartPos().x() << ","
            << v.GetTrack_muCandidate().GetStartPos().y() << ","
            << v.GetTrack_muCandidate().GetStartPos().z() << ","
            << v.GetTrack_pCandidate().GetStartPos().x() << ","
            << v.GetTrack_pCandidate().GetStartPos().y() << ","
            << v.GetTrack_pCandidate().GetStartPos().z() << ","
            << v.GetTrack_muCandidate().GetEndPos().x() << ","
            << v.GetTrack_muCandidate().GetEndPos().y() << ","
            << v.GetTrack_muCandidate().GetEndPos().z() << ","
            << v.GetTrack_pCandidate().GetEndPos().x() << ","
            << v.GetTrack_pCandidate().GetEndPos().y() << ","
            << v.GetTrack_pCandidate().GetEndPos().z() << ",";

            
            // reconstructed distance between the tracks
            for ( auto d : v.Get_distances_ij()){
                vertices_file << d ; //  for more than 2 tracks add << ";" ;
            }
            
            
            // finish
            vertices_file << endl;
        }
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::HeaderTrajectoryPairsInCSV(){
    // Mar-13, 2018
    truth_pairs_ctr = 0;
    truth_pairs_ctr_reconstructed = 0;
    
    truth_pairs_file
    << "run" << "," << "subrun" << "," << "event" << ","
    << "vertex_id" << ","
    << "pdg_1" << ","
    << "pdg_2" << ","
    << "MCParticle_1" << ","
    << "MCParticle_2" << ","
    << "MCParticle_1_length" << ","
    << "MCParticle_2_length" << ","
    << "MCParticle_1_mother" << ","
    << "MCParticle_2_mother" << ","
    << "MCParticle_1_process" << ","
    << "MCParticle_2_process" << ",";
    
    // position
    truth_pairs_file
    << "position_x" << "," << "position_y" << "," << "position_z" << ","
    << "MCParticle_1_startx" << "," << "MCParticle_1_starty" << "," << "MCParticle_1_startz" << ","
    << "MCParticle_1_endx" << "," << "MCParticle_1_endy" << "," << "MCParticle_1_endz" << ","
    << "MCParticle_2_startx" << "," << "MCParticle_2_starty" << "," << "MCParticle_2_startz" << ","
    << "MCParticle_2_endx" << "," << "MCParticle_2_endy" << "," << "MCParticle_2_endz" << ",";
    
    // reconstruction
    truth_pairs_file
    << "IsVertexReconstructable" << ","
    << "IsVertexReconstructed" << ","
    << "IsVertexContained" << ","
    << "distance" ;
    
    // finish
    truth_pairs_file << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::StreamTrajectoryPairsToCSV(){
    // Mar-13, 2018
    for (auto pair : truth_pairs) {
        truth_pairs_ctr ++;
        if (pair.GetIsVertexReconstructable())  truth_pairs_ctr_reconstructable++;
        if (pair.GetIsVertexReconstructed())    truth_pairs_ctr_reconstructed++;
        if (pair.GetIsVertexContained())        truth_pairs_ctr_contained++;
        
        truth_pairs_file
        << pair.GetRun() << "," << pair.GetSubrun() << "," << pair.GetEvent() << ","
        << pair.GetVertexID() << ","
        << pair.GetTracks().at(0).GetMCpdgCode() << ","
        << pair.GetTracks().at(1).GetMCpdgCode() << ","
        << pair.GetTracks().at(0).GetTruthMCParticle() << ","
        << pair.GetTracks().at(1).GetTruthMCParticle() << ","
        << pair.GetTracks().at(0).GetTruthLength() << ","
        << pair.GetTracks().at(1).GetTruthLength() << ","
        << pair.GetTracks().at(0).GetTruthMother() << ","
        << pair.GetTracks().at(1).GetTruthMother() << ","
        << pair.GetTracks().at(0).GetTruthProcess() << ","
        << pair.GetTracks().at(1).GetTruthProcess() << ","
        ;
        
        // position
        truth_pairs_file
        << pair.GetPosition().X() << ","
        << pair.GetPosition().Y() << ","
        << pair.GetPosition().Z() << ","
        << pair.GetTracks().at(0).GetStartPos().X() << ","
        << pair.GetTracks().at(0).GetStartPos().Y() << ","
        << pair.GetTracks().at(0).GetStartPos().Z() << ","
        << pair.GetTracks().at(0).GetEndPos().X() << ","
        << pair.GetTracks().at(0).GetEndPos().Y() << ","
        << pair.GetTracks().at(0).GetEndPos().Z() << ","
        << pair.GetTracks().at(1).GetStartPos().X() << ","
        << pair.GetTracks().at(1).GetStartPos().Y() << ","
        << pair.GetTracks().at(1).GetStartPos().Z() << ","
        << pair.GetTracks().at(1).GetEndPos().X() << ","
        << pair.GetTracks().at(1).GetEndPos().Y() << ","
        << pair.GetTracks().at(1).GetEndPos().Z() << ","
        ;

        // reconstruction
        truth_pairs_file
        << pair.GetIsVertexReconstructable() << ","
        << pair.GetIsVertexReconstructed() << ","
        << pair.GetIsVertexContained() << ","
        << pair.GetClosestDistance() ;
        
        // finish
        truth_pairs_file << endl;
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::PrintInformation(){
    
    PrintXLine();
    Printf( "processed so far %d events", (int)(fTree->GetEntries()) );
    SHOW3( run , subrun , event );
    
    
    if(!truth_trajectories.empty()){
        PrintHeader(Form("%d truth-trajectories",(int)truth_trajectories.size()));
        for (auto t: truth_trajectories) {
            if (DoPrintTruthTrajectories) t.Print( true );
            else cout << t.GetTrackID() << "\t";
        }
        cout << endl;
    } else { PrintHeader("no truth-trajectories");    }
    
    if(!truth_pairs.empty()){
        PrintHeader(Form("%d truth-trajectorirs pair-vertices",(int)truth_pairs.size()));
        for (auto pair: truth_pairs) {
            pair.Print(false, false);
            SHOW2(pair.GetIsVertexReconstructable(),pair.GetIsVertexReconstructed());
        }
    } else { PrintHeader("no truth-trajectorirs pair-vertices");}

    if(!tracks.empty()){
        PrintHeader(Form("%d pandoraNu-tracks",(int)tracks.size()));
        for (auto t: tracks) {
            if (DoPrintPandoraNuTracks) t.Print( true );
            else cout << t.GetTrackID() << "\t";
        }
        cout << endl;
    } else { PrintHeader("no pandoraNu-tracks");    }

    if(!vertices.empty()){
        PrintHeader(Form("%d pandoraNu pair-vertices",(int)vertices.size()));
        for (auto v: vertices) {
            v.Print(true, true);
            SHOW(v.GetIsBrokenTrajectory());
        }
    } else { PrintHeader("no pandoraNu pair-vertices"); }

    
    //    if(!cosmic_tracks.empty()){
    //        PrintHeader(Form("%d pandoraCosmic-tracks",(int)cosmic_tracks.size()));
    //        for (auto t: cosmic_tracks) {
    //            if (DoPrintCosmicTracks) t.Print( true );
    //            else cout << t.GetTrackID() << "\t";
    //        }
    //        cout << endl;
    //    } else { PrintHeader("no pandoraCosmic-tracks");    }
    //    if(!cosmic_vertices.empty()){
    //        PrintHeader(Form("%d pandoraCosmic pair-vertices",(int)cosmic_vertices.size()));
    //        for (auto v: cosmic_vertices) {
    //            v.Print(false, false);
    //        }
    //    } else { PrintHeader("no pandoraCosmic pair-vertices");}
    

    
    
    // time stamp
    PrintLine();
    end_ana_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
    std::cout << "\033[33m"
    << "finished analysis of this event at " << std::ctime(&end_time)
    << "time elapsed: " << elapsed_seconds.count() << "s"
    << "\033[31m" << endl;

    cout << "wrote " << vertices_counter << " vertices to output file " << endl;
    EndEventBlock();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::beginJob(){
    
    // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
    for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
        NwiresBox[i_box_size] = MinNwiresBox + i_box_size * dNwiresBox;
        NticksBox[i_box_size] = MinNticksBox + i_box_size * dNticksBox;
    }

    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of genie interactions");
    fTree->Branch("run"             ,&run           ,"run/I");
    fTree->Branch("subrun"          ,&subrun        ,"subrun/I");
    fTree->Branch("event"           ,&event         ,"event/I");
    fTree->Branch("Ntracks"         ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("NCosmicTracks"   ,&NCosmicTracks ,"NCosmicTracks/I");
    fTree->Branch("Nhits"           ,&Nhits_stored  ,"Nhits/I");
    
    // my objects
    fTree->Branch("cosmic_tracks"       ,&cosmic_tracks);
    fTree->Branch("cosmic_vertices"     ,&cosmic_vertices);
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("vertices"            ,&vertices);
    fTree->Branch("hits"                ,&hits);
    fTree->Branch("truth_trajectories"  ,&truth_trajectories);

    
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("pot",&pot,"pot/D");
    fPOTTree->Branch("run",&run,"run/I");
    fPOTTree->Branch("subrun",&subrun,"subrun/I");

    // output csv file
    //    vertices_file.open("/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+fDataSampleLabel+"_vertices.csv");
    vertices_file.open(fDataSampleLabel+"_pandora_vertices.csv");
    cout << "opened vertices file: "+fDataSampleLabel+"_pandora_vertices.csv" << endl;
    HeaderVerticesInCSV();
    
    truth_pairs_file.open(fDataSampleLabel+"_truth_trajectories_pairs.csv");
    cout << "opened vertices file: "+fDataSampleLabel+"_truth_trajectories_pairs.csv" << endl;
    HeaderTrajectoryPairsInCSV();
    
    
    pot_total = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::endJob(){
    Debug(3,"ub::CosmicTracksAnalyzer::endJob()");
    
    std::time_t now_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    summary_file.open(fDataSampleLabel+"_summary.csv");
    summary_file << "time" << ","
    << "Nevents" << ","
    << "NPandoraNuTracks" << ","
    << "NPandoraCosmicTracks" << ","
    << "Nvertices" << ","
    << "Nhits" << ","
    << "NMCParticles" << ","
    << "Ntruth_pairs" << ","
    << "Ntruth_pairs_reconstructable" << ","
    << "Ntruth_pairs_reconstructed" << ","
    << "Ntruth_pairs_contained"
    << endl;
    
    std::string sTimeS = std::ctime(&now_time);
    
    summary_file << sTimeS.substr(0,sTimeS.length()-1) << ","
    << fTree->GetEntries() << ","
    << Ntracks_total << ","
    << NCosmicTracks_total << ","
    << vertices_counter << ","
    << Nhits_tot << ","
    << NMCParticles_total << ","
    << truth_pairs_ctr << ","
    << truth_pairs_ctr_reconstructable << ","
    << truth_pairs_ctr_reconstructed << ","
    << truth_pairs_ctr_contained
    << endl;
    
    summary_file.close();
    
    cout << "Ended job. Wrote summary into file: "+fDataSampleLabel+"_summary.csv" << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::endSubRun(const art::SubRun& sr){
    
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
    SHOW2( NCosmicTracks_total, Ntracks_total );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel        = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel    = p.get< std::string >("GenieGenModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
    fDataSampleLabel        = p.get< std::string >("DataSampleLabel");
    fPOTModuleLabel         = p.get< std::string >("POTModuleLabel");
    fFlashModuleLabel       = p.get< std::string >("FlashModuleLabel");
    debug                   = p.get< int >("VerbosityLevel");
    fHitParticleAssnsModuleLabel = p.get< std::string >("HitParticleAssnsModuleLabel");
    fG4ModuleLabel          = p.get< std::string >("G4ModuleLabel","largeant");
    fMCTrackModuleLabel     = p.get< std::string >("MCTrackModuleLabel","mcreco");
    fCosmicTrackModuleLabel = p.get< std::string >("CosmicTrackModuleLabel");
    DoPrintTruthTrajectories= p.get< bool >("DoPrintTruthTrajectories");
    DoPrintCosmicTracks     = p.get< bool >("DoPrintCosmicTracks");
    DoPrintPandoraNuTracks  = p.get< bool >("DoPrintPandoraNuTracks");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::CosmicTracksAnalyzer::ResetVars(){
    
    MCmode = false;
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = 0 ;
    pot = 0;
    tracks.clear();
    cosmic_tracks.clear();
    vertices.clear();
    cosmic_vertices.clear();
    hits.clear();
    flashes.clear();
    truth_trajectories.clear();
    truth_pairs.clear();
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::CosmicTracksAnalyzer)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -