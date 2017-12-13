////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEGENIE
// Plugin Type: analyzer (art v2_05_00)
// File:        ErezCCQEGENIE_module.cc
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
namespace ub { class ErezCCQEGENIE; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::ErezCCQEGENIE : public art::EDAnalyzer {
public:
    explicit ErezCCQEGENIE(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEGENIE(ErezCCQEGENIE const &) = delete;
    ErezCCQEGENIE(ErezCCQEGENIE &&) = delete;
    ErezCCQEGENIE & operator = (ErezCCQEGENIE const &) = delete;
    ErezCCQEGENIE & operator = (ErezCCQEGENIE &&) = delete;
    
    // Required functions.
    void                   analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                  beginJob () override;
    void               reconfigure (fhicl::ParameterSet const& p) override;
    
    void                 endSubRun (const art::SubRun& sr);

    
    // functionallity
    void       FindGENIECCVertices ();
    void          PrintInformation ();
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
    int     N_CC_interactions;
    int     CC_interactions_ctr;
    
    int     NwiresBox[N_box_sizes], NticksBox[N_box_sizes];
    
    double  pot, pot_total;

    // my objects
    std::vector<PandoraNuTrack>     tracks;
    std::vector<hit>                hits;
    std::vector<GENIEinteraction>   genie_interactions, genie_CC_interactions;
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
ub::ErezCCQEGENIE::ErezCCQEGENIE(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::analyze(art::Event const & evt){
    
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
    if (debug>5) SHOW(Nflashes);
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
    // tracks information
    // ----------------------------------------
    Ntracks = tracklist.size();
    if (debug>4){
        SHOW4( run , subrun , event , Ntracks );
        for (int i=0 ; i<Ntracks ; i++) cout << tracklist[i]->ID() << "\t";
        cout << endl;
    }
    for(int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i ){
        recob::Track::Point_t start_pos, end_pos;
        std::tie( start_pos, end_pos ) = tracklist[i]->Extent();
        
        Debug(5 , Form("analyzing track %d in this event",tracklist[i]->ID()) );
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
        if (!isdata&&fmth.isValid()){
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
            if (particle){
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
        
        Debug(5 , Form("adding track %d to list of tracks in this event",track.GetTrackID()) );
        tracks.push_back( track );
        if (debug>5){
            cout << "track list includes:" << endl;
            for (auto t: tracks) cout << t.GetTrackID() << "\t";
            cout << endl;
        }
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
            
            for( int mcevent_id = 0; (mcevent_id < mcevts_truth) && (mcevent_id < kMaxTruth) ; mcevent_id++ ){
                art::Ptr<simb::MCTruth> mctruth = mclist[mcevent_id];
                
                if (mctruth->Origin() == simb::kBeamNeutrino){
                    
                    GENIEinteraction genie_interaction( run , subrun , event , mcevent_id );
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
                                    // Printf("found a primary-track match! plugging mcevent_id=%d into track %d",mcevent_id,track.GetTrackID());
                                    genie_interaction.AddTrack ( track );
                                    track.SetMCeventID( mcevent_id );
                                }
                            }

                        } // for particle
                    }
                    genie_interaction.SortNucleons();
                    genie_interaction.ComputePmissPrec();
                    genie_interaction.SetTruthTopology();
                    genie_interaction.SetReconstructedTopology();

                    if (debug>6){
                        cout << "finished analyzing genie interaction " << mcevent_id << ", genie track list includes:" << endl;
                        for (auto t: genie_interaction.GetTracks()) cout << t.GetTrackID() << "\t";
                        cout << endl;
                    }

                    genie_interactions.push_back( genie_interaction );
                    
                }//mctruth->Origin()
            }
        }
    }//is neutrino

    if (debug>4){
        cout << "after MC truth information loop ended, track list includes:" << endl;
        for (auto t: tracks) cout << t.GetTrackID() << "\t";
        cout << endl;
    }
    FindGENIECCVertices();
    
    if (debug>4){
        cout << "after FindGENIECCVertices() ended, track list includes:" << endl;
        for (auto t: tracks) cout << t.GetTrackID() << "\t";
        cout << endl;
    }
    if(!genie_CC_interactions.empty() && !tracks.empty()){
        PrintInformation();
    } else {
        cout << "either genie_CC_interactions or tracks is empty in this event" << endl;
    }
    fTree -> Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::FindGENIECCVertices(){
    
    // out of all GENIE interactions in this event (there is typically 1 per event)
    // keep only ones that are CC interactions (to this end we use genie's "ccnc" flag: CC=0/NC=1
    // store them into genie_CC_interactions
    for (auto g:genie_interactions){
        if ( g.GetCCNC() == 0 ) {
            genie_CC_interactions.push_back(g);
        }
    }
    if (debug>4){
        cout << "after generating genie_CC_interactions, track list includes:" << endl;
        for (auto t: tracks) cout << t.GetTrackID() << "\t";
        cout << endl;
    }
    N_CC_interactions = (int)genie_CC_interactions.size();
    
    // output to csv file
    StreamVerticesToCSV();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::HeaderVerticesInCSV(){
    
    CC_interactions_ctr = 0;
    
    vertices_file
    << "run" << "," << "subrun" << "," << "event" << "," ;
    
    // truth topology
    vertices_file
    << "IsCC_Np_200MeVc" << ","
    << "IsCC_1p_200MeVc" << ","
    << "IsCC_1p_200MeVc_0pi" << ","
    << "IsCCQE" << ",";
    
    // reconstructed topology
    vertices_file
    << "IsVertexContained" << ","
    << "IsInActiveVolume" << ","
    << "Is1mu1p" << ","
    << "Is_mu_TrackReconstructed" << ","
    << "Is_mu_TrackInFV" << ","
    << "Is_p_TrackReconstructed" << ","
    << "Is_p_TrackInFV" << ","
    << "IsVertexReconstructed" << ","
    << "IsVertexInFV" << ",";
    
    // general for CC events
    vertices_file
    << "truth_Pmu" << ","
    << "truth_Pmu_theta" << ","
    << "truth_Pp" << ","
    << "truth_Pp_theta" << ",";
    
    // relevant truth-information
    vertices_file
    << "truth_Ev" << ","
    << "truth_Q2" << ",";
    
    
    // only for 1mu-1p vertices
    vertices_file
    << "reconstructed mu-p distance" ;
    
    
    // finish
    vertices_file << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::StreamVerticesToCSV(){
    // Oct-22, 2017
    // whatever you add here - must add also in header - ub::ErezCCQEGENIE::HeaderVerticesInCSV()
    for (auto g:genie_CC_interactions){
        
        CC_interactions_ctr++;
        
        vertices_file
        << g.GetRun() << "," << g.GetSubrun() << "," << g.GetEvent() << "," ;
        
        
        // truth topology
        vertices_file
        << g.GetIsCC_Np_200MeVc() << ","
        << g.GetIsCC_1p_200MeVc() << ","
        << g.GetIsCC_1p_200MeVc_0pi() << ","
        << g.GetIsCCQE() << ",";
        
        // reconstructed topology
        vertices_file
        << g.GetVertexContained() << ","
        << g.GetIsInActiveVolume() << ","
        << g.GetIs1mu1p() << ","
        << g.GetIs_mu_TrackReconstructed() << ","
        << g.GetIs_mu_TrackInFV() << ","
        << g.GetIs_p_TrackReconstructed() << ","
        << g.GetIs_p_TrackInFV() << ","
        << g.GetIsVertexReconstructed() << ","
        << g.GetIsVertexInFV() << ",";
        
        
        // general for CC events
        vertices_file
        << g.GetPmu().P() << ","
        << g.GetPmu().Theta() << ","
        << g.GetPp().P() << ","
        << g.GetPp().Theta() << ",";
        
        // relevant truth-information
        vertices_file
        << g.GetEv() << ","
        << g.GetQ2() << ",";
        
        // only for 1mu-1p vertices
        vertices_file << g.GetReco_mu_p_distance();
        
        
        
        // finish
        vertices_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::PrintInformation(){
    
    PrintXLine();
    Printf( "processed so far %d events", (int)(fTree->GetEntries()) );
    SHOW3( run , subrun , event );
    

    if(!genie_CC_interactions.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << genie_CC_interactions.size() << " genie CC-interactions\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto g: genie_CC_interactions) {
            g.Print(true);
        }
    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no interactions\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
    
    
    //    if(!flashes.empty()){
    //        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << flashes.size() << " flashes\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
    //        for (auto f: flashes) {
    //            f.Print();
    //        }
    ////    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no flashes\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
    //    if(!tracks.empty()){
    //        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << tracks.size() << " pandoraNu tracks\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
    //        for (auto t: tracks) cout << t.GetTrackID() << "\t";
    //        cout << endl;
    //        for (auto t: tracks) {
    //            t.Print( true );
    //        }
    //    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no reco tracks\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}

    
    // time stamp
    PrintLine();
    end_ana_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
    std::cout << "\033[33m"
    << "finished analysis of this event at " << std::ctime(&end_time)
    << "time elapsed: " << elapsed_seconds.count() << "s"
    << "\033[31m" << endl;

    cout << "wrote " << CC_interactions_ctr << " vertices to output file " << endl;
    EndEventBlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::beginJob(){
    
    // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
    for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
        NwiresBox[i_box_size] = MinNwiresBox + i_box_size * dNwiresBox;
        NticksBox[i_box_size] = MinNticksBox + i_box_size * dNticksBox;
    }

    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of genie interactions");
    fTree->Branch("run"     ,&run           ,"run/I");
    fTree->Branch("subrun"  ,&subrun        ,"subrun/I");
    fTree->Branch("event"   ,&event         ,"event/I");
    fTree->Branch("Ntracks" ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("Nhits"   ,&Nhits_stored  ,"Nhits/I");
    fTree->Branch("N_CC_interactions",&N_CC_interactions    ,"N_CC_interactions/I");
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("hits"                ,&hits);
    fTree->Branch("genie_interactions"  ,&genie_interactions);
    fTree->Branch("genie_CC_interactions"  ,&genie_CC_interactions);

    
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("pot",&pot,"pot/D");
    fPOTTree->Branch("run",&run,"run/I");
    fPOTTree->Branch("subrun",&subrun,"subrun/I");

    // output csv file
    //    vertices_file.open("/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+fDataSampleLabel+"_vertices.csv");
    vertices_file.open(fDataSampleLabel+"_genie.csv");
    cout << "opened vertices file: "+fDataSampleLabel+"_genie.csv" << endl;
    HeaderVerticesInCSV();
    
    pot_total = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEGENIE::endSubRun(const art::SubRun& sr){
    
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
void ub::ErezCCQEGENIE::reconfigure(fhicl::ParameterSet const & p){
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
void ub::ErezCCQEGENIE::ResetVars(){
    
    MCmode = false;
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = N_CC_interactions = 0 ;
    pot = 0;
    tracks.clear();
    hits.clear();
    flashes.clear();
    genie_interactions.clear();
    genie_CC_interactions.clear();
    
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezCCQEGENIE)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -