////////////////////////////////////////////////////////////////////////
// Class:       ErezTracksAnalyzerNewMCTruthMatching
// Plugin Type: analyzer (art v2_05_00)
// File:        ErezTracksAnalyzerNewMCTruthMatching_module.cc
//
// Generated at Fri Oct 13 13:12:16 2017 by Erez Cohen using cetskelgen
// from cetlib version v1_21_00.
//
// The purpose of this module is only to collect tracks
// with reconstructed features
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
//#include "larsim/MCCheater/BackTracker.h"
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

// for the new MC truth matching (by Wes)
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"


// constants
constexpr int debug          = 1;
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTruth      = 100;
constexpr int kMaxNgenie     = 100;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class ErezTracksAnalyzerNewMCTruthMatching; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::ErezTracksAnalyzerNewMCTruthMatching : public art::EDAnalyzer {
public:
    explicit ErezTracksAnalyzerNewMCTruthMatching(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezTracksAnalyzerNewMCTruthMatching(ErezTracksAnalyzerNewMCTruthMatching const &) = delete;
    ErezTracksAnalyzerNewMCTruthMatching(ErezTracksAnalyzerNewMCTruthMatching &&) = delete;
    ErezTracksAnalyzerNewMCTruthMatching & operator = (ErezTracksAnalyzerNewMCTruthMatching const &) = delete;
    ErezTracksAnalyzerNewMCTruthMatching & operator = (ErezTracksAnalyzerNewMCTruthMatching &&) = delete;
    
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
    void         HeaderTracksInCSV ();
    void         StreamTracksToCSV ();
    
    
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
    
    // For keeping track of the replacement backtracker
    std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;

    Int_t run, subrun, event;
    
    short   isdata;
    
    bool    MCmode;
    
    int     Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    int     tracks_ctr;
    
    
    double  pot, pot_total;

    // my objects
    std::vector<PandoraNuTrack>     tracks;
    std::vector<hit>                hits;
    
    
    //Module labels to get data products
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
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
    
    // output csv file of tracks
    ofstream tracks_file;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::ErezTracksAnalyzerNewMCTruthMatching::ErezTracksAnalyzerNewMCTruthMatching(fhicl::ParameterSet const & p):EDAnalyzer(p){
    Debug(0,"ErezTracksAnalyzerNewMCTruthMatching::ErezTracksAnalyzerNewMCTruthMatching();");
    reconfigure(p);
    Debug(3,"reconfigure(p);");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::analyze(art::Event const & evt){
    
    Debug(2,"ub::ErezTracksAnalyzerNewMCTruthMatching::analyze(art::Event const & evt)");
    ResetVars();
    Debug(3,"ResetVars()");
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //    art::ServiceHandle<cheat::BackTracker> bt;
    isdata = evt.isRealData();
    bool isMC = !evt.isRealData();
    
    // * run data
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    Debug(3,"run = evt.run(); subrun = evt.subRun(); event = evt.id().event()");
    
    // * hits
    Debug(3,"// * hits");
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
    
    // * tracks
    Debug(3,"// * tracks");
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
    
    // * associations
    Debug(3,"// * associations");
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    Debug(3,"art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);");


    // * new truth matching from /uboone/app/users/wketchum/dev_areas/mcc8_4_drop/gallery_macros/TruthMatchTracks.C
    auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    auto const& trk_vec(*trk_handle);
    
    std::cout << "\tThere are " << trk_vec.size() << " tracks in this event." << std::endl;
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, fTrackModuleLabel);

    // * MCTruth information
    // get the particles from the event
    art::Handle<std::vector<simb::MCParticle>> pHandle;
    evt.getByLabel(fG4ModuleLabel, pHandle);
    art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
    
    
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
    Debug(3,"after for (int i = 0; i < Nhits_stored ; ++i)");
    


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
        
        
        // new MC-truth matching
        bool DoNewMCtruthMatching = true;
        if (isMC && DoNewMCtruthMatching){
            std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i);
            Debug(3,Form( "\tThere are %d associated hits." ,(int)trk_hits_ptrs.size() ));
            
            std::unordered_map<int,double> trkide;
            double maxe=-1, tote=0;
            art::Ptr< simb::MCParticle > maxp_me; //pointer for the particle match we will calculate
            
            Debug(3,"art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit();");
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle , evt , fHitParticleAssnsModuleLabel);
            std::vector<art::Ptr<simb::MCParticle>> particle_vec;
            std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
            
            
            Debug(3,"before for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h)");
            if (debug>3) SHOW(trk_hits_ptrs.size());
            //loop only over our hits
            for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h){
                
                // for each hit, ask how many particles match this hit
                Debug(4,Form("i_h: %d, particle_vec.clear(); match_vec.clear();",(int)i_h));
                particle_vec.clear(); match_vec.clear();
                Debug(4,"particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);");
                particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection
                Debug(4,Form("There are %d particles matched to hit %d" , (int)particle_vec.size() ,(int)i_h ));
                
                if (particle_vec.size()>0){
                    Debug(5,Form("%d",(int)particle_vec.size()));
                    //loop over particles
                    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                        Debug(5,Form("i_p: %d, particle_vec[i_p]->TrackId(): %d...",(int)i_p,(int)particle_vec[i_p]->TrackId()));
                        float Edep_particle = match_vec[i_p]->energy;  // energy deposited by ionization by this track ID [MeV]
                        float Edep_particle_frac = match_vec[i_p]->ideFraction; // fraction of energy in hit from this particle
                        Debug(5, Form("Edep_particle:%f, Edep_particle_frac:%f",Edep_particle,Edep_particle_frac));
                        
                        trkide[ particle_vec[i_p]->TrackId() ] += Edep_particle; //store energy per track id
                        tote += Edep_particle; //calculate total energy deposited
                        
                        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                            maxe = trkide[ particle_vec[i_p]->TrackId() ];
                            maxp_me = particle_vec[i_p];
                        }
                    }//end loop over particles per hit
                }
                Debug(4,"after if (particle_vec.size()>0)");
            }
            // We want to add a completeness variable,
            // defined as the ratio of the energy of the particle that contributed most to this track to the total true energy associated to this trackID
            
            // Now have matched the truth information - plug into the track object
            const art::Ptr< simb::MCParticle > particle = maxp_me;
            if (!particle.isNull()){
                track.SetMCpdgCode( particle->PdgCode() );
                track.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                track.SetTruthDirection();
                track.SetTruthLength();
                track.SetTruthMomentum( particle -> Momentum() );
                track.SetTruthMother( particle -> Mother() );
                track.SetTruthProcess( particle -> Process() );
                
                // * MC-truth information
                // To this end, we need to back-track the MCTruth/MCParticle association.
                // we do this with art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
                // but pHandle must get the correct particle Key,
                // which, since we defined it as an art::Ptr, is obtained from key() method
                int particle_key = (int)particle.key();
                if (debug>5){
                    SHOW(particle.key());
                    SHOW( particle->PdgCode() );
                    SHOW3( particle->Vx() , particle->Vy() , particle->Vz() );
                }
                if( fo.isValid() ){
                    Debug(5,"fo.isValid()");
                    auto pHandle_at_particle_key = pHandle->at(particle_key);
                    if (debug>5){
                        Printf("particle handle at particle_key %d",particle_key);
                        SHOW( pHandle_at_particle_key.PdgCode() );
                        SHOW3( pHandle_at_particle_key.Vx(),pHandle_at_particle_key.Vy(),pHandle_at_particle_key.Vz() );
                    }
                    art::Ptr<simb::MCTruth> mc_truth = fo.at(particle_key);
                    if (mc_truth->Origin() == simb::kBeamNeutrino)      track.SetOrigin( "beam neutrino" );
                    else if (mc_truth->Origin() == simb::kCosmicRay)    track.SetOrigin( "cosmic ray" );
                    SHOW( mc_truth -> Origin() );
                }// end if fo.isValid()
                else {
                    Debug(5,"fo.isValid() is false");
                }
            }
        }
        tracks.push_back( track );
    }
    Debug(3,"after for (int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i )");
    
    fTree -> Fill();
    StreamTracksToCSV();
    PrintInformation();
    Debug(3,"PrintInformation()");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::HeaderTracksInCSV(){
    
    tracks_ctr = 0;
    
    tracks_file
    << "run"        << "," << "subrun" << "," << "event" << ","
    << "track_id"   << ",";
    
    // reconstructed features
    tracks_file
    << "start_x"<< "," << "start_y" << "," << "start_z" << ","
    << "end_x"  << "," << "end_y"   << "," << "end_z"   << ","
    << "theta"  << "," << "phi"     << ","
    << "PIDa"   << "," << "length"  << ",";
    
    // truth information
    tracks_file
    << "pdg" << ","
    << "origin";
    
    // finish header
    tracks_file << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::StreamTracksToCSV(){

    // whatever you add here - must add also in header
    // i.e. in
    // ub::ErezTracksAnalyzerNewMCTruthMatching::HeaderTracksInCSV()
    for (auto t:tracks){
        
        tracks_ctr++;
        
        tracks_file
        << t.GetRun()       << "," << t.GetSubrun() << "," << t.GetEvent() << ","
        << t.GetTrackID()   << ",";
        
        // reconstructed features
        tracks_file
        << t.GetStartPos().x()  << "," << t.GetStartPos().y()   << "," << t.GetStartPos().z() << ","
        << t.GetEndPos().x()    << "," << t.GetEndPos().y()     << "," << t.GetEndPos().z() << ","
        << t.GetTheta()         << "," << t.GetPhi()            << ","
        << t.GetPIDa()          << "," << t.GetLength()         << ",";
        
        // truth information
        tracks_file
        << t.GetMCpdgCode() << ","
        << t.GetOrigin();
        
        // finish track features
        tracks_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::PrintInformation(){
    
    PrintXLine();
    Printf( "processed so far %d events", (int)(fTree->GetEntries()) );
    SHOW3( run , subrun , event );
    
    if(!tracks.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << tracks.size() << " pandoraNu tracks\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto t: tracks) {
            t.Print( true );
        }
    } else {cout << "\033[33m" << "xxxxxxxxxxxxxx\n" << "no reco tracks\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;}
    

    // time stamp
    PrintLine();
    end_ana_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
    std::cout << "\033[33m"
    << "finished analysis of this event at " << std::ctime(&end_time)
    << "time elapsed: " << elapsed_seconds.count() << "s"
    << "\033[31m" << endl;

    cout << "wrote " << tracks_ctr << " tracks to output file " << endl;
    EndEventBlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::beginJob(){
    
    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of events");
    fTree->Branch("run"     ,&run           ,"run/I");
    fTree->Branch("subrun"  ,&subrun        ,"subrun/I");
    fTree->Branch("event"   ,&event         ,"event/I");
    fTree->Branch("Ntracks" ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("Nhits"   ,&Nhits_stored  ,"Nhits/I");
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("hits"                ,&hits);

    
    fPOTTree = tfs->make<TTree>("pottree","pot tree");
    fPOTTree->Branch("pot",&pot,"pot/D");
    fPOTTree->Branch("run",&run,"run/I");
    fPOTTree->Branch("subrun",&subrun,"subrun/I");

    // output csv file
    tracks_file.open(fDataSampleLabel+"_tracks.csv");
    cout << "opened tracks file: "+fDataSampleLabel+"_tracks" << endl;
    HeaderTracksInCSV();
    
    pot_total = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::endSubRun(const art::SubRun& sr){
    
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
void ub::ErezTracksAnalyzerNewMCTruthMatching::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel        = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel    = p.get< std::string >("GenieGenModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
    fDataSampleLabel        = p.get< std::string >("DataSampleLabel");
    fPOTModuleLabel         = p.get< std::string >("POTModuleLabel");
    debug = p.get< int >("VerbosityLevel");
    fHitParticleAssnsModuleLabel = p.get< std::string >("HitParticleAssnsModuleLabel");
    fG4ModuleLabel          = p.get< std::string >("G4ModuleLabel","largeant");

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezTracksAnalyzerNewMCTruthMatching::ResetVars(){
    
    MCmode = false;
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = 0 ;
    pot = 0;
    tracks.clear();
    hits.clear();
    
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezTracksAnalyzerNewMCTruthMatching)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -