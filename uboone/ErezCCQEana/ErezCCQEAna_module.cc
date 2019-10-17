////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEAna
// Plugin Type: analyzer (art v2_05_00)
// File:        ErezCCQEAna_module.cc
//
// Generated at Wed May 2 2018 by Erez Cohen as a copy of "ErezCCQEAnalyzerNewTruthMatching_module.cc"
// The main differences are (during creation)
// (1) cleaning of the code,
// (2) adding Marco' flash-matching from NeutrinoFlashMatch_module.cc
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
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/PxUtils.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// for SwT emulation
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
// for MC truth matching
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
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"
#include "uboone/ErezCCQEana/MyObjects/pairVertex.h"
#include "uboone/ErezCCQEana/MyObjects/TruncMean.h"
//// flash matching from Marco
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
// MC event weight
#include "uboone/EventWeight/MCEventWeight.h"


// constants
constexpr int debug          = 1;
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTruth      = 100;
constexpr int kMaxNgenie     = 100;
constexpr float kMaxInterTrackDistance = 11; // 11 cm between tracks - maximal distance for clustering
constexpr float kMin_dQ_inTruthMatchedHits = 0.1; // the minimal fraction of hits-charge for MC-truth matching
constexpr float EPSILON      = 0.1;   // tollerance for equations

// charge deposition around the vertex in a box of N(wires) x N(time-ticks)
constexpr int N_box_sizes    = 30;
constexpr int MinNwiresBox   = 5;
constexpr int dNwiresBox     = 5;
constexpr int MinNticksBox   = 10;
constexpr int dNticksBox     = 10;

constexpr int N_r_around_vertex = 40;
constexpr float Min_r_around_vertex = 0;
constexpr float dr_around_vertex  = 0.5;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class ErezCCQEAna; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::ErezCCQEAna : public art::EDAnalyzer {
public:
    explicit ErezCCQEAna(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEAna(ErezCCQEAna const &) = delete;
    ErezCCQEAna(ErezCCQEAna &&) = delete;
    ErezCCQEAna & operator = (ErezCCQEAna const &) = delete;
    ErezCCQEAna & operator = (ErezCCQEAna &&) = delete;
    
    // Required functions.
    void                        analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                       beginJob () override;
    void                         endJob () override;
    void                    reconfigure (fhicl::ParameterSet const& p) override;
    void                      endSubRun (const art::SubRun& sr);

    
    // functionallity
    void                 CollectGlobals (art::Event const & evt);
    void                    CollectHits (art::Event const & evt);
    void                 CollectFlashes (art::Event const & evt);
    void                  CollectTracks (art::Event const & evt);
    void           CollectMCinformation (art::Event const & evt);
    void            CollectEventWeights (art::Event const & evt);
    void                  RunFlashMatch (art::Event const & evt);
    void                   OpenCSVFiles ();
    void                    StreamToCSV ();
    void                   PrintAndFill ();
    
    void              ConstructVertices ();
    void        ClusterTracksToVertices ();
    void                AnalyzeVertices ();
    void         FilterGoodPairVertices ();
    void                    TagVertices ();
    void               PrintInformation (bool DoPrintTracksFull=false);
    bool         TrackAlreadyInVertices (int ftrack_id);
    void              HeaderEventsInCSV ();
    void              StreamEventsToCSV ();
    void            HeaderVerticesInCSV ();
    void            StreamVerticesToCSV ();
    void              HeaderTracksInCSV ();
    void              StreamTracksToCSV ();
    void               HeaderGENIEInCSV ();
    void               StreamGENIEToCSV ();

    bool ParticleAlreadyMatchedInThisHit ( std::vector<int> ,int );
    double    GetRdQInSphereAroundVertex ( pairVertex v, int plane, float r);

    flashana::QCluster_t     GetQCluster (std::vector<art::Ptr<recob::Track>>);
    void                GetFlashLocation (std::vector<double>, double&, double&, double&, double&);

    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=5;
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

    Int_t run, subrun, event;
    
    short   isdata;
    
    bool    CSVfilesAlreadyOpened=false;
    bool    DoAnalyze=true;
    bool    MCmode=false;
    bool    DoFillOutputTree=false; // fill output-tree with all the events (consumes time and resources)
    bool    DoWriteVerticesInformation=false; // save also all the vertices information to a csv file
    bool    DoWriteTracksInformation=false; // save also all the tracks information to a csv file
    bool    DoAddTracksEdep=false; // add tracks dE/dx information
    bool    DoWriteGENIEInformation=false; // save also all the genie information to a csv file
    bool    DoWriteEventsInformation=false; // save also all the genie information to a csv file
    bool    EventPassedSwTrigger=false;
    bool    DoOnlySwT=false;
    bool    DoDebugRSE=false; // run only a specific set of r/s/e given as a fhicl parameter for debugging
    
    int     Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    int     Nflashes;
    int     Nvertices;
    int     vertices_ctr=0, tracks_ctr=0, genie_interactions_ctr=0,events_ctr=0;
    int     CC1pVertices_ctr=0 , CosmicVertices_ctr=0;
    int     NwiresBox[N_box_sizes], NticksBox[N_box_sizes];
    float   r_around_vertex[N_r_around_vertex];
    float   TruncMeanRad;
    
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
    std::string fCalibCaloModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDataSampleLabel;
    std::string fPOTModuleLabel;
    std::string fHitParticleAssnsModuleLabel;
    std::string fG4ModuleLabel;
    std::string fSwTModuleLabel;
    std::string fSwTAlgoModuleLabel;
    std::string fPIDModuleLabel;
    std::string fCaliPIDModuleLabel;
    std::string fTPCobjectProducer;
    std::string fgenieEventWeightModuleLabel;
    std::string fluxEventWeightModuleLabel;
    
    //mctruth information
    Int_t    mcevts_truth;    //number of neutrino Int_teractions in the spill
    
    // time stamp
    std::chrono::time_point<std::chrono::system_clock> start_ana_time, end_ana_time;
    
    // output csv file of vertices
    ofstream vertices_file, summary_file, tracks_file, genie_file, events_file;
    
    
    //flux information
    Int_t       ptype_flux;        //Parent GEANT code particle ID
    Int_t       tptype_flux;     //Type of parent particle leaving BNB/NuMI target
    Int_t       pntype_flux;      //oscillated neutrino type
    
    Float_t     pdpx_flux;        //Parent X momentum at decay point (GeV)
    Float_t     pdpy_flux;        //Parent Y momentum at decay point (GeV)
    Float_t     pdpz_flux;        //Parent Z momentum at decay point (GeV)
    Float_t     ppvx_flux;        //Parent production vertex X (cm)
    Float_t     ppvy_flux;        //Parent production vertex Y (cm)
    Float_t     ppvz_flux;        //Parent production vertex Z (cm)
    Float_t     pppz_flux;        //Parent Z momentum at production (GeV)
    Float_t     vx_flux;          //X position of hadron/muon decay (cm)
    Float_t     vy_flux;          //Y position of hadron/muon decay (cm)
    Float_t     vz_flux;          //Z position of hadron/muon decay (cm)
    Float_t     muparpx_flux;     //Muon neutrino parent production vertex X (cm)
    Float_t     muparpy_flux;     //Muon neutrino parent production vertex Y (cm)
    Float_t     muparpz_flux;     //Muon neutrino parent production vertex Z (cm)
    Float_t     mupare_flux;      //Muon neutrino parent energy (GeV)
    Float_t     tprivx_flux;      //Primary particle interaction vertex X (cm)
    Float_t     tprivy_flux;      //Primary particle interaction vertex Y (cm)
    Float_t     tprivz_flux;      //Primary particle interaction vertex Z (cm)
    
    
    // Marco' flash matching
    Float_t     FlashRangeStart, FlashRangeEnd;
    
    std::vector<::flashana::Flash_t>    beam_flashes;
    bool DoOpDetSwap;                 ///< If true swaps reconstructed OpDets according to _opdet_swap_map
    
    double                              FlashMatch_xfixed_chi2, FlashMatch_xfixed_ll;
    
    std::vector<int>                    OpDetSwapMap;    ///< The OpDet swap map for reco flashes
    std::vector<double>                 NumcFlashSpec;
    std::vector<double>                 BeamFlashSpec;
    std::vector<std::vector<double>>    HypothesisFlashSpec;
    std::vector<double>                 FlashMatch_score, FlashMatch_t0;
    std::vector<double>                 FlashMatch_qll_xmin, FlashMatch_tpc_xmin;
    std::vector<double>                 FlashMatch_xfixed_hypo_spec;
    std::vector<ubana::FlashMatch>      flashMatchTrackVector;

    ::flashana::FlashMatchManager       FlashMatch_mgr;
    std::vector<flashana::FlashMatch_t> FlashMatch_result;
    
    
    // event-weights
    std::vector<std::string>            event_weight_names;
    std::vector<float>                  event_weight_values;
    std::vector<std::vector<int>>       DebugRSEarray;
    
    // Handlers
    // * MCTruth information

    std::vector< art::Ptr<simb::MCTruth> >          mclist;
    std::vector< art::Ptr<simb::MCFlux> >           fluxlist;
    std::vector< art::Ptr<recob::Track> >           tracklist;
    std::vector< art::Ptr<recob::OpFlash> >         flashlist;
    std::vector< art::Ptr<recob::Hit> >             hitlist;

    
    art::ServiceHandle<geo::Geometry> geom;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::ErezCCQEAna::ErezCCQEAna(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::analyze(art::Event const & evt){
    
    ResetVars();
    CollectGlobals(evt);
    if (!DoAnalyze) return;
    
    CollectHits( evt );
    CollectFlashes( evt );
    std::cout<<"before CollectTracks in analyze "<<std::endl;
    CollectTracks( evt );
    std::cout<<"after CollectTracks in analyze "<<std::endl;
    
    // if this is an MC event, collect all MC information
    if (MCmode) {
        CollectMCinformation( evt );
        CollectEventWeights(evt);
    }
    if (DoWriteVerticesInformation) {
        ConstructVertices();
        RunFlashMatch(evt);
    }
    StreamToCSV();
    PrintAndFill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectGlobals(art::Event const & evt){
    
    Debug(4,"ub::ErezCCQEAna::CollectGlobals(art::Event const & evt)");
    events_ctr++;
    // global variabls (of the event) and the SwT (emulation)
    // * r/s/e
    isdata = evt.isRealData();
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    if (DoDebugRSE) {
        DoAnalyze=false;
        SHOW3(run,subrun,event);
        for (auto DebugRSE : DebugRSEarray) {
            SHOW3(DebugRSE.at(0),DebugRSE.at(1),DebugRSE.at(2));
            if (run==DebugRSE.at(0) && subrun==DebugRSE.at(1) && event==DebugRSE.at(2)){
                DoAnalyze = true;
            }
        }
    }
    // * SwT
    art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    evt.getByLabel(fSwTModuleLabel, softwareTriggerHandle);
    if (!softwareTriggerHandle.isValid()){ Debug(2,"SwT Handle doesn't exist!"); }
    if (softwareTriggerHandle.isValid()) {
        Debug(3,"Got SwT Handle!, number of algorithms: %",softwareTriggerHandle->getNumberOfAlgorithms());
        std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
        
        for (int i = 0; i < int(algoNames.size()); i++) {
            Debug(3,"algoNames[%] = %",i,algoNames[i]);
            // This is the algorithm that we need
            if (algoNames[i] == fSwTAlgoModuleLabel) {
                EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[0]) ? true : false;
            }
        }
    }
    if (DoOnlySwT && MCmode==false) {
        StreamEventsToCSV();
        DoAnalyze = false;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectHits(art::Event const & evt){
    Debug(2,"ub::ErezCCQEAna::CollectHits(art::Event & evt)");
    // ----------------------------------------
    // hits information
    // ----------------------------------------
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
        art::fill_ptr_vector(hitlist, hitListHandle);
    
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectFlashes(art::Event const & evt){
    Debug(2,"ub::ErezCCQEAna::CollectFlashes(art::Event & evt)");
    
    // ----------------------------------------
    // flash information
    // ----------------------------------------
    art::Handle< std::vector<recob::OpFlash> >      flashListHandle;
    if (evt.getByLabel(fFlashModuleLabel, flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);
    Nflashes = flashlist.size();
    for ( int f = 0; f < std::min(Nflashes,kMaxHits); f++ ) {
        
        if (debug>3){ SHOW2(flashlist[f]->Time() , flashlist[f]->TotalPE() ) };
        flash fflash(
                     flashlist[f]->Time(),         // flash time
                     flashlist[f]->TimeWidth(),    // flash time width
                     flashlist[f]->ZCenter(),      // flash z-center
                     flashlist[f]->ZWidth(),       // flash z-width
                     flashlist[f]->YCenter(),      // flash y-center
                     flashlist[f]->YWidth(),       // flash y-width
                     flashlist[f]->TotalPE()       // flash total PE
                     );
        // Katherine:
        // I chose a flash time between 0 and 10 as a very wide beam window that would definitely include anything beam related.
        // I use a much smaller window now similar to what the cc-inclusive analysis uses.
        // The thing to be careful about is that the different MC and data files have the beam window at different times.
        // See Ariana's slide 2 in the MCC8 validation slides. These are up to date.
        // The reason for a 6.5 P.E. cut is because that is the threshold for the online data trigger.
        // So, I didn't want to include MC with smaller flashes if they are thrown out in data.
        // I also wanted as low as possible of a flash threshold because single protons give a small amount of scintillation light.
        // I know that the cc-inclusive use much higher to get a higher purity. I would guess that you probably also want a higher P.E. cut.
        if (debug>3){ SHOW2(fflash.GetTime() , fflash.GetTotalPE()) };
        if( (0.0 < fflash.GetTime()) && (fflash.GetTime() < 10.0) && (6.5 < fflash.GetTotalPE()) ){
            flashes.push_back( fflash );
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectTracks(art::Event const & evt){
    
    // * tracks
    art::Handle< std::vector<recob::Track> >        trackListHandle;
    // * truth matching
    art::Handle< std::vector<simb::MCParticle> >    pHandle;
    art::Handle< std::vector<simb::MCTruth> >       mctruthListHandle;
    art::Handle< std::vector<simb::MCFlux> >        mcfluxListHandle;

    Debug(4,"ub::ErezCCQEAna::CollectTracks(art::Event & evt)");
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
        art::fill_ptr_vector(tracklist, trackListHandle);
    art::Handle< std::vector<recob::OpFlash> >      flashListHandle;

    
    // * associations
    Debug(4,"// * hits-tracks association");
    art::FindManyP<recob::Hit>          fmth(trackListHandle, evt, fTrackModuleLabel);
    //std::cout<<"after setting fmcalipid"<<std::endl;
    // delete by Sep-30
    //    // uncalibrated calorimetry
    //    Debug(4,"// * uncalibrated calorimetry");
    //    art::FindMany<anab::Calorimetry>    fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    // calibrated calorimetry
        Debug(4,"// * calibrated calorimetry");
        art::FindMany<anab::Calorimetry>    fmcalical(trackListHandle, evt, fCalibCaloModuleLabel);
    // PIDa from PandoraNu
    Debug(4,"// * calibrated PID");
    art::FindMany<anab::ParticleID>     fmcalipid(trackListHandle, evt, fCaliPIDModuleLabel );
    //    Debug(4,"// * PIDa from PandoraNu");
    //    art::FindMany<anab::ParticleID>     fmpid(trackListHandle, evt, fPIDModuleLabel );
        // PIDaCali from PandoraNu
    
    

    auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, fTrackModuleLabel);
    

    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    trkf::TrackMomentumCalculator TrackMomCalc;

    // MCParticle from largeant
    int Nparticles=0;
    art::Handle< std::vector<simb::MCParticle> > MCParticleListHandle;
    std::vector<art::Ptr<simb::MCParticle> > mcparticlelist;
    if (MCmode){
        Debug(2,"// * MCParticle information");
        mcparticlelist.clear();
        if (evt.getByLabel("largeant",MCParticleListHandle))
            art::fill_ptr_vector(mcparticlelist, MCParticleListHandle);
        Nparticles = (int)mcparticlelist.size();
        Debug(2,"Nparticles: %",Nparticles);
        evt.getByLabel(fG4ModuleLabel, pHandle);
        art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
        
        mclist.clear();
        if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
            art::fill_ptr_vector(mclist, mctruthListHandle);
        
        if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
            art::fill_ptr_vector(fluxlist, mcfluxListHandle);
    }

    // MCS fit
    art::ValidHandle< std::vector <recob::MCSFitResult> > MCSMuHandle = evt.getValidHandle < std::vector <recob::MCSFitResult> > ("pandoraNuMCSMu");

    
    // ----------------------------------------
    // step through the tracks to collect tracks information
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
        
        track.SetPandoraTheta( tracklist[i]->Theta() );
        track.SetPandoraPhi( tracklist[i]->Phi() );
        const recob::MCSFitResult& mcsMu = MCSMuHandle->at(i);
        Debug(4 , "particle ID hypothesis used for MCS fit: %",mcsMu.particleIdHyp());
        track.SetMomentumMCS_fwd( mcsMu.fwdMomentum() );
        track.SetMomentumMCS_bwd( mcsMu.bwdMomentum() );
        
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
        
        //        // calorimetric KE from calibrated calorimetery
        //        if (fmcalical.isValid()){
        //            unsigned maxnumhits = 0;
        //            std::vector<const anab::Calorimetry*> calicalos = fmcalical.at(i);
        //            for (auto const& calicalo : calicalos){
        //                if (calicalo->PlaneID().isValid){
        //                    int plane = calicalo->PlaneID().Plane;
        //                    if (plane<0||plane>2) continue;
        //
        //                    // get the calorimetric kinetic energy of the track
        //                    track.SetCaloKEPerPlane( plane , calicalo->KineticEnergy() );
        //
        //                    // select the best plane as the one with the maximal number of charge deposition points
        //                    if (calicalo->dEdx().size() > maxnumhits){
        //                        maxnumhits = calicalo->dEdx().size();
        //                        track.SetBestPlaneCali ( plane );
        //                        track.SetMaxNHitsCali ( maxnumhits );
        //                    }
        //                }
        //            }
        //        }
        
        // PIDaCali and chi2 for particles hypotheses from PandoraNuCali
        if (fmcalipid.isValid()) {
            std::vector<const anab::ParticleID*> pids = fmcalipid.at(i);
            for (auto const& pid : pids){
                if (pid->PlaneID().isValid){
                    int plane = pid->PlaneID().Plane;
                    if (plane<0||plane>2) continue;
                    
                    track.SetPandoraNuCaliPID( plane
                                              , pid->Pdg()
                                              , pid->MinChi2()
                                              , pid->Chi2Proton()
                                              , pid->Chi2Kaon()
                                              , pid->Chi2Pion()
                                              , pid->Chi2Muon()
                                              , pid->PIDA()
                                              );
                    Debug(3,"% PIDa (plane %) for track %: %",fCaliPIDModuleLabel,plane,track.GetTrackID(),track.GetCaliPID_PIDA(plane));
                }
        
            }
	}
	else { std::cout<<"fmcalipid is NOT Valid "<<std::endl;}

	if (fmcalical.isValid()) {
            //adi adding dQdx dEdx and res range
            unsigned maxnumhits = 0;
            std::vector<const anab::Calorimetry*> calos = fmcalical.at(i);
            for (auto const& calo : calos){
                   if (calo->PlaneID().isValid){
                       int plane = calo->PlaneID().Plane;

                       // get the calorimetric kinetic energy of the track
                       track.SetCaloKEPerPlane( plane , calo->KineticEnergy() );
                       // select the best plane as the one with the maximal number of charge deposition points
                       if (calo->dEdx().size() > maxnumhits){
                       if (debug>3) SHOW3(track.GetTrackID(),plane,calo->dEdx().size());
                       if (debug>3) SHOW3(track.GetTrackID(),plane,calo->dQdx().size());
                       maxnumhits = calo->dEdx().size();
                       track.SetBestPlane ( plane );
                       track.SetMaxNHits ( maxnumhits );
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
                           track.FilldQdxHit( plane , calo->dQdx()[ip] );
                       }
                       if (used_trkres) pida /= used_trkres;
                       track.SetPIDaPerPlane( plane , pida );
                   }
               }
               track.SetPIDa();
               Debug(3 , "track.GetPIDa(): %",track.GetPIDa());
           }

               TruncMean truncmean( debug , TruncMeanRad );
               //std::cout<<"TruncMean (radius "<<TruncMeanRad<<")"<<std::endl;
               for (int plane=0; plane<3; plane++){
                   std::vector<float> trkdedx = track.GetdEdx(plane);
                   std::vector<float> trkresrg = track.GetResRange(plane);
                   std::vector<float> trkdedx_truncated;
//                   
                   Debug(3, "TruncMean (dE/dx) for track % plane % (radius=%), track.GetdEdx(%).size(): % \n -----------------"
                         , track.GetTrackID(), plane, TruncMeanRad ,plane,track.GetdEdx(plane).size());
//                   
                   truncmean.CalcTruncMean( trkresrg , trkdedx , trkdedx_truncated );
                   track.SetdEdxTrunc( plane , trkdedx_truncated );
                   if (debug>5) {
                       SHOWstdVector( track.GetdEdx(plane) );
                       SHOWstdVector( track.GetdQdx(plane) );
                       SHOWstdVector( track.GetdEdxTrunc(plane) );
                   }
               }
       }
       else {std::cout<<"fmcali not valid"<<std::endl;}
        
        // momentum from momentum calculator
        track.SetMomCalc(
                         TrackMomCalc.GetTrackMomentum(track.GetLength(),13),  // momentum for muon hypothesis
                         TrackMomCalc.GetTrackMomentum(track.GetLength(),2212) // momentum for proton hypothesis
                         );
        
        // flash - matching: find the closest flash to the track
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
        
        // MC-truth mathching for the tracks
        Debug(4,"before if (MCmode=%) on track %",MCmode,track.GetTrackID());
        if (MCmode){
            
            Debug(5,"inside (MCmode)");
            
            evt.getByLabel(fG4ModuleLabel, pHandle);
            art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
            
            
            std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i);
            Debug(6,"\tThere are % associated hits to track %." ,(int)trk_hits_ptrs.size(), track.GetTrackID());
            
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
                Debug(6,"--------------------------------");
                float dQinHit = trk_hits_ptrs[i_h]->Integral();
                dQinAllHits += dQinHit;
                Debug(6,"in track %, hit %, dQinAllHits: %",track.GetTrackID(),(int)i_h,dQinAllHits);
                // for each hit, ask how many particles match this hit
                particle_vec.clear(); match_vec.clear();
                particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
                //the .key() gives us the index in the original collection
                Debug(6,"in track %, There are % particles matched to hit %" ,track.GetTrackID(), (int)particle_vec.size() ,(int)i_h );
                
                // to avoid from matching the same particle more than once
                // we introduce a vector of matched TrackId-s
                // and require that the matched particle has not been mathced already for this hit
                std::vector <int> ParticlesMatchedInThisHit;
                if (particle_vec.size()>0){
                    //loop over particles that match this hit and ask which one deposited the most energy
                    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                        Debug(7,"i_p: %, particle_vec[i_p]->TrackId(): %...",(int)i_p,(int)particle_vec[i_p]->TrackId());
                        float Edep_particle = match_vec[i_p]->energy;  // energy deposited by ionization by this track ID [MeV]
                        float Edep_particle_frac = match_vec[i_p]->ideFraction; // fraction of energy in hit from this particle
                        Debug(7, "Edep_particle:%, Edep_particle_frac:%",Edep_particle,Edep_particle_frac);
                        
                        trkide[ particle_vec[i_p]->TrackId() ] += Edep_particle; //store energy [MeV] deposited by track id
                        
                        tote += Edep_particle; //calculate total energy deposited
                        if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
                            maxe = trkide[ particle_vec[i_p]->TrackId() ];
                            maxp_me = particle_vec[i_p];
                            Debug(7,"changing maxp_me to particle_vec[%], pdg:%, deposited maxe:%f",(int)i_p,maxp_me->PdgCode(),maxe);
                            
                            if (!ParticleAlreadyMatchedInThisHit( ParticlesMatchedInThisHit , (int)particle_vec[i_p]->TrackId()) ) {
                                ParticlesMatchedInThisHit.push_back( (int)particle_vec[i_p]->TrackId() );
                                trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ] += dQinHit; // store the integral on the hit by the track id
                                Debug(7, "dQinHit:%, trackid_dQinTruthMatchedHits[%]:%",dQinHit,particle_vec[i_p]->TrackId(),trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ]);
                                max_dQinTruthMatchedHits = trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ];
                                Debug(7,"max_dQinTruthMatchedHits:%, dQinAllHits:%",max_dQinTruthMatchedHits,dQinAllHits);
                            } else {
                                Debug(7,"particle of TrackID % was already matched in this hit",(int)particle_vec[i_p]->TrackId());
                            }
                        }
                    }//end loop over particles per hit
                }
                Debug(6,"after if (particle_vec.size()>0); maxe=%",maxe);
                purity = (fabs(tote)>0) ? maxe/tote : 0.0;
            }
            
            const art::Ptr< simb::MCParticle > particle = maxp_me;
            // Now have matched the truth information - plug into the track object
            // first, check completeness of the track MC-truth matching:
            // if a particle is to be matched with a track, require that it deposited more than some fraction kMin_dQ_inTruthMatchedHits in the charge deposited in the track
            // keep track of the ratio max_dQinTruthMatchedHits / dQinAllHits,
            // in order to select the best value to cut on...
            track.SetMaxdQinTruthMatchedHits(max_dQinTruthMatchedHits);
            track.SetdQinAllHits(dQinAllHits);
            if (debug>5) SHOW3( max_dQinTruthMatchedHits , dQinAllHits, track.GetRatiodQinTruthMatchedHits() );
            if ( (!particle.isNull()) && (track.GetRatiodQinTruthMatchedHits() >= kMin_dQ_inTruthMatchedHits )) {
                if (debug>5) {
                    cout << "matched track " << track.GetTrackID() << " with particle " << particle->PdgCode() << endl;
                    SHOW2( particle -> Mother() , particle -> Process() );
                }
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
                
                
                // from which particle did it come in the largeant list
                for( int i_mcparticle = 0; (i_mcparticle < std::min(Nparticles,kMaxTruth)) ; i_mcparticle++ ){
                    art::Ptr<simb::MCParticle> mcparticle = mcparticlelist.at(i_mcparticle);
                    if (mcparticle == particle) {
                        track.SetTruthMCParticle( i_mcparticle );
                    }
                }
                
            }else {
                Debug(6,"particle.isNull() = true!...");
                if (!maxp_me.isNull()) {
                    Debug(6,  "Un-Matching particle % since it deposited too few of hit-charge in track (max_dQinTruthMatchedHits/dQinAllHits = %)" , maxp_me->PdgCode() , max_dQinTruthMatchedHits / dQinAllHits );
                }
            }
            
            
        }//MC
        tracks.push_back( track );
    }
    return;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectMCinformation(art::Event const & evt){
    
    Debug(3,"ub::ErezCCQEAna::CollectMCinformation()");
    art::Handle< std::vector<simb::MCParticle> >    pHandle;
    art::Handle< std::vector<simb::MCTruth> >       mctruthListHandle;
    art::Handle< std::vector<simb::MCFlux> >        mcfluxListHandle;
    art::Handle< std::vector<recob::Track> >        trackListHandle;

    
    // MCParticle from largeant
    Debug(2,"// * MCParticle information");
    art::Handle< std::vector<simb::MCParticle> > MCParticleListHandle;
    std::vector<art::Ptr<simb::MCParticle> > mcparticlelist;    
    if (evt.getByLabel("largeant",MCParticleListHandle))
        art::fill_ptr_vector(mcparticlelist, MCParticleListHandle);
    auto Nparticles = (int)mcparticlelist.size();
    Debug(2,"Nparticles: %",Nparticles);
    evt.getByLabel(fG4ModuleLabel, pHandle);
    art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
    
    mclist.clear();
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
            art::fill_ptr_vector(mclist, mctruthListHandle);
        
    if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
            art::fill_ptr_vector(fluxlist, mcfluxListHandle);
    
    
    

    // ----------------------------------------
    // MC information
    // ----------------------------------------
    // ----------------------------------------
    // MC flux information
    // ----------------------------------------
    Debug(5,"fluxlist size: %",(int)fluxlist.size());
    if (fluxlist.size()){
        ptype_flux  = fluxlist[0]->fptype;
        pdpx_flux   = fluxlist[0]->fpdpx;
        pdpy_flux   = fluxlist[0]->fpdpy;
        pdpz_flux   = fluxlist[0]->fpdpz;
        pntype_flux = fluxlist[0]->fntype;
        //  position of neutrino interaction point in the beamline
        vx_flux     = fluxlist[0]->fvx;
        vy_flux     = fluxlist[0]->fvy;
        vz_flux     = fluxlist[0]->fvz;
        if (debug>5) {
            SHOW4(ptype_flux,pdpx_flux,pdpy_flux,pdpz_flux);
            SHOW4(pntype_flux,vx_flux,vy_flux,vz_flux);
        }
    }
    
    // ----------------------------------------
    // MC truth information
    // ----------------------------------------
    mcevts_truth = mclist.size();
    Debug(4,"before if (mcevts_truth), mclist.size(): %",mclist.size());
    if (mcevts_truth){
        MCmode = true;
        Debug(4,"MCmode = true;");
        
        for( int mc_evend_id = 0; (mc_evend_id < mcevts_truth) && (mc_evend_id < kMaxTruth) ; mc_evend_id++ ){
            art::Ptr<simb::MCTruth> mctruth = mclist[mc_evend_id];
            Debug(5,"art::Ptr<simb::MCTruth> mctruth = mclist[%];",mc_evend_id);
            if (mctruth->Origin() == simb::kBeamNeutrino){
                Debug(5,"mctruth->Origin() == simb::kBeamNeutrino");
                
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
                    Debug(5,"inside if ( NgenieParticles )");
                    for( int iPart = 0; iPart < std::min( NgenieParticles , kMaxNgenie ); iPart++ ){
                        const simb::MCParticle& part( mctruth->GetParticle(iPart) );
                        // add a primary to genie-interaction
                        genie_interaction.AddPrimary( part.PdgCode()    // pdg code
                                                     ,part.Momentum()   // 4-momentum
                                                     ,part.StatusCode() // status code
                                                     ,part.Mother()     // mother
                                                     ,part.Process()    // process
                                                     ) ;
                        if (part.StatusCode()==1) {
                            // try to match the primary particle with a track
                            for (auto & track : tracks){
                                if (debug>6){
                                    cout << "in for (auto & track : tracks) of iPart " << iPart << ", which is a particle with PDG code " << part.PdgCode() << endl;
                                    SHOW2(track.GetTrackID(),track.GetMCpdgCode());
                                    SHOW2(part.Momentum().P(),track.GetTruthMomentum().P());
                                    SHOW2(part.Momentum().M(),track.GetTruthMomentum().M());
                                }
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
                                    Debug(6,"adding track % to genie interaction % as a particle %",track.GetTrackID(),mc_evend_id,part.PdgCode());
                                }
                            } // for (auto & track : tracks)
                        } // end if (part.StatusCode()==1)
                    } // for particle
                }
                genie_interaction.SortNucleons();
                genie_interaction.ComputePmissPrec();
                genie_interaction.SetTruthTopology();
                genie_interaction.SetReconstructedTopology();
                
                // Jan-28, 2018: adding the neutrino origin in the beam coordinates and interaction point the detector system (from Zarko)
                if ((int)fluxlist.size() >= (int)mc_evend_id) {
                    //  position of neutrino interaction point in the beamline
                    genie_interaction.SetNuIntInBeam(TVector3(fluxlist[mc_evend_id]->fvx
                                                              ,fluxlist[mc_evend_id]->fvy
                                                              ,fluxlist[mc_evend_id]->fvz));
                    
                    
                }
                
                genie_interactions.push_back( genie_interaction );
            }//mctruth->Origin()
        } //  for( int mc_evend_id = 0; (mc_evend_id < mcevts_truth))
    }//is neutrino
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::ErezCCQEAna::ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs
                                         ,int cTrackID ){
    // to avoid from matching the same particle more than once
    // we introduce a vector of matched TrackId-s for each hit
    // and require that the matched particle has not been mathced already for this hit
    if (debug>5) {
        cout << "ub::ErezCCQEAna::ParticleAlreadyMatchedInThisHit()" << endl;
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
void ub::ErezCCQEAna::ConstructVertices(){
    Debug(4,"ub::ErezCCQEAna::ConstructVertices()");
    // cluster all tracks at close proximity to vertices
    ClusterTracksToVertices();
    // analyze these vertices: inter-tracks distances, angles...
    AnalyzeVertices();
    // retain only vertices with pairs of 2-tracks at close proximity
    FilterGoodPairVertices();
    // if its a MC event, tag the vertex by their MC information
    if (MCmode) TagVertices();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::ClusterTracksToVertices(){
    
    // Jan-18, 2018
    // cluster all tracks at close proximity to vertices, and fix the position of each vertex
    bool    FoundCloseTracks=false;
    float   closest_distance_ij;
    TVector3 vertex_position;
    
    for (int i=0; i < Ntracks; i++){
        if (!tracks[i].IsTrackContainedSoft()) continue; // 3<x<257, -115<y<115, 5<z<1045
        pairVertex vertex( run, subrun, event , vertices.size() );
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
        if (FoundCloseTracks) {
            Debug( 5 , "if (MCmode)");
            if (MCmode){ // if this is an MC event, match a GENIE interaction to the vertex
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
                        vertex.SetGENIEinfo( genie_interactions.at( mc_id_t0 ) );
                        if (debug>6){
                            Debug( 8 , "vertex.SetGENIEinfo( genie_interactions.at( mc_id_t0 ) );");
                            Printf("matched GENIE information");
                            vertex.GetGENIEinfo().Print(); // PRINTOUT
                        }
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
bool ub::ErezCCQEAna::TrackAlreadyInVertices(int ftrack_id){
    for (auto v:vertices){
        if ( v.IncludesTrack( ftrack_id ) ) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::AnalyzeVertices(){
    if (vertices.size()>0){
        for (auto & v:vertices){
            // after fixing the vertext position, remove far tracks
            Debug( 4 , "v.RemoveFarTracks( kMaxInterTrackDistance );");
            if (v.GetTracks().size()<2) continue;
            
            v.RemoveFarTracks( kMaxInterTrackDistance );
            
            // now sort the tracks
            Debug( 4 ,"SortTracksByPIDA;" );
            v.SortTracksByPIDA ();
            
            Debug( 4 , "SortTracksByChi2Proton;" );
            v.SortTracksByChi2Proton ();
            
            
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
void ub::ErezCCQEAna::FilterGoodPairVertices(){
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    
    // for momentum reconstuction from range / otherwise
    trkf::TrackMomentumCalculator trkm;
    
    // this funciton also kiils all vertices which are not up to standard.
    // to this end, we create a temporary vertices array and clear all vertices
    std::vector<pairVertex> tmp_vertices = vertices;
    vertices.clear();
    
    Debug(3 , "ub::ErezCCQEAna::FilterGoodPairVertices()");
    for (auto & v:tmp_vertices) {
        // vertices with only two tracks at close proximity and nothing else
        if (
            // we are looking for clusters of only two tracks that are fully contained
            v.GetNtracks() == 2
            // if there is a semi-contained track, with start/end point too close to the vertex, we don't want the vertex...
            &&  v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size() == 0
            ){
            
            // assign muon and proton tracks by chi2-proton
            auto Track_muCandidate = v.GetLargeChi2ProtonTrack();
            auto PmuFromRange = trkm.GetTrackMomentum( Track_muCandidate.GetLength()  , 13  );
            v.AssignMuonTrack( Track_muCandidate  );
            
            auto Track_pCandidate = v.GetSmallChi2ProtonTrack();
            auto PpFromRange = trkm.GetTrackMomentum( Track_pCandidate.GetLength() , 2212);
            v.AssignProtonTrack( Track_pCandidate );
            
            Debug(4,"set muon and proton candidates: track % is muon and track % is proton candidate",Track_muCandidate.GetTrackID(),Track_pCandidate.GetTrackID());
            
            v.FixTracksDirections ();
            // after fixing the direction of the muon track we set its MCS momentum.
            // (the direction fixation is needed in order to determine
            // forward or backward MCS fit.)
            v.SetMCSMuMomentum ( Track_muCandidate.GetMomentumMCS() );
            // now we can set the reco. kinematics
            // which also include mcs-based kinematics
            v.SetReconstructedFeatures ( PmuFromRange , PpFromRange );
            
            // set vertex position in the three wire planes
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
            // plug hits associated with the proton and the muon tracks
            v.AssociateHitsToTracks( hits );
            
            vertices.push_back( v );
        }
        // if the vertex is not with only two tracks at close proximity and nothing else, kill it
        else { Debug( 3 , Form("erasing vertex %d, since N(tracks)=%d, N(close semi-contained tracks):%d"
                               , v.GetVertexID(), v.GetNtracks(), (int)(v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size()) ) );}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::TagVertices(){
    
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
        
        PandoraNuTrack t1 = v.GetTrack_muCandidate();
        PandoraNuTrack t2 = v.GetTrack_pCandidate();
        
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
                ){
                if (v.GetClosestGENIE().GetIsCC_1p_200MeVc()==true) {    // the closest GENIE is a CC1p
                    v.SetAsCC1p();
                    if (v.GetClosestGENIE().GetIsCC_1p_200MeVc_0pi()==true) {    // the closest GENIE is a CC1p0π
                        v.SetAsCC1p0pi();
                    }
                }
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
double ub::ErezCCQEAna::GetRdQInSphereAroundVertex(pairVertex v, int plane, float r){
    Debug(3,"ub::ErezCCQEAna::GetRdQInSphereAroundVertex(vertex %, plane %, radius % cm)",v.GetVertexID(),plane,r);
    // Feb-18,2018
    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a sphere of radius r [cm] around the vertex in plane i=0,1,2
    // a mal-function will return -1 (if plane is not 0,1,2) or -9999 (if Qtotal=0)
    // input:
    // plane, r [cm] for the sphere, hits in event
    //
    // This method is implemented here (and not in pairVertex.cxx)
    // since we need lar::util::PxPoint,
    // which can not be easily implemented outside a LArSoft module
    auto const& geom = ::util::GeometryUtilities::GetME();
    
    float Qtotal=0.0, Qmuon=0.0, Qproton=0.0;
    // vertex PxPoint in this plane
    //    util::PxPoint *vPxPoint = new util::PxPoint( plane , v.GetVertexWire(plane), v.GetVertexTime(plane) );
    Debug(3,"vertex wire: %, time: % \n ------ --- -- -- -- --- -- -- -- ----- -",v.GetVertexWire(plane), v.GetVertexTime(plane));
    
    // PxPoint of all hits in this plane
    for (auto hit: hits) {
        if (hit.GetPlane()==plane) {
            //            util::PxPoint *hitPxPoint = new util::PxPoint( hit.GetPlane(), hit.GetWire(), hit.GetPeakT() );
            //            double distance_hit_vertex = geom->Get2DDistance( hitPxPoint , vPxPoint );
            double distance_hit_vertex = geom->Get2DDistance( v.GetVertexWire(plane), v.GetVertexTime(plane) , hit.GetWire(), hit.GetPeakT() );
            // check if the hit & the vertex are close enough
            if (distance_hit_vertex < r) {
                // if yes, aggregate the charge
                Qtotal += hit.GetCharge();
            }
            if (debug>7) {
                hit.Print();
                cout << "distance to vertex: " << distance_hit_vertex << " cm" << endl;
                if (distance_hit_vertex < r) {
                    Debug(0, "added hit charget to Qtotal: %",Qtotal);
                }
            }
        }
    }
    // PxPoint of the vertex-muon hits in this plane
    for (auto hit: v.GetMuonHits(plane)) {
        if (hit.GetPlane()==plane) {
            
            //            util::PxPoint *hitPxPoint = new util::PxPoint( hit.GetPlane(), hit.GetWire(), hit.GetPeakT() );
            //            double distance_hit_vertex = geom->Get2DDistance( hitPxPoint , vPxPoint );
            double distance_hit_vertex = geom->Get2DDistance( v.GetVertexWire(plane), v.GetVertexTime(plane) , hit.GetWire(), hit.GetPeakT() );
            if (distance_hit_vertex < r) {
                Qmuon += hit.GetCharge();
            }
        }
    }
    // PxPoint of the vertex-proton hits in this plane
    for (auto hit: v.GetProtonHits(plane)) {
        if (hit.GetPlane()==plane) {
            
            //            util::PxPoint *hitPxPoint = new util::PxPoint( hit.GetPlane(), hit.GetWire(), hit.GetPeakT() );
            //            double distance_hit_vertex = geom->Get2DDistance( hitPxPoint , vPxPoint );
            double distance_hit_vertex = geom->Get2DDistance( v.GetVertexWire(plane), v.GetVertexTime(plane) , hit.GetWire(), hit.GetPeakT() );
            
            if (distance_hit_vertex < r) {
                Qproton += hit.GetCharge();
            }
        }
    }
    Debug(3,"completed ub::ErezCCQEAna::GetRdQInSphereAroundVertex(), Qtotal: %, Qmuon: %, Qproton: %",Qtotal,Qmuon,Qproton);
    if (Qmuon>Qtotal || Qproton>Qtotal) {
        Debug(1,"Qtotal=% < Qmuon=% or Qproton=% !? ",Qtotal,Qmuon,Qproton);
    }
    if (fabs(Qtotal)>0) return ((Qmuon+Qproton)/Qtotal);
    else return -1;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::RunFlashMatch(art::Event const & evt){
    Debug(5, "ub::ErezCCQEAna::RunFlashMatch()");
    // ----------------------------------------
    // Neutrino Flash match from Marco
    //
    // We emplace the tracks-vector from
    // each candidate vertex (muon and proton candidates),
    // and create a vector of matchings
    // to the candidate vertices.
    // This way we can ask for example about
    // the matching score for each vertex,
    // the distance of the matched flash
    // to that vertex, the N(PE) etc.
    // ----------------------------------------
    art::Handle< std::vector<recob::OpFlash> >      flashListHandle;
    evt.getByLabel(fFlashModuleLabel, flashListHandle);
    Debug(5, "flashListHandle->size() is %",flashListHandle->size());
    int nBeamFlashes = 0;
    for (size_t n = 0; n < flashListHandle->size(); n++) {
        auto const& flash = (*flashListHandle)[n];
        Debug( 2 , "[NeutrinoFlashMatch] Flash time from % : %",  fFlashModuleLabel,flash.Time());
        // keep only flashes in the veto time range
        if(flash.Time() < FlashRangeStart || FlashRangeEnd < flash.Time()) continue;
        nBeamFlashes++;
        
        // Construct a Flash_t
        ::flashana::Flash_t f;
        f.x = f.x_err = 0;
        f.pe_v.resize(geom->NOpDets());
        f.pe_err_v.resize(geom->NOpDets());
        for (unsigned int i = 0; i < f.pe_v.size(); i++) {
            unsigned int opdet = geom->OpDetFromOpChannel(i);
            if (DoOpDetSwap && evt.isRealData()) {
                opdet = OpDetSwapMap.at(opdet);
            }
            f.pe_v[opdet] = flash.PE(i);
            f.pe_err_v[opdet] = sqrt(flash.PE(i));
        }
        double Ycenter, Zcenter, Ywidth, Zwidth;
        GetFlashLocation(f.pe_v, Ycenter, Zcenter, Ywidth, Zwidth);
        f.y = Ycenter;
        f.z = Zcenter;
        f.y_err = Ywidth;
        f.z_err = Zwidth;
        f.time = flash.Time();
        f.idx = nBeamFlashes-1;
        beam_flashes.resize(nBeamFlashes);
        beam_flashes[nBeamFlashes-1] = f;
    } // flash loop
    // If more than one beam flash, take the one with more PEs
    if (nBeamFlashes > 1) {
        Debug( 4, "More than one beam flash in this event. Taking beam flash with more PEs.");
        // Sort flashes by length
        std::sort(beam_flashes.begin(), beam_flashes.end(),
                  [](::flashana::Flash_t a, ::flashana::Flash_t b) -> bool
                  {return a.TotalPE() > b.TotalPE();});
    }
    Debug(3,"nBeamFlashes: %",nBeamFlashes);
    // Emplace the 'best' flash to Flash Matching Manager
    if (nBeamFlashes==0) {
        Debug(4, "no beam flashes, returning...");
        return;
    }
    ::flashana::Flash_t f = beam_flashes[0];
    BeamFlashSpec.resize(f.pe_v.size());
    BeamFlashSpec = f.pe_v;
    Debug(4,"f.pe_v.size(): %",f.pe_v.size());
    FlashMatch_mgr.Emplace(std::move(f));
    Debug(5, "Done Neutrino Flash match from Marco");
    
    
    // match
    // -------
    // The class method GetQCluster(std::vector<art::Ptr<recob::Track>>)
    // takes in input a vector of tracks and returns a QCluster_t,
    // which is the data product needed by the flash matching
    // In Marco' original module (NeutrinoFlashMatch) this is done by passing TPCobject to the manager
    // which is a product of the entire UBXsec chain.
    // Here,
    // since we do not want to impose the UBXsec chain analysiss
    // we work aroud this by giving the pandoraNu track-list (tracklist) as an input to the manager
    Debug(4,"match: len(vertices): %",vertices.size());
    for (auto & v: vertices) {
        // Build the QCluster for this vertex
        flashana::QCluster_t qcluster;
        // to get the right tracks for this vertex we use the track-id's for this event
        art::Ptr<recob::Track> recobTrack_muCandidate = tracklist[v.GetTrack_muCandidate().GetTrackID()];
        art::Ptr<recob::Track> recobTrack_pCandidate = tracklist[v.GetTrack_pCandidate().GetTrackID()];
        std::vector<art::Ptr<recob::Track> > vertex_tracks = {recobTrack_muCandidate,recobTrack_pCandidate};
        qcluster = this->GetQCluster( vertex_tracks );
        FlashMatch_mgr.Emplace(std::move(qcluster));
    }
    
    // Now the matching is done by the manager 1 time
    // and it scores the matching for all objects
    // (candidate pairs in our case)
    Debug(4,"FlashMatch_mgr.Match!...");
    FlashMatch_result = FlashMatch_mgr.Match();
    
    // save the results
    // i.e. the matched flashs
    // to each of the vertices
    // ------------------------
    Debug(4,"save the results: FlashMatch_result: size %",FlashMatch_result.size());
    HypothesisFlashSpec.resize( FlashMatch_result.size() );
    for( int _matchid=0; _matchid < (int)(FlashMatch_result.size()); ++_matchid) {
        
        auto const& match = FlashMatch_result[_matchid];
        Debug(4,"match.flash_id: %, match.score: %", match.flash_id, match.score);
        
        FlashMatch_score.push_back( match.score );
        
        auto const& matched_flash = FlashMatch_mgr.FlashArray()[match.flash_id];
        FlashMatch_t0.push_back( matched_flash.time );
        
        double Ycenter, Zcenter, Ywidth, Zwidth;
        GetFlashLocation(matched_flash.pe_v, Ycenter, Zcenter, Ywidth, Zwidth);
        
        flash matched_fflash(
                             matched_flash.time,    // flash time
                             matched_flash.time,    // flash time width, that I do not know, and is unimportant here
                             Zcenter,      // flash z-center
                             Zwidth,       // flash z-width
                             Ycenter,      // flash y-center
                             Ywidth,       // flash y-width
                             matched_flash.TotalPE()       // flash total PE
                             );
        auto & v = vertices.at(_matchid);
        v.SetMatchedFlash( matched_fflash , match.score );
        if (debug>4){
            matched_fflash.Print();
            Debug(4,"t0[%]: %, HypothesisFlashSpec[%].size(): %", _matchid, matched_flash.time, _matchid, HypothesisFlashSpec[_matchid].size());
            for(size_t pmt=0; pmt < HypothesisFlashSpec[_matchid].size(); ++pmt) {
                Debug(4,"match.hypothesis[pmt %]: %", pmt, match.hypothesis[pmt]);
            };
        }
    }
    // ----------------------------------------
    // end Marco' flash matching
    // ----------------------------------------
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::OpenCSVFiles(){
    // output csv files
    Debug( 2, Form( "DoWriteEventsInformation: %d",DoWriteEventsInformation) );
    if (DoWriteEventsInformation) {
        events_file.open(fDataSampleLabel+"_events.csv");
        cout << "opened events file: "+fDataSampleLabel+"_events.csv" << endl;
        HeaderEventsInCSV();
    }
    
    
    Debug( 2, Form( "DoWriteGNEIEInformation: %d",DoWriteGENIEInformation) );
    if (DoWriteGENIEInformation) {
        genie_file.open(fDataSampleLabel+"_genie.csv");
        cout << "opened genie file: "+fDataSampleLabel+"_genie.csv" << endl;
        HeaderGENIEInCSV();
    }
    
    Debug( 2, Form( "DoWriteTracksInformation: %d",DoWriteTracksInformation) );
    if (DoWriteTracksInformation) {
        tracks_file.open(fDataSampleLabel+"_tracks.csv");
        cout << "opened tracks file: "+fDataSampleLabel+"_tracks.csv" << endl;
        HeaderTracksInCSV();
    }
    //    vertices_file.open("/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+fDataSampleLabel+"_vertices.csv");
    Debug( 2, Form( "DoWriteVerticesInformation: %d",DoWriteVerticesInformation) );
    if (DoWriteVerticesInformation) {
        vertices_file.open(fDataSampleLabel+"_vertices.csv");
        cout << "opened vertices file: "+fDataSampleLabel+"_vertices.csv" << endl;
        HeaderVerticesInCSV();
    }
    CSVfilesAlreadyOpened = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::StreamToCSV(){
    
    // ----------------------------------------
    // write information to CSV files
    // ----------------------------------------
    // if this is the first event - open and write a header into the csv files
    if (!CSVfilesAlreadyOpened) {
        OpenCSVFiles();
    }
    
    Debug(4,"ub::ErezCCQEAna::StreamToCSV()");
    if (DoWriteEventsInformation) {
        StreamEventsToCSV();
        Debug(5,"StreamEventsToCSV()");
    }
    if (DoWriteTracksInformation){
        StreamTracksToCSV();
        Debug(5,"StreamTracksToCSV()");
    }
    if (DoWriteGENIEInformation){
        StreamGENIEToCSV();
        Debug(5,"StreamGENIEToCSV()");
    }
    // ----------------------------------------
    // write ccqe-candidate pairs to CSV files
    // ----------------------------------------
    if (DoWriteVerticesInformation){
        StreamVerticesToCSV();
        Debug(5,"StreamVerticesToCSV();");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::HeaderEventsInCSV(){
    
    events_ctr = 0;
    
    events_file
    << "run"        << "," << "subrun" << "," << "event" << ",";
    
    // SwT
    events_file
    << "passed_SwT" << ",";
    
    // neutrino interactions
    events_file
    << "N(v-interactions)" << ",";
    
    // first neutrino interaction
    events_file
    << "E(v-interaction)" << ","
    << "x(v-interaction)" << ","
    << "y(v-interaction)" << ","
    << "z(v-interaction)" ;
    
    
    // finish header
    events_file << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::StreamEventsToCSV(){
    
    // whatever you add here - must add also in header
    // i.e. in
    // ub::ErezCCQEAna::HeaderEventsInCSV()
    events_file
    << run        << "," << subrun << "," << event << ",";
    
    // SwT
    events_file
    << EventPassedSwTrigger << ",";
    
    // neutrino interactions
    if (!genie_interactions.empty()){
        events_file
        << genie_interactions.size() << ",";
        
        // first neutrino interaction
        events_file
        << genie_interactions.at(0).GetEv() << ","
        << genie_interactions.at(0).GetVertexPosition().x() << ","
        << genie_interactions.at(0).GetVertexPosition().y() << ","
        << genie_interactions.at(0).GetVertexPosition().z() ;
        
    } else {
        
        events_file
        << 0 << ",";
        
        // first neutrino interaction
        events_file
        << -999 << ","
        << -999 << ","
        << -999 << ","
        << -999 ;
    }
    
    
    // finish header
    events_file << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::HeaderGENIEInCSV(){
    Debug(0,"ub::ErezCCQEAna::HeaderGENIEInCSV()");
    genie_interactions_ctr = 0;
    
    genie_file
    << "run" << "," << "subrun" << "," << "event" << "," ;
    
    // truth topology
    genie_file
    << "IsCC_Np_200MeVc" << ","
    << "IsCC_1p_200MeVc" << ","
    << "IsCC_1p_200MeVc_0pi" << ","
    << "IsCCQE" << ",";
    
    // reconstructed topology
    genie_file
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
    genie_file
    << "truth_Pmu" << ","
    << "truth_Pmu_theta" << "," << "truth_Pmu_phi" << ","
    << "truth_Pmu_x" << "," << "truth_Pmu_y" << "," << "truth_Pmu_z" << ","
    << "truth_Pmu_cos_theta" << ","
    << "truth_Pp" << ","
    << "truth_Pp_theta" << "," << "truth_Pp_phi" << ","
    << "truth_Pp_x" << "," << "truth_Pp_y" << "," << "truth_Pp_z" << ","
    << "truth_Pp_cos_theta" << ",";
    
    // relevant truth-information
    genie_file
    << "truth_Ev" << ","
    << "truth_Q2" << ",";
    
    genie_file
    << "truth_Pv_x" << ","
    << "truth_Pv_y" << ","
    << "truth_Pv_z" << ","
    << "truth_Pv_theta" << ",";
    
    
    // only for 1mu-1p vertices
    genie_file
    << "reconstructed mu-p distance" << "," ;
    
    // v-interaction point in detector system
    genie_file
    << "truth_x" << ","
    << "truth_y" << ","
    << "truth_z" << ",";

    // event weight
    for (auto name: event_weight_names) {
        Debug(3,"for (auto name: event_weight_names), name=%",name);
        genie_file << name << ",";
    }
    
    // v-interaction point in beam coordinates system
    genie_file
    << "truth_x_beamCoordinates" << ","
    << "truth_y_beamCoordinates" << ","
    << "truth_z_beamCoordinates" << ",";
    
    
    // application of phase-space and kinematical cuts
    genie_file
    << "delta_phi"              << "," << "Pt" << "," << "theta_12";
    
    // finish
    genie_file << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::StreamGENIEToCSV(){
    // Jan-26, 2018
    // whatever you add here - must add also in header - ub::ErezCCQEGENIENewTruthMatching::HeaderVerticesInCSV()
    for (auto g:genie_interactions){
        
        genie_interactions_ctr++;
        
        genie_file
        << g.GetRun() << "," << g.GetSubrun() << "," << g.GetEvent() << "," ;
        
        
        // truth topology
        genie_file
        << g.GetIsCC_Np_200MeVc() << ","
        << g.GetIsCC_1p_200MeVc() << ","
        << g.GetIsCC_1p_200MeVc_0pi() << ","
        << g.GetIsCCQE() << ",";
        
        // reconstructed topology
        genie_file
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
        genie_file
        << g.GetPmu().P() << ","
        << g.GetPmu().Theta() << "," << g.GetPmu().Phi() << ","
        << g.GetPmu().Px() << "," << g.GetPmu().Py() << "," << g.GetPmu().Pz() << ","
        << TMath::Cos(g.GetPmu().Theta()) << ","
        << g.GetPp().P() << ","
        << g.GetPp().Theta() << "," << g.GetPp().Phi() << ","
        << g.GetPp().Px() << "," << g.GetPp().Py() << "," << g.GetPp().Pz() << ","
        << TMath::Cos(g.GetPp().Theta()) << ",";
        
        // relevant truth-information
        genie_file
        << g.GetEv() << ","
        << g.GetQ2() << ",";
        
        genie_file
        << g.GetPv().Px() << ","
        << g.GetPv().Py() << ","
        << g.GetPv().Pz() << ","
        << g.GetPv().Theta() << ",";
        
        
        // only for 1mu-1p vertices
        genie_file << g.GetReco_mu_p_distance() << ",";
        
        // v-interaction point in detector system
        genie_file
        << g.GetVertexPosition().x() << ","
        << g.GetVertexPosition().y() << ","
        << g.GetVertexPosition().z() << ",";
        
        // event weight
        for (auto weight: event_weight_values) {
            genie_file << weight << ",";
        }

        
        // v-interaction point in beam coordinates system
        genie_file
        << g.GetNuIntInBeam().x() << ","
        << g.GetNuIntInBeam().y() << ","
        << g.GetNuIntInBeam().z() << ",";
        

        // application of phase-space and kinematical cuts
        genie_file
        << r2d*(g.GetPp().Phi()-g.GetPmu().Phi())       << ","
        << (g.GetPp().Vect() + g.GetPmu().Vect()).Pt()  << ","
        << r2d*(g.GetPp().Vect().Angle(g.GetPmu().Vect()));

        // finish
        genie_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::HeaderTracksInCSV(){
    
    tracks_ctr = 0;
    
    tracks_file
    << "run"        << "," << "subrun" << "," << "event" << ","
    << "track_id"   << ",";
    
    // reconstructed features
    tracks_file
    << "start_x"<< "," << "start_y" << "," << "start_z" << ","
    << "end_x"  << "," << "end_y"   << "," << "end_z"   << ","
    << "theta"  << "," << "phi"     << "," << "length"  << ","
    << "PIDa"   << "," << "PIDaCali"<< ","
    << "pid_PIDaYplane"             << ","
    << "calipid_PIDaYplane"         << ","
    << "calipid_Chi2ProtonYplane"   << ","
    << "calipid_Chi2MuonYplane"     << "," ;
    
    // truth information
    tracks_file
    << "pdg" << ","
    << "origin" << ","
    << "purity" << ",";
    
    // completeness of the track MC-truth matching
    tracks_file
    << "max_dQinTruthMatchedHits" << ","
    << "dQinAllHits" << ","
    << "Ratio_max_dQinTruthMatchedHits_dQinAllHits" << ",";
    
    
    // dE/dx
    if (DoAddTracksEdep) {
        tracks_file
        << "ResRange_0" << ","
        << "dEdx_0" << ","  //  << ","
        << "dQdx_0" << ","  //  << ","
        << "ResRange_1" << ","
        << "dEdx_1" << "," //  << ","
        << "dQdx_1" << "," //  << ","
        << "ResRange_2" << ","
        << "dEdx_2" << "," //  << ","
        << "dQdx_2" ; //  << ","
        //        << "Truncated_dEdx_Y";
    }
    
    // finish header
    tracks_file << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::StreamTracksToCSV(){
    
    // whatever you add here - must add also in header
    // i.e. in
    // ub::ErezSimpleTracksAnalyzer::HeaderTracksInCSV()
    for (auto t:tracks){
        
        tracks_ctr++;
        
        tracks_file
        << t.GetRun()       << "," << t.GetSubrun() << "," << t.GetEvent() << ","
        << t.GetTrackID()   << ",";
        
        // reconstructed features
        tracks_file
        << t.GetStartPos().x()  << "," << t.GetStartPos().y()   << "," << t.GetStartPos().z()   << ","
        << t.GetEndPos().x()    << "," << t.GetEndPos().y()     << "," << t.GetEndPos().z()     << ","
        << t.GetTheta()         << "," << t.GetPhi()            << "," << t.GetLength()         << ","
        << t.GetPIDa()          << "," << t.GetPIDaCali()       << ","
        << t.GetPID_PIDA(2)     << ","
        << t.GetCaliPID_PIDA(2)         << ","
        << t.GetCaliPID_Chi2Proton(2)   << ","
        << t.GetCaliPID_Chi2Muon(2)     << ",";
        
        
        // truth information
        tracks_file
        << t.GetMCpdgCode() << ","
        << t.GetOrigin() << ","
        << t.GetTruthPurity() << ",";
        
        // completeness of the track MC-truth matching
        tracks_file
        << t.GetMaxdQinTruthMatchedHits() << ","
        << t.GetdQinAllHits() << ","
        << t.GetRatiodQinTruthMatchedHits() << ",";

        // dE/dx
        if (DoAddTracksEdep) {
	  for (unsigned int iplane = 0 ; iplane < 3 ; iplane++) {
  	    // residual range
            std::vector<float> ResRange_Y = t.GetResRange( iplane ); // collection plane (Y) = 2
            tracks_file << "\"[";
            for (auto ResRange_Y_hit:ResRange_Y) {
                tracks_file << ResRange_Y_hit << ",";
            }
            tracks_file << "]\"" << ",";
            
            // dE/dx
            std::vector<float> dEdx_Y = t.GetdEdx( iplane ); // collection plane (Y) = 2
            tracks_file << "\"[";
            for (auto dEdx_Y_hit:dEdx_Y) {
                tracks_file << dEdx_Y_hit << ",";
            }
            tracks_file << "]\"" << ",";
            
	    // dQ/dx
            std::vector<float> dQdx_Y = t.GetdQdx( iplane ); // collection plane (Y) = 2
            tracks_file << "\"[";
            for (auto dQdx_Y_hit:dQdx_Y) {
                tracks_file << dQdx_Y_hit << ",";
            }
            tracks_file << "]\"" << ",";
          }
            /*
  	    // residual range
            std::vector<float> ResRange_Y = t.GetResRange( 2 ); // collection plane (Y) = 2
            tracks_file << "\"[";
            for (auto ResRange_Y_hit:ResRange_Y) {
                tracks_file << ResRange_Y_hit << ",";
            }
            tracks_file << "]\"" << ",";
            
            // dE/dx
            std::vector<float> dEdx_Y = t.GetdEdx( 2 ); // collection plane (Y) = 2
            tracks_file << "\"[";
            for (auto dEdx_Y_hit:dEdx_Y) {
                tracks_file << dEdx_Y_hit << ",";
            }
            tracks_file << "]\"" << ",";
            */
//            // trunctated dE/dx
//            std::vector<float> dEdxTrunc_Y = t.GetdEdxTrunc( 2 ); // collection plane (Y) = 2
//            tracks_file << "\"[";
//            for (auto dEdx_Y_hit:dEdxTrunc_Y) {
//                tracks_file << dEdx_Y_hit << ",";
//            }
//            tracks_file << "]\"";
        }
        
        // finish track features
        tracks_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::HeaderVerticesInCSV(){
    
    vertices_ctr = 0;
    
    vertices_file
    << "run"    << "," << "subrun"  << "," << "event"   << "," << "vertex_id" << ","
    << "x"      << "," << "y"       << "," << "z"       << ","
    << "track_id" << ","
    
    // µ/p assigned tracks
    << "l_muCandidate"  << ","
    << "l_pCandidate"   << ","
    << "l_mu-l_p"       << ","
    
    // pandoraNu pid object
    << "pidcali_Chi2ProtonYplane_muCandidate"   << ","
    << "pidcali_Chi2ProtonYplane_pCandidate"    << ","
    << "pidcali_Chi2MuonYplane_muCandidate"     << ","
    << "pidcali_Chi2MuonYplane_pCandidate"      << ","
    
    // start/end points, for FV cuts
    << "startx_muCandidate" << ","  << "truth_startx_muCandidate"   << ","
    << "starty_muCandidate" << ","  << "truth_starty_muCandidate"   << ","
    << "startz_muCandidate" << ","  << "truth_startz_muCandidate"   << ","
    << "startx_pCandidate"  << ","  << "truth_startx_pCandidate"    << ","
    << "starty_pCandidate"  << ","  << "truth_starty_pCandidate"    << ","
    << "startz_pCandidate"  << ","  << "truth_startz_pCandidate"    << ","
    << "endx_muCandidate"   << ","  << "truth_endx_muCandidate"     << ","
    << "endy_muCandidate"   << ","  << "truth_endy_muCandidate"     << ","
    << "endz_muCandidate"   << ","  << "truth_endz_muCandidate"     << ","
    << "endx_pCandidate"    << ","  << "truth_endx_pCandidate"      << ","
    << "endy_pCandidate"    << ","  << "truth_endy_pCandidate"      << ","
    << "endz_pCandidate"    << ","  << "truth_endz_pCandidate"      << ",";
    

    // CC1p0π reconstructed featues
    vertices_file
    << "distance" << "," << "delta_phi" << "," << "delta_theta" << ","
    << "theta_12" << ","
    
    // reconstructed kinematics
    << "reco_Ev"    << ","  << "reco_Q2" << ","
    << "reco_Xb"    << ","  << "reco_y" << ","
    << "reco_W2"    << ","
    << "reco_Pt"    << ","  << "reco_theta_pq" << ","
    << "reco_Ev_mcs"<< ","  << "reco_Q2_mcs"<< "," << "reco_Pt_mcs"<< ","
    // muon from range
    << "reco_Emu"       << ","
    << "reco_Pmu"       << ","
    << "reco_Pmu_x"     << "," << "reco_Pmu_y" << "," << "reco_Pmu_z" << ","
    << "reco_Pmu_theta" << "," << "reco_Pmu_phi" << ","
    // muon from mcs
    << "reco_Emu_mcs"       << ","
    << "reco_Pmu_mcs"       << ","
    << "reco_Pmu_mcs_x"     << ","<< "reco_Pmu_mcs_y"   << "," << "reco_Pmu_mcs_z" << ","
    << "reco_Pmu_mcs_theta" << ","<< "reco_Pmu_mcs_phi" << ","
    << "reco_Pmu_cos_theta" << ","
    // proton
    << "reco_Ep"            << ","
    << "reco_Pp"            << ","
    << "reco_Pp_x"          << "," << "reco_Pp_y"   << "," << "reco_Pp_z" << ","
    << "reco_Pp_theta"      << "," << "reco_Pp_phi" << ","
    << "reco_Pp_cos_theta"  << ","
    // virtual boson
    << "reco_q"         << ","  << "reco_omega"   << ","
    
    
    // missing momentum (reconstructed struck neutron)
    << "reco_Pmiss"     << ","
    << "reco_Pmiss_x"   << ","  << "reco_Pmiss_y" << ","  << "reco_Pmiss_z" << ","
    << "reco_Pmiss_t"   << ",";
    
    // calorimetry
    if (DoAddTracksEdep) { 
	vertices_file
        << "trackmu_resRange_0" << "," 
        << "trackmu_dQdx_0"     << ","
        << "trackmu_dEdx_0"     << ","
        << "trackmu_resRange_1" << "," 
        << "trackmu_dQdx_1"     << ","
        << "trackmu_dEdx_1"     << ","
        << "trackmu_resRange_2" << "," 
        << "trackmu_dQdx_2"     << ","
        << "trackmu_dEdx_2"     << ","
        << "trackp_resRange_0"  << "," 
        << "trackp_dQdx_0"      << ","
        << "trackp_dEdx_0"      << ","
        << "trackp_resRange_1"  << "," 
        << "trackp_dQdx_1"      << ","
        << "trackp_dEdx_1"      << ","
        << "trackp_resRange_2"  << "," 
        << "trackp_dQdx_2"      << ","
        << "trackp_dEdx_2"      << ",";
    }
 
    // truth MC information
    vertices_file
    << "truth_l_muCandidate"        << "," << "truth_l_pCandidate" << ","
    << "truth_purity_muCandidate"   << "," << "truth_purity_pCandidate" << ","
    
    << "truth_Emu"          << ","
    << "truth_Pmu"          << ","
    << "truth_Pmu_x"        << "," << "truth_Pmu_y"     << "," << "truth_Pmu_z" << ","
    << "truth_Pmu_theta"    << "," << "truth_Pmu_phi"   << ","
    << "truth_Pmu_cos_theta"<< ","
    << "truth_Ep"       << ","
    << "truth_Pp"       << ","
    << "truth_Pp_x"     << "," << "truth_Pp_y" << "," << "truth_Pp_z" << ","
    << "truth_Pp_theta" << "," << "truth_Pp_phi" << ","
    << "truth_Pp_cos_theta" << ","
    

    
    // mathing genie interaction
    << "genie_distance" << ","
    << "truth_Ev"       << "," << "truth_Q2"    << ","
    << "truth_Xb"       << "," << "truth_y"     << ","
    << "truth_W2"       << ","
    << "truth_Pt"       << "," << "truth_theta_pq" << ","
    << "genie_mode"     << ","
    << "truth_q"        << "," << "truth_omega" << ","

    
    // closest genie (e.g. a proton was detected in a µp event, which is not the original proton in a CC interaction, since the real proton rescattered)
    << "closest_genie_distance" << ","
    << "closest_genie_Ev"       << "," << "closest_genie_Q2"    << ","
    << "closest_genie_Xb"       << "," << "closest_genie_y"     << ","
    << "closest_genie_W2"       << ","
    << "closest_genie_Pt"       << "," << "closest_genie_theta_pq" << ","
    << "closest_genie_mode"     << ","
    << "closest_genie_Nprotons" << ","
    << "closest_genie_Nneutrons"<< ","
    << "closest_genie_Npions"   << ","
    << "closest_genie_Npi0"     << ","
    << "closest_genie_CCNC"     << ","
    
    
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
    vertices_file
    << "Nflashes"                   << ","
    << "MatchedFlash_Ydistance"     << ","
    << "MatchedFlash_Zdistance"     << ","
    << "MatchedFlash_YZdistance"    << ","
    << "MatchedFlash_TotalPE"       << ","
    << "MatchedFlash_Time"          << ","
    << "MatchedFlash_Score"         << ",";

    
    
    // broken trajectories
    vertices_file
    << "isBrokenTrajectory" << ",";
    
    // event weight
    for (auto name: event_weight_names) {
        vertices_file << name << ",";
    }

    
    // vertex truth-topology in MC
    vertices_file
    << "1mu-1p" << "," << "CC1p" << "," << "other-pairs" << "," << "cosmic";
    
    // finish
    vertices_file << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::StreamVerticesToCSV(){
    // July-25, 2017
    // whatever you add here - must add also in header - ub::ErezCCQEAna::HeaderVerticesInCSV()
    for (auto v:vertices){
        auto trk_mu = v.GetTrack_muCandidate();
        auto trk_p = v.GetTrack_pCandidate();
        
        vertices_ctr++;
        if (v.GetIsGENIECC_1p_200MeVc()) { CC1pVertices_ctr++;}
        else if (v.GetIsCosmic()){ CosmicVertices_ctr++;}
        
        vertices_file
        << v.GetRun()           << "," << v.GetSubrun() << "," << v.GetEvent() << "," << v.GetVertexID() << ","
        << v.GetPosition().x()  << "," << v.GetPosition().y() << "," << v.GetPosition().z() << ",";
        
        // tracks id
        for ( auto track : v.GetTracks()){
            vertices_file << track.GetTrackID() << ";" ;
        }
        vertices_file << ",";
        
        
        // µ/p assigned tracks
        vertices_file
        << trk_mu.GetLength() << ","
        << trk_p.GetLength()  << ","
        << trk_mu.GetLength() - trk_p.GetLength() << ","
        
        // pandoraNu pid object
        << trk_mu.GetCaliPID_Chi2Proton(2)    << ","
        << trk_p.GetCaliPID_Chi2Proton(2)     << ","
        << trk_mu.GetCaliPID_Chi2Muon(2)      << ","
        << trk_p.GetCaliPID_Chi2Muon(2)       << ",";
        
        // start/end points, for FV cuts
        vertices_file
        << trk_mu.GetStartPos().x() << ","  << trk_mu.GetTruthStartPos().x()  << ","
        << trk_mu.GetStartPos().y() << ","  << trk_mu.GetTruthStartPos().y()  << ","
        << trk_mu.GetStartPos().z() << ","  << trk_mu.GetTruthStartPos().z()  << ","
        << trk_p.GetStartPos().x()  << ","  << trk_p.GetTruthStartPos().x()   << ","
        << trk_p.GetStartPos().y()  << ","  << trk_p.GetTruthStartPos().y()   << ","
        << trk_p.GetStartPos().z()  << ","  << trk_p.GetTruthStartPos().z()   << ","
        << trk_mu.GetEndPos().x()   << ","  << trk_mu.GetTruthEndPos().x()    << ","
        << trk_mu.GetEndPos().y()   << ","  << trk_mu.GetTruthEndPos().y()    << ","
        << trk_mu.GetEndPos().z()   << ","  << trk_mu.GetTruthEndPos().z()    << ","
        << trk_p.GetEndPos().x()    << ","  << trk_p.GetTruthEndPos().x()     << ","
        << trk_p.GetEndPos().y()    << ","  << trk_p.GetTruthEndPos().y()     << ","
        << trk_p.GetEndPos().z()    << ","  << trk_p.GetTruthEndPos().z()     << ",";
        
        // CC1p0π reconstructed featues
        // distances between the tracks
        for ( auto d : v.Get_distances_ij()) vertices_file << d ; //  for more than 2 tracks add << ";" ;
        vertices_file << ",";
        // delta_phi
        for ( auto delta_phi : v.Get_delta_phi_ij()) vertices_file << delta_phi ; // in degrees
        vertices_file << ",";
        // delta_theta
        for ( auto delta_theta : v.Get_delta_theta_ij()) vertices_file << delta_theta ; // in degrees
        vertices_file << ",";
        // theta_12 (the 3D angle between the two tracks)
        vertices_file << v.GetAngleBetween2tracks() << ","; // in degrees
        
        // reconstructed kinematics
        vertices_file
        << v.GetRecoEv()    << "," << v.GetRecoQ2()    << ","
        << v.GetRecoXb()    << "," << v.GetRecoY()     << ","
        << v.GetRecoW2()    << ","
        << v.GetRecoPt()    << "," << v.GetReco_theta_pq()  << ","
        << v.GetRecoEv_mcs()<< "," << v.GetRecoQ2_mcs()     << "," << v.GetRecoPt_mcs() << ",";
        // muon from range
        vertices_file
        << v.GetRecoPmu().E()       << ","
        << v.GetRecoPmu().P()       << ","
        << v.GetRecoPmu().Px()      << "," << v.GetRecoPmu().Py()   << "," << v.GetRecoPmu().Pz() << ","
        << v.GetRecoPmu().Theta()   << "," << v.GetRecoPmu().Phi()  << ","
        // muon from mcs
        << v.GetRecoPmu_mcs().E()       << ","
        << v.GetRecoPmu_mcs().P()       << ","
        << v.GetRecoPmu_mcs().Px()      << "," << v.GetRecoPmu_mcs().Py()   << "," << v.GetRecoPmu_mcs().Pz() << ","
        << v.GetRecoPmu_mcs().Theta()   << "," << v.GetRecoPmu_mcs().Phi()  << ","
        << TMath::Cos(v.GetRecoPmu_mcs().Theta()) << ",";
        
        // proton
        vertices_file
        << v.GetRecoPp().E()    << ","
        << v.GetRecoPp().P()    << ","
        << v.GetRecoPp().Px()   << "," << v.GetRecoPp().Py()    << "," << v.GetRecoPp().Pz() << ","
        << v.GetRecoPp().Theta()<< "," << v.GetRecoPp().Phi()   << ","
        << TMath::Cos(v.GetRecoPp().Theta()) << ",";
        
        // virtual boson
        vertices_file
        << v.GetReco_q()        << "," << v.GetReco_omega() << ",";
        
        // missing momentum (reconstructed struck neutron)
        vertices_file
        << v.GetRecoPmiss().P()     << ","
        << v.GetRecoPmiss().Px()    << ","  << v.GetRecoPmiss().Py() << ","  << v.GetRecoPmiss().Pz() << ","
        << v.GetRecoPmiss().Pt()    << ",";
      
        // calorimetry
        if (DoAddTracksEdep) { 
         for (unsigned int iplane = 0 ; iplane < 3 ; iplane++) {

	     std::vector<float> muResRange = v.GetTrack_muCandidate().GetResRange( iplane ); // mu res range 
             vertices_file << "\"[";
             for (auto ResRange_hit:muResRange) {
                vertices_file << ResRange_hit << ",";
             }
             vertices_file << "]\"" << ",";

	     std::vector<float> mudQdx = v.GetTrack_muCandidate().GetdQdx( iplane ); // mu dQdx 
             vertices_file << "\"[";
             for (auto dQdx_hit:mudQdx) {
                vertices_file << dQdx_hit << ",";
             }
             vertices_file << "]\"" << ",";
	     
	     std::vector<float> mudEdx = v.GetTrack_muCandidate().GetdEdx( iplane ); // mu dEdx 
             vertices_file << "\"[";
             for (auto dEdx_hit:mudEdx) {
                vertices_file << dEdx_hit << ",";
             }
             vertices_file << "]\"" << ",";
         }
         for (unsigned int iplane = 0 ; iplane < 3 ; iplane++) {
	     
	     std::vector<float> pResRange = v.GetTrack_pCandidate().GetResRange( iplane ); // p res range 
             vertices_file << "\"[";
             for (auto ResRange_hit:pResRange) {
                vertices_file << ResRange_hit << ",";
             }
             vertices_file << "]\"" << ",";

	     std::vector<float> pdQdx = v.GetTrack_pCandidate().GetdQdx( iplane ); // p dQdx 
             vertices_file << "\"[";
             for (auto dQdx_hit:pdQdx) {
                vertices_file << dQdx_hit << ",";
             }
             vertices_file << "]\"" << ",";

	     std::vector<float> pdEdx = v.GetTrack_pCandidate().GetdEdx( iplane ); // p dEdx 
             vertices_file << "\"[";
             for (auto dEdx_hit:pdEdx) {
                vertices_file << dEdx_hit << ",";
             }
             vertices_file << "]\"" << ",";
         }
        }
 
        // truth MC information
        vertices_file << trk_mu.GetTruthLength() << "," << trk_p.GetTruthLength() << ",";
        vertices_file << trk_mu.GetTruthPurity() << "," << trk_p.GetTruthPurity() << ",";
        
        vertices_file
        << trk_mu.GetTruthMomentum().E()    << ","
        << trk_mu.GetTruthMomentum().P()    << ","
        << trk_mu.GetTruthMomentum().Px()   << ","
        << trk_mu.GetTruthMomentum().Py()   << ","
        << trk_mu.GetTruthMomentum().Pz()   << ","
        << trk_mu.GetTruthMomentum().Theta()<< ","
        << trk_mu.GetTruthMomentum().Phi()  << ","
        << TMath::Cos(trk_mu.GetTruthMomentum().Theta())  << ",";
        
        vertices_file
        << trk_p.GetTruthMomentum().E()     << ","
        << trk_p.GetTruthMomentum().P()     << ","
        << trk_p.GetTruthMomentum().Px()    << ","
        << trk_p.GetTruthMomentum().Py()    << ","
        << trk_p.GetTruthMomentum().Pz()    << ","
        << trk_p.GetTruthMomentum().Theta() << ","
        << trk_p.GetTruthMomentum().Phi()   << ","
        << TMath::Cos(trk_p.GetTruthMomentum().Theta())  << ",";
        
        
        

        
        // mathing genie interaction
        vertices_file
        << v.GetDistanceToGENIE()       << ","
        << v.GetGENIEinfo().GetEv()     << "," << v.GetGENIEinfo().GetQ2()  << ","
        << v.GetGENIEinfo().GetXb()     << "," << v.GetGENIEinfo().GetY()   << ","
        << v.GetGENIEinfo().GetW2()     << ","
        << v.GetGENIEinfo().GetPt()     << "," << v.GetGENIEinfo().Get_theta_pq() << ","
        << v.GetGENIEinfo().GetMode()   << ","
        << v.GetGENIEinfo().Get_q()     << ","  << v.GetGENIEinfo().Get_omega() << ",";
        
        
        
        // closest genie (e.g. a proton was detected in a µp event, which is not the original proton in a CC interaction, since the real proton rescattered)
        vertices_file
        << v.GetDistanceToClosestGENIE()        << ","
        << v.GetClosestGENIE().GetEv()          << "," << v.GetClosestGENIE().GetQ2()       << ","
        << v.GetClosestGENIE().GetXb()          << "," << v.GetClosestGENIE().GetY()        << ","
        << v.GetClosestGENIE().GetW2()          << ","
        << v.GetClosestGENIE().GetPt()          << "," << v.GetClosestGENIE().Get_theta_pq()<< ","
        << v.GetClosestGENIE().GetMode()        << ","
        << v.GetClosestGENIE().GetNprotons()    << ","
        << v.GetClosestGENIE().GetNneutrons()   << ","
        << v.GetClosestGENIE().GetNpions()      << ","
        << v.GetClosestGENIE().GetNpi0()        << ","
        << v.GetClosestGENIE().GetCCNC()        << ",";
        

        
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
        vertices_file
        << flashes.size()                   << ","
        << v.GetYDis2MatchedFlash()         << ","
        << v.GetZDis2MatchedFlash()         << ","
        << v.GetDis2MatchedFlash()          << ","
        << v.GetMatchedFlash().GetTotalPE() << ","
        << v.GetMatchedFlash().GetTime()    << ","
        << v.GetMatchedFlashScore()         << ",";

        
        // broken trajectories
        vertices_file
        << v.GetIsBrokenTrajectory() << ",";

        
        // event weight
        for (auto weight: event_weight_values) {
            vertices_file << weight << ",";
        }
        
        // vertex truth-topology in MC
        vertices_file
        << v.GetIs1mu1p()               << ","
        << v.GetIsGENIECC_1p_200MeVc()  << ","
        << v.GetIsNon1mu1p()            << ","
        << v.GetIsCosmic();
        
        // finish
        vertices_file << endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::PrintAndFill(){
    Debug(4,"ub::ErezCCQEAna::PrintAndFill()");
    // ----------------------------------------
    // print and finish
    // ----------------------------------------
    if(DoFillOutputTree){
        fTree -> Fill();
    }
    if (debug>0) {
        PrintInformation( (debug>0) ? true : false );
    }
    if (debug>0) PrintXLine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::PrintInformation(bool DoPrintTracksFull){
    
    PrintXLine();
    SHOW3( run , subrun , event );
    
    if(MCmode && !genie_interactions.empty()){
        PrintHeader(Form("%d genie-interactions",(int)genie_interactions.size()));
        for (auto g: genie_interactions) {
            g.Print();
        }
    } else { PrintHeader("no genie-interactions"); }
    
    if(!flashes.empty()){
        PrintHeader(Form("%d flashes",(int)flashes.size()));
        for (auto f: flashes) {
            f.Print();
        }
    } else { PrintHeader("no flashes"); }
    
    if(!tracks.empty()){ PrintHeader(Form("%d pandoraNu-tracks",(int)tracks.size()));
        for (auto t: tracks) {
            if (DoPrintTracksFull) t.Print( true );
            else cout << t.GetTrackID() << "\t";
        }
        cout << endl;
    } else { PrintHeader("no pandoraNu-tracks");    }
    
    if(!vertices.empty()){
        PrintHeader(Form("%d pair-vertices",(int)vertices.size()));
        for (auto v: vertices) {
            v.Print((debug>2) ? true : false        // DoPrintTracks
                    ,(debug>0) ? true : false       // DoPrintFull
                    );
        }
    } else { PrintHeader("no pair-vertices"); }
    
    
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
void ub::ErezCCQEAna::beginJob(){
    
    // charge deposition around the vertex in a box of N(wires) x N(time-ticks)
    for (int i_box_size=0 ; i_box_size < N_box_sizes ; i_box_size++){
        NwiresBox[i_box_size] = MinNwiresBox + i_box_size * dNwiresBox;
        NticksBox[i_box_size] = MinNticksBox + i_box_size * dNticksBox;
    }
    // charge deposition around the vertex in a sphere of radius r [cm]
    for (int i_r_around_vertex=0 ; i_r_around_vertex < N_r_around_vertex ; i_r_around_vertex++){
        r_around_vertex[i_r_around_vertex] = Min_r_around_vertex + i_r_around_vertex * dr_around_vertex;
        Debug(10,"r_around_vertex[%]:%",i_r_around_vertex,r_around_vertex[i_r_around_vertex]);
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
    
    SHOW(debug);
    pot_total = 0;
    
    Printf("FlashMatch_mgr Configuration:");
    FlashMatch_mgr.PrintConfig();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::endJob(){
    Debug(3,"ub::ErezCCQEAna::endJob()");
    
    std::time_t now_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    summary_file.open(fDataSampleLabel+"_summary.csv");
    summary_file << "time" << ","
    << "POT" << ","
    << "Nevents" << ","
    << "Ntracks" << ","
    << "Ngenie_interactions" << ","
    << "Nvertices" << ","
    << "NCosmicVertices" << ","
    << "NCC1pVertices"
    << endl;
    
    std::string sTimeS = std::ctime(&now_time);
    
    summary_file << sTimeS.substr(0,sTimeS.length()-1) << ","
    << pot_total << ","
    << events_ctr << ","
    << tracks_ctr << ","
    << genie_interactions_ctr << ","
    << vertices_ctr << ","
    << CosmicVertices_ctr << ","
    << CC1pVertices_ctr
    << endl;
    
    summary_file.close();
    
    cout << "Ended job. Wrote summary into file: "+fDataSampleLabel+"_summary.csv" << endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::endSubRun(const art::SubRun& sr){
   

    art::Handle< sumdata::POTSummary > potListHandle;
    
    if(
       ((fPOTModuleLabel=="generator") && (sr.getByLabel(fPOTModuleLabel,potListHandle))) // for MC-BNB
        ||
       ((fPOTModuleLabel=="beamdata") && (sr.getByLabel(fPOTModuleLabel,"bnbETOR860",potListHandle))) // for BNB, Marco Del Tutto, Aug-27, 2017
        )
        pot = potListHandle->totpot - pot_total; // AA AP adding the POT bug fix 
        //pot = potListHandle->totpot;
    else
        pot = 0.;
    if (fPOTTree) fPOTTree->Fill();
    pot_total += pot;
    
    
    if (debug>3){
        Printf("end subrun %d",subrun);
        SHOW2( pot , pot_total );
    }
    Debug(4,"POT from this subrun: %" ,pot );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::ResetVars(){
    
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = Nvertices = 0 ;
    pot = 0;
    tracks.clear();
    hits.clear();
    flashes.clear();
    genie_interactions.clear();
    vertices.clear();
    EventPassedSwTrigger=false;
    DoAnalyze=true;
    BeamFlashSpec.clear();
    FlashMatch_result.clear();
    HypothesisFlashSpec.clear();
    NumcFlashSpec.clear();
    FlashMatch_score.clear();
    FlashMatch_t0.clear();
    FlashMatch_qll_xmin.clear();
    FlashMatch_tpc_xmin.clear();
    FlashMatch_xfixed_hypo_spec.clear();
    FlashMatch_mgr.Reset();
    flashMatchTrackVector.clear();
    event_weight_names.clear();
    event_weight_values.clear();

    mclist.clear();
    fluxlist.clear();
    flashlist.clear();
    tracklist.clear();
    hitlist.clear();
    
    NumcFlashSpec.clear();
    BeamFlashSpec.clear();
    HypothesisFlashSpec.clear();
    FlashMatch_score.clear();
    FlashMatch_t0.clear();
    FlashMatch_qll_xmin.clear();
    FlashMatch_tpc_xmin.clear();
    FlashMatch_xfixed_hypo_spec.clear();
    flashMatchTrackVector.clear();
    FlashMatch_result.clear();
    event_weight_values.clear();
    DebugRSEarray.clear();
    

    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::GetFlashLocation(std::vector<double> pePerOpDet,
                                          double& Ycenter,
                                          double& Zcenter,
                                          double& Ywidth,
                                          double& Zwidth)
{
    
    // Reset variables
    Ycenter = Zcenter = 0.;
    Ywidth  = Zwidth  = -999.;
    double totalPE = 0.;
    double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
    
    for (unsigned int opdet = 0; opdet < pePerOpDet.size(); opdet++) {
        
        if (opdet > 31 && opdet < 200){
            //  std::cout << "Ignoring channel " << opch << " as it's not a real channel" << std::endl;
            continue;
        }
        
        // Get physical detector location for this opChannel
        double PMTxyz[3];
        ::art::ServiceHandle<geo::Geometry> geo;
        geo->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
        
        // Add up the position, weighting with PEs
        sumy    += pePerOpDet[opdet]*PMTxyz[1];
        sumy2   += pePerOpDet[opdet]*PMTxyz[1]*PMTxyz[1];
        sumz    += pePerOpDet[opdet]*PMTxyz[2];
        sumz2   += pePerOpDet[opdet]*PMTxyz[2]*PMTxyz[2];
        
        totalPE += pePerOpDet[opdet];
    }
    
    Ycenter = sumy/totalPE;
    Zcenter = sumz/totalPE;
    
    // This is just sqrt(<x^2> - <x>^2)
    if ( (sumy2*totalPE - sumy*sumy) > 0. )
        Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;
    
    if ( (sumz2*totalPE - sumz*sumz) > 0. ) 
        Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
flashana::QCluster_t ub::ErezCCQEAna::GetQCluster(std::vector<art::Ptr<recob::Track>> track_v) {
    
    flashana::QCluster_t summed_qcluster;
    summed_qcluster.clear();
    
    for (unsigned int trk = 0; trk < track_v.size(); trk++) {
        
        art::Ptr<recob::Track> trk_ptr = track_v.at(trk);
        
        ::geoalgo::Trajectory track_geotrj;
        track_geotrj.resize(trk_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));
        
        for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
            auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
            track_geotrj[pt_idx][0] = pt[0];
            track_geotrj[pt_idx][1] = pt[1];
            track_geotrj[pt_idx][2] = pt[2];
        }
        
        auto qcluster = ((flashana::LightPath*)(FlashMatch_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(track_geotrj);
        summed_qcluster += qcluster;
        
    } // track loop
    
    return summed_qcluster;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel        = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel    = p.get< std::string >("GenieGenModuleLabel");
    fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel");
    fCalibCaloModuleLabel   = p.get< std::string >("CalibratedCalorimetryModuleLabel");
    
    fDataSampleLabel        = p.get< std::string >("DataSampleLabel");
    fPOTModuleLabel         = p.get< std::string >("POTModuleLabel");
    fFlashModuleLabel       = p.get< std::string >("FlashModuleLabel");
    debug                   = p.get< int >("VerbosityLevel");
    MCmode                  = p.get< bool >("MCmodeLabel",false);
    fHitParticleAssnsModuleLabel = p.get< std::string >("HitParticleAssnsModuleLabel");
    fG4ModuleLabel          = p.get< std::string >("G4ModuleLabel","largeant");
    DoFillOutputTree        = p.get< bool >("DoFillOutputTree",false);
    DoWriteTracksInformation= p.get< bool >("DoWriteTracksInfo",false);
    DoWriteVerticesInformation = p.get< bool >("DoWriteVerticesInfo",true);
    DoAddTracksEdep         = p.get< bool >("DoAddTracksEdep",false);
    DoWriteGENIEInformation = p.get< bool >("DoWriteGENIEInfo",true);
    DoWriteEventsInformation= p.get< bool >("DoWriteEventsInfo",true);
    DoOnlySwT               = p.get< bool >("DoOnlySwT",false);
    fSwTModuleLabel         = p.get< std::string >("SwTModuleLabel");
    fSwTAlgoModuleLabel     = p.get< std::string >("SwTAlgoModuleLabel");    
    TruncMeanRad            = p.get< float >("TruncMeanRadLabel",1.0);
    
    fPIDModuleLabel         = p.get< std::string >("PIDModuleLabel");
    fCaliPIDModuleLabel     = p.get< std::string >("CaliPIDModuleLabel");
    
    // For Marco' Flash-matching
    FlashRangeStart         = p.get<double>     ("FlashVetoTimeStart",    3);
    FlashRangeEnd           = p.get<double>     ("FlashVetoTimeEnd",      5);
    DoOpDetSwap             = p.get<bool>       ("DoOpDetSwap", false);
    OpDetSwapMap            = p.get<std::vector<int> >("OpDetSwapMap");
    FlashMatch_mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));
    
    // event weight
    fgenieEventWeightModuleLabel    = p.get< std::string >("EventWeightModuleLabel");
    fluxEventWeightModuleLabel      = p.get< std::string >("fluxEventWeightModuleLabel","");
    DoDebugRSE              = p.get<bool>                           ("DoDebugRSE",false);
    DebugRSEarray           = p.get<std::vector<std::vector<int>> > ("DebugRSE",{});
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAna::CollectEventWeights (art::Event const & evt){
    int fDebug = 3;
    
    // superceeded by the event-weights used for flux variation analysis
    // June-14, 2018
    // reweighting for different mA values,
    // to extract mA from CCQE-like events
    Debug ( fDebug , "in ub::ErezCCQEAna::CollectEventWeights()" );
    
    art::Handle< std::vector <evwgh::MCEventWeight> > MCEventWeightlistHandle;
    std::vector< art::Ptr <evwgh::MCEventWeight> > MCEventWeightlist;
    //    if (evt.getByLabel( fgenieEventWeightModuleLabel , MCEventWeightlistHandle))
    if (evt.getByLabel( fluxEventWeightModuleLabel , MCEventWeightlistHandle))
        art::fill_ptr_vector( MCEventWeightlist , MCEventWeightlistHandle );

    
    Debug ( fDebug , "MCEventWeightlist size is %" , MCEventWeightlist.size());
    if (MCEventWeightlist.size() > 0) {
        Debug ( fDebug , "stepping through the first element ...");
        art::Ptr<evwgh::MCEventWeight> evt_wgt = MCEventWeightlist.at(0); // Just for the first nu interaction
        std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
        for(auto it : evtwgt_map) {
            Debug(fDebug,"for(auto it : evtwgt_map), it.first: %, it.second.size(): %",it.first,it.second.size());
            event_weight_names.push_back(it.first);
            event_weight_values.push_back(it.second.at(0));
        }
    }
    if ( debug > fDebug ) {
        for (size_t i=0; i<event_weight_names.size(); i++) {
            SHOW2(event_weight_names.at(i),event_weight_values.at(i));
        }
    }
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezCCQEAna)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -












