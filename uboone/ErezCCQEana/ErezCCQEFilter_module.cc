
////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEFilter
// Plugin Type: filter (art v2_05_01)
// File:        ErezCCQEFilter_module.cc
//
// Generated at Tue Jun 19 10:46:33 2018 by Erez Cohen using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// my objects...
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/flash.h"
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"
#include "uboone/ErezCCQEana/MyObjects/pairVertex.h"
#include "uboone/ErezCCQEana/MyObjects/TruncMean.h"

#include <memory>


constexpr int   kMaxTrack = 1000;  //maximum number of tracks
constexpr float kMaxInterTrackDistance = 11; // maximal distance for clustering


namespace ub {
    class ErezCCQEFilter;
}


class ub::ErezCCQEFilter : public art::EDFilter {
public:
    explicit ErezCCQEFilter(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEFilter(ErezCCQEFilter const &) = delete;
    ErezCCQEFilter(ErezCCQEFilter &&) = delete;
    ErezCCQEFilter & operator = (ErezCCQEFilter const &) = delete;
    ErezCCQEFilter & operator = (ErezCCQEFilter &&) = delete;

    
    
    
    // Required functions.
    bool      filter (art::Event & e) override;
    
    // Optional functions.
    void    beginJob () override;
    void      endJob () override;
    
    
    
    
    // My methods
    void               ResetVars ();
    void           CollectTracks (art::Event & evt);
    void         AnalyzeVertices ();
    void ClusterTracksToVertices ();
    void  FilterGoodPairVertices ();
    void       ConstructVertices ();
    void             reconfigure (fhicl::ParameterSet const & p);
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=0;
    void Debug(Int_t verobosity_level, const char* format)     {
        if ( debug < verobosity_level ) return;
        std::cout << format << std::endl;
    }
    template<typename T, typename... Targs>
    void Debug(Int_t verobosity_level, const char* format, T value, Targs... Fargs) {
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
    
    Int_t   run=0, subrun=0, event=0;

    Int_t   ctr_events=0, ctr_events_in=0, ctr_events_out=0;
    Int_t   Ntracks=0;

    
    
    std::string fTrackModuleLabel;
    
    std::vector<PandoraNuTrack>     tracks;
    std::vector<pairVertex>         vertices;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::ErezCCQEFilter::ErezCCQEFilter(fhicl::ParameterSet const & p){
    // Call appropriate produces<>() functions here.
    reconfigure(p);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::ResetVars(){
    
    run = subrun = event = -9999;
    Ntracks = 0 ;
    tracks.clear();
    vertices.clear();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::ErezCCQEFilter::filter(art::Event & evt)
{
    // In this routine we will filter out events which have no candidate pair-vertices.
    // In MC-overlay this happens in ~ 0.8 of the times
    // and in data about 0.88 of the times.
    // We do this by collecting pair vertices in our standard manner,
    // and then asking if the size of the vertices vector is greater than zero
    ctr_events++;
    
    ResetVars();
    CollectTracks(evt);
    ConstructVertices();
    
    
    if (vertices.size()>0) {
        Debug(0,"keeping run %/sub %/event % , since the size of vertices vector is %",run,subrun,event,vertices.size());
        ctr_events_in ++;
        return true;
    }

    ctr_events_out ++;
    Debug(0, "filtering out run %/sub %/event %, since there are no CCQE candidate vertices ...",run,subrun,event);
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::beginJob(){
    // Implementation of optional member function here.
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::endJob(){
    // Implementation of optional member function here.
    Debug(0,"ending job of ub::ErezCCQEFilter");
    Debug(0,"stepped through % events",ctr_events);
    Debug(0,"number of events kept: % (%%)",ctr_events_in,100.*(float)ctr_events_in/ctr_events);
    Debug(0,"number of events filtered out: % (%%)",ctr_events_out,100.*(float)ctr_events_out/ctr_events);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::CollectTracks(art::Event & evt){
    // ----------------------------------------
    // tracks information
    // ----------------------------------------
    Debug(2,"// * ub::ErezCCQEFilter::CollectTracks()");
    run = evt.run();  subrun = evt.subRun();    event = evt.id().event();
    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
        art::fill_ptr_vector(tracklist, trackListHandle);

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
        tracks.push_back( track );
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::ConstructVertices(){
    // cluster all tracks at close proximity to vertices
    ClusterTracksToVertices();
    // analyze these vertices: inter-tracks distances, angles...
    AnalyzeVertices();
    // retain only vertices with pairs of 2-tracks at close proximity
    FilterGoodPairVertices();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::ClusterTracksToVertices(){
    
    // June-20, 2018
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
        if (FoundCloseTracks) { vertices.push_back( vertex ); }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::AnalyzeVertices(){
    if (vertices.size()>0){
        for (auto & v:vertices){
            // after fixing the vertext position, remove far tracks
            Debug( 4 , "v.RemoveFarTracks( kMaxInterTrackDistance );");
            if (v.GetTracks().size()<2) continue;
            v.RemoveFarTracks( kMaxInterTrackDistance );
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::FilterGoodPairVertices(){
    
    // this funciton also kiils all vertices which are not up to standard.
    // to this end, we create a temporary vertices array and clear all vertices
    std::vector<pairVertex> tmp_vertices = vertices;
    vertices.clear();
    
    Debug(3 , "ub::FilterGoodPairVertices::FilterGoodPairVertices()");
    for (auto & v:tmp_vertices) {
        // vertices with only two tracks at close proximity and nothing else
        if (// we are looking for clusters of only two tracks that are fully contained
            v.GetNtracks() == 2
            // if there is a semi-contained track, with start/end point too close to the vertex, we don't want the vertex...
            &&  v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size() == 0
            ){
            vertices.push_back( v );
        }
        // if the vertex is not with only two tracks at close proximity and nothing else, kill it
        else { Debug( 3 , Form("erasing vertex %d, since N(tracks)=%d, N(close semi-contained tracks):%d"
                               , v.GetVertexID(), v.GetNtracks(), (int)(v.CloseSemiContainedTracks( tracks , kMaxInterTrackDistance ).size()) ) );}
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEFilter::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel       = p.get< std::string >  ("TrackModuleLabel","pandoraNu");
}


DEFINE_ART_MODULE(ub::ErezCCQEFilter)
