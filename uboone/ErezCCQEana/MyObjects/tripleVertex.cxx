#ifndef tripleVertex_CXX
#define tripleVertex_CXX

#include "tripleVertex.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
tripleVertex::tripleVertex(Int_t frun,Int_t fsubrun,Int_t fevent,Int_t fID):
run(frun),
subrun(fsubrun),
event(fevent),
vertex_id(fID)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void tripleVertex::AddTrack (PandoraNuTrack ftrack){
    
    tracks.push_back(ftrack);
    Ntracks=(int)tracks.size();
    AddTrackID( ftrack.GetTrackID() );
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool tripleVertex::IncludesTrack (Int_t ftrack_id){
    for (auto track:tracks){
        if (track.GetTrackID() == ftrack_id) return true;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<size_t> tripleVertex::sort_l(const vector<PandoraNuTrack> &v) {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1].GetLength() > v[i2].GetLength();});
    return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool tripleVertex::SortTracksByLength(){
    std::vector<PandoraNuTrack> tmp_tracks = tracks;
    for (auto i:sort_l( tmp_tracks )){
        tracks_lengthsorted.push_back( tmp_tracks.at(i) );
    }
    ShortestTrack = tracks_lengthsorted.at( tmp_tracks.size() - 1 );
    LongestTrack = tracks_lengthsorted.at( 0 );
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void tripleVertex::SetTracksRelations(){
    // March 2018
    // the angles (delta_phi and delta_theta) are calculated in degrees (!)
    for(int i = 0; i < Ntracks ; i++){
        
        tracks_dis_from_vertex.push_back( tracks[i].DistanceFromPoint( position ) );
        
        std::vector<float> distances_track_i;
        std::vector<float> delta_phi_track_i;
        std::vector<float> delta_theta_track_i;
        
        for(int j = 0; j < Ntracks ; j++){
            distances_track_i.push_back( tracks[i].ClosestDistanceToOtherTrack(tracks[j]) );
            delta_phi_track_i.push_back( r2d*fabs(tracks[i].GetPhi() - tracks[j].GetPhi()) );
            delta_theta_track_i.push_back( r2d*fabs(tracks[i].GetTheta() - tracks[j].GetTheta() ) );
            
            
            if (i==0 && j!=i) {
                distances_ij.push_back( distances_track_i.back() );
                delta_phi_ij.push_back( delta_phi_track_i.back() );
                delta_theta_ij.push_back( delta_theta_track_i.back() );
            }
        }
        tracks_distances.push_back( distances_track_i );
        tracks_delta_phi.push_back( delta_phi_track_i );
        tracks_delta_theta.push_back( delta_theta_track_i );
        
        
        distances_track_i.clear();
        delta_phi_track_i.clear();
        delta_theta_track_i.clear();
    }
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<PandoraNuTrack> tripleVertex::RemoveTrackFromVector( vector<PandoraNuTrack> TracksVector , PandoraNuTrack TrackToBeRemoved ){
    std::vector<PandoraNuTrack> tmp_tracks = TracksVector;
    TracksVector.clear();
    for (auto t : tmp_tracks){
        if ( t != TrackToBeRemoved ){
            TracksVector.push_back(t);
        }
    }
    return TracksVector;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<PandoraNuTrack> tripleVertex::RemoveTracksFromVector( vector<PandoraNuTrack> TracksVector , vector<PandoraNuTrack> TracksToBeRemoved ){
    // loop over tracks to be removed
    // and remove each of them
    for (auto t:TracksToBeRemoved) TracksVector = RemoveTrackFromVector( TracksVector , t );
    return TracksVector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void tripleVertex::Print(bool DoPrintTracks,bool DoPrintFull) const {
    
    cout << "\033[35m" << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n vertex " << vertex_id ;
    if (DoPrintFull==false){
        cout << ", tracks: ";
        for (auto t: tracks) {
            cout
            << t.GetTrackID() << "(pdg "<< t.GetMCpdgCode()<<", MCParticle " << t.GetTruthMCParticle()
            <<"),  ";
        }
        cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< "\033[0m" << endl;
        return;
    }
    cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< "\033[0m" << endl;
    
    SHOWTVector3( position );
    
    for (auto t: tracks) {
        Printf("track %d (pdg %d, MCParticle %d), distance %.1f cm",t.GetTrackID(),t.GetMCpdgCode(),t.GetTruthMCParticle(),t.DistanceFromPoint( position ));
    }
    
    // tracks
    if (!tracks.empty()){
        
        if (DoPrintTracks) {
            cout << "\033[33m" << tracks.size() << " tracks in vertex " << "\033[31m" << endl;
            for (auto t: tracks) {
                if (t.GetTrackID()!=-100)   t.Print( true );
                else                    Printf("unidentified/non-reconstructed track. continuing...");
            }
        }
        
        // inter-tracks distances
        cout << "\033[33m" << tracks.size()*tracks.size() << " inter-tracks distances " << "\033[31m" << endl;
        for(auto vec : tracks_distances) {
            for(auto x : vec) std::cout << setprecision(2) << x << "\t";
            std::cout << std::endl;
        }
    }
    else {
        cout << "\033[33m" << "no reconstructed tracks in vertex " << "\033[31m" << endl;
    }
    
    PrintLine();
    cout << "truth topology: " << TruthTopologyString << endl;
}


#endif
