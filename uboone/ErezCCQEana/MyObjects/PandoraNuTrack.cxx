#ifndef PANDORANUTRACK_CXX
#define PANDORANUTRACK_CXX

#include "PandoraNuTrack.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PandoraNuTrack::PandoraNuTrack( Int_t frun, Int_t fsubrun, Int_t fevent
                               ,Int_t ftrack_id
                               ,Float_t flength
                               ,TVector3 fstart_pos, TVector3 fend_pos):
run(frun),
subrun(fsubrun),
event(fevent),
track_id(ftrack_id),
length(flength),
start_pos(fstart_pos),
end_pos(fend_pos)
{
    SetRecDirection();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::Print( bool DoPrintPandoraNuFeatures ) const{
    
    cout << "\033[31m" << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl
    << "track " << track_id << endl << "-------------------"    << "\033[0m" << endl;
    
    SHOWTVector3(start_pos);
    SHOWTVector3(end_pos);
    SHOWTVector3(rec_dir);

    if (DoPrintPandoraNuFeatures){
        
        
        SHOW( PIDa ); // SHOW3( PIDaPerPlane[0] , PIDaPerPlane[1] , PIDaPerPlane[2] );
        PrintPhys(length,"cm");
        Printf("rec theta=%.1f, phi=%.1f deg.", r2d*rec_dir.Theta(), r2d*rec_dir.Phi() );
        SHOW(IsFlipped);
        
    }
    if (MCpdgCode!=-9999){
        cout << "........................" << endl << "MC information " << endl ;
        SHOW(MCpdgCode);
        SHOW(mcevent_id);
        SHOWTVector3(truth_start_pos);
        SHOWTVector3(truth_end_pos);
        SHOWTVector3(truth_dir);
        SHOWTLorentzVector(truth_momentum);
        SHOW(truth_mother);
        SHOW(truth_process);
        SHOW(truth_origin);
        Printf("truth theta=%.1f, phi=%.1f deg.", r2d*truth_dir.Theta(), r2d*truth_dir.Phi() );
        if ( r2d*fabs(truth_dir.Theta()-rec_dir.Theta()) > 90 ) {
            Printf("theta - truth-theta = %.1f(!)",r2d*fabs(truth_dir.Theta()-rec_dir.Theta()));
        }
        cout << "........................" << endl;
    }
    cout << "closest-flash:" << endl;
    ClosestFlash.Print();
    cout << "\033[31m" << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << "\033[0m" << endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::SetStartEndPlane(Int_t plane ,
                                      Int_t start_wire, Int_t start_time ,
                                      Int_t end_wire, Int_t end_time ){
    
    roi[plane] = box( start_wire , start_time , end_wire , end_time );
    switch (plane) {
        case 0:
            start_wire_u = start_wire;
            start_time_u = start_time;
            end_wire_u = end_wire;
            end_time_u = end_time;
            break;
        case 1:
            start_wire_v = start_wire;
            start_time_v = start_time;
            end_wire_v = end_wire;
            end_time_v = end_time;
            break;
        case 2:
            start_wire_y = start_wire;
            start_time_y = start_time;
            end_wire_y = end_wire;
            end_time_y = end_time;
            break;
            
        default:
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartWire (int plane) const{
    switch (plane) {
        case 0:
            return start_wire_u;
            break;
        case 1:
            return start_wire_v;
            break;
        case 2:
        default:
            return start_wire_y;
            break;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartTime (int plane) const{
    switch (plane) {
        case 0:
            return start_time_u;
            break;
        case 1:
            return start_time_v;
            break;
        case 2:
        default:
            return start_time_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndWire (int plane) const{
    switch (plane) {
        case 0:
            return end_wire_u;
            break;
        case 1:
            return end_wire_v;
            break;
        case 2:
        default:
            return end_wire_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndTime (int plane) const{
    switch (plane) {
        case 0:
            return end_time_u;
            break;
        case 1:
            return end_time_v;
            break;
        case 2:
        default:
            return end_time_y;
            break;
    }
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::DistanceFromPoint( TVector3 position, std::string * fStartOrEnd  ){
    Float_t DistanceStart, DistanceEnd , distance = 1000;
    std::string StartOrEnd = "None";
    
    DistanceStart = ( start_pos - position).Mag();
    DistanceEnd = ( end_pos - position).Mag();
    if ( DistanceStart < DistanceEnd ){
        StartOrEnd = "Start";
        distance = DistanceStart;
    }
    else{
        StartOrEnd = "End";
        distance = DistanceEnd;
    }
    
    return distance;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::ClosestDistanceToOtherTrack( PandoraNuTrack other_track, std::string * fStartOrEnd ){
    Float_t MinDistanceToOtherTrack = 10000;
    std::string StartOrEnd = "None";
    Float_t DistanceStartStart = (start_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartStart)     {MinDistanceToOtherTrack = DistanceStartStart; StartOrEnd = "Start";}
    
    Float_t DistanceStartEnd = (start_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartEnd)       {MinDistanceToOtherTrack = DistanceStartEnd; StartOrEnd = "Start";}
    
    Float_t DistanceEndStart = (end_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndStart)       {MinDistanceToOtherTrack = DistanceEndStart; StartOrEnd = "End";}
    
    Float_t DistanceEndEnd = (end_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndEnd)         {MinDistanceToOtherTrack = DistanceEndEnd; StartOrEnd = "End";}
    
    
    if (fStartOrEnd!=nullptr) *fStartOrEnd = StartOrEnd;
    
    return MinDistanceToOtherTrack;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::FlipTrack(){
    
    // flip start and end positions
    TVector3 tmp_pos = start_pos;
    start_pos   = end_pos;
    end_pos     = tmp_pos;
    
    // change angles
    rec_dir = -rec_dir;
//    theta       = 3.1416 - theta;
//    phi         = (phi > 0 ? phi-3.1416 : phi+3.1416);
    
    IsFlipped  = !IsFlipped;
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::GetDis2Flash (flash Flash) const {
    float TrackZcenter = (start_pos.z() + end_pos.z())/2.;
    float FlashZcenter = Flash.GetZcenter();
    float Zdis = TrackZcenter - FlashZcenter;
    float TrackYcenter = (start_pos.y() + end_pos.y())/2.;
    float FlashYcenter = Flash.GetYcenter();
    float Ydis = TrackYcenter - FlashYcenter;
    float YZdistance = sqrt( Zdis*Zdis + Ydis*Ydis );
    return YZdistance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::GetDis2ClosestFlash () const {
    return GetDis2Flash(ClosestFlash);
}

#endif

























