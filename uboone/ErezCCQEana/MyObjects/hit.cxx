#ifndef hit_CXX
#define hit_CXX

#include "hit.h"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
hit::hit(Int_t fplane, Int_t fwire, Int_t fid, float fpeakT, float fcharge):
hit_plane(fplane),
hit_wire(fwire),
hit_id(fid),
hit_peakT(fpeakT),
hit_charge(fcharge)
{}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool hit::InBox(box fbox) const {
    // check if hit is inside a box
    if (
        ( fbox.GetStartWire() <= hit_wire && hit_wire <= fbox.GetEndWire() )
        &&
        ( fbox.GetStartTime() <= hit_peakT && hit_peakT <= fbox.GetEndTime() )
        )
    {
        return true;
    }
    return false;
}



//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//bool hit::InSphere(float cwire, float ctime, float r , float f_timeticks2cm, float f_wire2cm ) const {
//    // check if hit is inside a sphere of radius r [cm]
//    // localized around (cwire,ctime)
//    d2_hit_center = (hit_wire - cwire)*(hit_wire - cwire)*f_wire2cm + (hit_time - ctime)*(hit_time - ctime)*f_timeticks2cm;
//    r2 = r*r;
//    // no nead to use the expensive sqrt method here, we can use the squared distance for this check
//    if ( d2_hit_center < r2 )        return true;
//    return false;
//}
//




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void hit::Print() const{
    std::cout << "hit " << hit_id << ", plane " << hit_plane << "/ wire " << hit_wire << "/ peakT " << hit_peakT
    << ", " << "\t trkKey: " << hit_trkKey << std::endl;
}




#endif
