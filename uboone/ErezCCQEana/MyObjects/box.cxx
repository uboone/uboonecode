#ifndef box_CXX
#define box_CXX

#include "box.h"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
box::box(Int_t fstart_wire, Int_t fstart_time, Int_t fend_wire, Int_t fend_time):
start_wire(fstart_wire),
start_time(fstart_time),
end_wire(fend_wire),
end_time(fend_time)
{

    // order the box start and end positions
    if (start_wire > end_wire){
        
        int tmp_wire= end_wire;
        end_wire   = start_wire;
        start_wire = tmp_wire;
        
    }
    if (start_time > end_time){
        
        int tmp_time = end_time;
        end_time    = start_time;
        start_time  = tmp_time;
        
    }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void box::Print() const{
    
    std::cout
    << "\033[34m"
    << ": (" << start_wire << "," << start_time
    << ") => ("
    << end_wire << "," << end_time << ")"
    << "\033[0m"
    << std::endl;
}


#endif
