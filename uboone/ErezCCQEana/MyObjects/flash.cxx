#ifndef flash_CXX
#define flash_CXX

#include "flash.h"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
flash::flash(Float_t ftime, Float_t ftimewidth,
             Float_t fZcenter , Float_t fZwidth,
             Float_t fYcenter , Float_t fYwidth,
             Float_t ftotalPE):
time(ftime),
timewidth(ftimewidth),
Zcenter(fZcenter),
Zwidth(fZwidth),
Ycenter(fYcenter),
Ywidth(fYwidth),
totalPE(ftotalPE)
{

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void flash::Print() const{
    
    std::cout << "\033[34m" << "flash"
    << "time: " << time << ",timewidth: " << timewidth << std::endl
    << "Zcenter: " << Zcenter << ",Zwidth: " << Zwidth << std::endl
    << "Ycenter: " << Ycenter << ",Ywidth: " << Ywidth << std::endl
    << "totalPE: " << totalPE
    << "\033[31m" << std::endl;
}


#endif
