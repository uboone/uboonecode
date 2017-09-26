#ifndef bernfebdaq_Utilities_AuxFunctions_hh
#define bernfebdaq_Utilities_AuxFunctions_hh


#include <stdio.h>
#include <map> 
#include <vector> 

#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
#include "artdaq-core/Data/Fragments.hh"
//#include "uboone/CRT/CRTProducts/CRTHit.hh"
//#include "uboone/CRT/CRTProducts/CRTTrack.hh"

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>

namespace crt {
  
  namespace auxfunctions {
    
    int getFEBN(uint64_t febid); //return FEB ID in format (11 - 129)
    
    std::string filePos;
    std::map<int, std::vector<double> > sensor_pos;
    void FillPos(std::string filePos,   std::map <int, std::vector<double> >& sensor_pos); //key = FEB*100+ch //read SiPMs positions from a file a fill a map with access key
    std::vector<double> getPos(int ID, std::map <int, std::vector<double> >& sensor_pos);//return position of a SiPM.
    
    std::string fileFEBDel;
    std::map<int,double> FEBDel;
    void FillFEBDel(std::string fileFEBDel,   std::map <int,double >& FEBDel); //key = FEB //read FEB delay from a file and fill a map
    double getFEBDel(int ID, std::map <int, double >& FEBDel);//return FEB delay
    
    //
    std::string fileGain;
    std::map<int, std::pair<double,double> > sensor_gain;
    void FillGain(std::string fileGain,   std::map <int, std::pair<double,double> >& sensor_gain); //key = FEB*100+ch //read SiPMs gain from a file a fill a map with access key
    std::pair<double,double> getGain(int ID, std::map <int, std::pair<double,double> >& sensor_gain);//return gain and error of a SiPM.
    //
    
    //
    double S1, S2, L;
    std::vector<double> inter_X(double S1, std::vector<double>& pS1, double S2, std::vector<double>& pS2, double L);//return interaction position within the strip given relative intensities.
    double inter_X_error(double S1, double S2, double L);//return error por interaction position.
    
    double getTcorr(std::vector<double>& inpos1, std::vector<double>& inpos2 , double T);//return corrected time, along the fiber. 6.2ns/m 
    
  }
  
}

#endif
