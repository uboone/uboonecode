/**
 * @file   PMTRemapProvider.cxx
 * @brief  MicroBooNE service provider that provides a map to correct PMT OpChannel numbers
 * @author Brandon Eberly (eberly@slac.stanford.edu)
 */

#include "uboone/Utilities/PMTRemapProvider.h"

namespace util {

  PMTRemapProvider::PMTRemapProvider(fhicl::ParameterSet const& pset) {
  
    //hard-code the map
    fOriginal_to_corrected_map[26] = 27;
    fOriginal_to_corrected_map[27] = 28;
    fOriginal_to_corrected_map[28] = 29;
    fOriginal_to_corrected_map[29] = 30;
    fOriginal_to_corrected_map[30] = 31;
    fOriginal_to_corrected_map[31] = 26;
  }
  
  unsigned int PMTRemapProvider::CorrectedOpChannel(unsigned int orig) const {
    
    auto it = fOriginal_to_corrected_map.find(orig%100);
    unsigned int hundreds = orig/100;
    
    if ( it == fOriginal_to_corrected_map.end() ) {
      return orig;
    }
    else return hundreds*100 + it->second;
  }
  
  unsigned int PMTRemapProvider::OriginalOpChannel(unsigned int corr) const {
    
    unsigned int mod = corr%100;
    unsigned int hundreds = corr/100;
    
    for (auto it = fOriginal_to_corrected_map.begin(); it != fOriginal_to_corrected_map.end(); ++it) {
      if (it->second == mod) {
        return hundreds*100 + it->first;
      }
    }
    
    return corr;
  }
} //end namespace util
      
      
