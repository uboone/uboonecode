/**
 * \class CRTData
 *
 * \ingroup crt
 *
 * \brief CRT Track Info
 *
 * \author $Author: David Lorca $
 *
 */


#ifndef CRTtzero_hh_
#define CRTtzero_hh_

#include <cstdint>
#include <vector>
#include <map>

namespace crt {
  
  struct CRTTzero{
    uint32_t ts0_s;
    uint16_t ts0_s_err;
    uint32_t ts0_ns;
    uint16_t ts0_ns_err;
    int32_t ts1_ns; 
    uint16_t ts1_ns_err;                    

    int nhits[4];
    
    double pes[4];
    // double xpos[4];
    // double xerr[4];
    // double ypos[4];
    // double yerr[4];
    // double zpos[4];
    // double zerr[4];

       
    CRTTzero() {}
    
  };
  
  
}

#endif
