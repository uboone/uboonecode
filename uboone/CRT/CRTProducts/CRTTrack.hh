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


#ifndef CRTTrack_hh_
#define CRTTrack_hh_

#include <cstdint>
#include <vector>
#include <map>

namespace crt {
  
  struct CRTTrack{
    std::vector<uint8_t> feb_id;
    std::map< uint8_t, std::vector<std::pair<int,double> > > pesmap;
    double peshit;
    uint32_t ts0_s;
    uint16_t ts0_s_err;
    uint32_t ts0_ns;
    uint16_t ts0_ns_err;
    uint32_t ts1_ns; 
    uint16_t ts1_ns_err;                                                                                                                             
    double x1_pos;
    double x1_err;
    double y1_pos;
    double y1_err;
    double z1_pos;
    double z1_err;
    double x2_pos;
    double x2_err;
    double y2_pos;
    double y2_err;
    double z2_pos;
    double z2_err;
    double length;
    double thetaxy;
    double phizy;
    uint32_t ts0_ns_h1;
    uint16_t ts0_ns_err_h1;
    uint32_t ts0_ns_h2;
    uint16_t ts0_ns_err_h2;
       
    CRTTrack() {}
    
  };
  
  
}

#endif
