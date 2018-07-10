////////////////////////////////////////////////////////////////////////
// $Id: DAQHeaderTimeUBooNE.cxx,v 1.0 2017/08/10 19:34:20 brebel Exp $
//
// DAQHeaderTimeUBooNE class
//
// kirby@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RawData/DAQHeader.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

namespace raw{

  //----------------------------------------------------------------------
  // Default constructor.

  DAQHeaderTimeUBooNE::DAQHeaderTimeUBooNE() :
    fGPSTime(0),
    fGPSAdjTime(0),
    fNTPTime(0),
    fPPSsec(0),
    fPPSmicro(0),
    fPPSnano(0)
  {}

  //----------------------------------------------------------------------
  // Initializing constructor.
  DAQHeaderTimeUBooNE::DAQHeaderTimeUBooNE(time_t gps_time, time_t gps_adj_time, time_t ntp_time,
					   uint32_t pps_sec, uint32_t pps_micro,
					   uint32_t pps_nano) :
    fGPSTime(gps_time),
    fGPSAdjTime(gps_adj_time),
    fNTPTime(ntp_time),
    fPPSsec(pps_sec),
    fPPSmicro(pps_micro),
    fPPSnano(pps_nano)
  {}
}
