////////////////////////////////////////////////////////////////////////
// Name: DAQHeaderTimeUBooNE.h
//
// Purpose: Class to hold extended DAQ header information, in particular
//          GPS and NTP (host) event time stamp.
//
// Created: 10-Aug-2017
//
////////////////////////////////////////////////////////////////////////

#ifndef DAQHEADERTIMEUBOONE_H
#define DAQHEADERTIMEUBOONE_H

#include <time.h>

namespace raw {

  class DAQHeaderTimeUBooNE {
  public:

    // Constructors.

    DAQHeaderTimeUBooNE();
    DAQHeaderTimeUBooNE(time_t gps_time, time_t ntp_time, time_t gps_adj_time);

    //Set Methods
    void SetGPSTime(time_t t);
    void SetNTPTime(time_t t);
    void SetGPSAdjTime(time_t t);

    // Accessors.

    time_t gps_time() const {return fGPSTime;}
    time_t ntp_time() const {return fNTPTime;}
    time_t gps_adj_time() const {return fGPSAdjTime;}

  private:

    // Data members.

    time_t fGPSTime;    // (high, low)=(seconds, nanoseconds)
    time_t fNTPTime;    // (high, low)=(seconds, nanoseconds)
    time_t fGPSAdjTime; // (high, low)=(seconds, nanoseconds)
  };
}

inline void raw::DAQHeaderTimeUBooNE::SetGPSTime(time_t t)    { fGPSTime = t; }
inline void raw::DAQHeaderTimeUBooNE::SetNTPTime(time_t t)    { fNTPTime = t; }
inline void raw::DAQHeaderTimeUBooNE::SetGPSAdjTime(time_t t) { fGPSAdjTime = t; }

#endif // DAQHEADERTIMEUBOONE_H
