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

    //Set Methods

    void SetGPSTime(time_t t);
    void SetGPSAdjTime(time_t t);
    void SetNTPTime(time_t t);
    void SetPPSTime(uint32_t sec, uint32_t micro, uint32_t nano);
    void SetTrigTime(uint32_t frame, uint16_t sample, uint16_t div);
    void SetTrigPPSTime(uint32_t frame, uint16_t sample, uint16_t div);

    // Accessors.

    time_t gps_time() const {return fGPSTime;}
    time_t gps_adj_time() const {return fGPSAdjTime;}
    time_t ntp_time() const {return fNTPTime;}

    uint32_t pps_sec() const {return fPPSsec;}
    uint32_t pps_micro() const {return fPPSmicro;}
    uint32_t pps_nano() const {return fPPSnano;}

    uint32_t trig_frame() const {return fTrigFrame;}
    uint16_t trig_sample() const {return fTrigSample;}
    uint16_t trig_div() const {return fTrigDiv;}

    uint32_t trig_pps_frame() const {return fTrigPPSFrame;}
    uint16_t trig_pps_sample() const {return fTrigPPSSample;}
    uint16_t trig_pps_div() const {return fTrigPPSDiv;}

  private:

    // Data members.

    // Complete event times.

    time_t fGPSTime;    // (high, low)=(seconds, nanoseconds)
    time_t fGPSAdjTime; // (high, low)=(seconds, nanoseconds)
    time_t fNTPTime;    // (high, low)=(seconds, nanoseconds)

    // GPS PPS.

    uint32_t fPPSsec;   // PPS seconds.
    uint32_t fPPSmicro; // PPS microseconds.
    uint32_t fPPSnano;  // PPS nanoseconds.

    // Trigger event time.

    uint32_t fTrigFrame;  // Event frame (1.6 ms).
    uint16_t fTrigSample; // Event sample (2 MHz).
    uint16_t fTrigDiv;    // Event division (16 MHz).

    // Trigger PPS time.

    uint32_t fTrigPPSFrame;  // PPS frame (1.6 ms).
    uint16_t fTrigPPSSample; // PPS sample (2 MHz).
    uint16_t fTrigPPSDiv;    // PPS division (16 MHz).
  };
}

inline void raw::DAQHeaderTimeUBooNE::SetGPSTime(time_t t)    { fGPSTime = t; }
inline void raw::DAQHeaderTimeUBooNE::SetGPSAdjTime(time_t t) { fGPSAdjTime = t; }
inline void raw::DAQHeaderTimeUBooNE::SetNTPTime(time_t t)    { fNTPTime = t; }
inline void raw::DAQHeaderTimeUBooNE::SetPPSTime(uint32_t sec, uint32_t micro, uint32_t nano)
{
  fPPSsec=sec;
  fPPSmicro=micro;
  fPPSnano=nano;
}
inline void raw::DAQHeaderTimeUBooNE::SetTrigTime(uint32_t frame, uint16_t sample, uint16_t div)
{
  fTrigFrame=frame;
  fTrigSample=sample;
  fTrigDiv=div;
}
inline void raw::DAQHeaderTimeUBooNE::SetTrigPPSTime(uint32_t frame, uint16_t sample, uint16_t div)
{
  fTrigPPSFrame=frame;
  fTrigPPSSample=sample;
  fTrigPPSDiv=div;
}

#endif // DAQHEADERTIMEUBOONE_H
