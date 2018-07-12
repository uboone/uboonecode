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
    fNTPTime(0),
    fPPSsec(0),
    fPPSmicro(0),
    fPPSnano(0),
    fTrigFrame(0),
    fTrigSample(0),
    fTrigDiv(0),
    fTrigPPSFrame(0),
    fTrigPPSSample(0),
    fTrigPPSDiv(0)
  {}
}
