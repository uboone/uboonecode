/**
 * \file DqdxCorrection.h
 * 
 * \brief Class def header for a class DqdxCorrection
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef IOVDATA_DQDXCORRECTION_H
#define IOVDATA_DQDXCORRECTION_H

#include "larevt/CalibrationDBI/IOVData/ChData.h"
#include "larevt/CalibrationDBI/IOVData/CalibrationExtraInfo.h"

namespace lariov {
  /**
     \class DqdxCorrection   
  *
  */
  class DqdxCorrection : public ChData {
    
    public:
    
      /// Constructor
      DqdxCorrection(unsigned int bin) : 
        ChData(bin),
	fExtraInfo("DqdxCorrection") {}
      
      /// Default destructor
      ~DqdxCorrection() {}
            
      float Correction()    const { return fCorrection; }
      float CorrectionErr() const { return fCorrectionErr; }
      CalibrationExtraInfo const& ExtraInfo() const { return fExtraInfo; }
      
      void SetCorrection(float v)    { fCorrection   = v; }
      void SetCorrectionErr(float v) { fCorrectionErr = v; }
      void SetExtraInfo(CalibrationExtraInfo const& info) { fExtraInfo = info; }
      
    private:
    
      float fCorrection;
      float fCorrectionErr;
      CalibrationExtraInfo fExtraInfo;
      
  }; // end class
} // end namespace lariov

#endif
