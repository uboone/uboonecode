/**
 * \file PmtCalibrationContainer.h
 *
 * \ingroup IOVData
 * 
 * \brief Class def header for a class PmtCalibrationContainer
 *
 * @author eberly@slac.stanford.edu
 */

/** \addtogroup IOVData

    @{*/
#ifndef IOVDATA_PMTCALIBRATIONCONTAINER_H
#define IOVDATA_PMTCALIBRATIONCONTAINER_H 1

#include "CalibrationDBI/IOVData/ChData.h"

namespace lariov {
  /**
     \class PmtCalibrationContainer
  */
  class PmtCalibrationContainer : public ChData {
    
    public:
    
      /// Constructor
      PmtCalibrationContainer(unsigned int ch) : ChData(ch) {}
      
      /// Default destructor
      ~PmtCalibrationContainer() {}
            
      float                     Amplitude()     const { return fAmplitude; }
      float                     AmplitudeErr()  const { return fAmplitudeErr; }
      float                     Width()         const { return fWidth; }
      float                     WidthErr()      const { return fWidthErr; }
      float                     Area()          const { return fArea; }
      float                     AreaErr()       const { return fAreaErr; }
      const std::vector<float>& AvWaveForm()    const { return fAvWaveform; }
      const std::vector<float>& AvWaveFormErr() const { return fAvWaveformErr; }
      
      void SetAmplitude(float v)                   { fAmplitude     = v; }
      void SetAmplitudeErr(float v)                { fAmplitudeErr  = v; }
      void SetWidth(float v)                       { fWidth         = v; }
      void SetWidthErr(float v)                    { fWidthErr      = v; }
      void SetArea(float v)                        { fArea          = v; }
      void SetAreaErr(float v)                     { fAreaErr       = v; }
      void SetAvWaveform(std::vector<float>& v)    { fAvWaveform    = v; }
      void SetAvWaveformErr(std::vector<float>& v) { fAvWaveformErr = v; }
      
    private:
    
      float fAmplitude;
      float fAmplitudeErr;
      float fWidth;
      float fWidthErr;
      float fArea;
      float fAreaErr;
      std::vector<float> fAvWaveform;
      std::vector<float> fAvWaveformErr;
      
  }; // end class
} // end namespace lariov

#endif
/** @} */ // end of doxygen group 
