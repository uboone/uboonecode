/**
 * \class CRTData
 *
 * \ingroup crt
 *
 * \brief The basic structure of data coming from the CRT in both
 * Simulations and Data.
 *
 * \author $Author: Kevin Wierman<kevin.wierman@pnnl.gov> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2016/12/12 $
 *
 * Contact: kevin.wierman@pnnl.gov
 *
 * Created on: Tuesday, December 13, 2016
 *
 */


#ifndef CRTSimData_hh_
#define CRTSimData_hh_

#include <cstdint>

namespace crt {

  class CRTSimData {

    /// The channel number which can be referenced to module, strip via the channel map
    uint32_t fChannel;
    /// Time relative to GPS
    uint32_t fT0;
    /// Time relative to beam trigger
    uint32_t fT1;
    /// ADC value returned by CRT
    uint32_t fADC;
   public:
    /// Default constructor
    CRTSimData();
    /// To be used when constructing from existing data
    CRTSimData(uint32_t channel, uint32_t t0, uint32_t t1, uint32_t adc);
    virtual ~CRTSimData();
    /// Channel getter
    uint32_t Channel();
    /// T0 Getter
    uint32_t T0();
    /// T1 Getter
    uint32_t T1();
    /// ADC Getter
    uint32_t ADC();
  };

}

#endif
