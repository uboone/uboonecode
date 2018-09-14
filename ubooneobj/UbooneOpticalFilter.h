/** ****************************************************************************
 * @file UbooneOpticalFilter.h
 * @brief Definition of basic information for the common optical filter
 * @author wketchum@fnal.gov
 * 
 * ****************************************************************************/

#ifndef UBOONEOBJ_UBOONEOPTICALFILTER_H
#define UBOONEOBJ_UBOONEOPTICALFILTER_H

#include <stdint.h>

namespace uboone {

  class UbooneOpticalFilter {

  public:
  UbooneOpticalFilter():
    fPE_Beam(-999),fPE_Veto(-999),fPMT_MaxFraction(-999) {} 
    
#ifndef __GCCXML__
  public:
      
    UbooneOpticalFilter(float pe_b, float pe_veto, float pmt_maxfrac)
      { fPE_Beam = pe_b; fPE_Veto = pe_veto; fPMT_MaxFraction = pmt_maxfrac; }

    float PE_Beam() const;         /// integrated pe in "beam" window
    float PE_Veto() const;         /// integrated pe in "veto" region
    float PMT_MaxFraction() const; /// max fraction of pe in single PMT

#endif // !__GCCXML__
  private:

    float fPE_Beam;
    float fPE_Veto;
    float fPMT_MaxFraction;

  }; // class UbooneOpticalFilter()
  
#ifndef __GCCXML__
  inline float uboone::UbooneOpticalFilter::PE_Beam() const { return fPE_Beam; }
  inline float uboone::UbooneOpticalFilter::PE_Veto() const { return fPE_Veto; }
  inline float uboone::UbooneOpticalFilter::PMT_MaxFraction() const { return fPMT_MaxFraction; }
#endif // !__GCCXML__
  
} // namespace uboone


#endif // UBOONEOBJ_UBOONEOPTICALFILTER_H

////////////////////////////////////////////////////////////////////////
