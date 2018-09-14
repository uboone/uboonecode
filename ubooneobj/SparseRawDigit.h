//==============================================================================
//
// Name: SparseRawDigit.h
//
// Purpose: Header for data product class SparseRawDigit.
//          This class is similar to the regular RawDigit class, except the data
//          member representing the ADC wavefore is stored as sparse_vector<short>
//          instead of vector<short>.  It is intended that the sparse regions
//          of this class will match those in a corresponding Wire object.
//
// Created: 17-Nov-2017  H. Greenlee
//
//==============================================================================

#ifndef RAW_SPARSE_RAWDIGIT_H
#define RAW_SPARSE_RAWDIGIT_H

#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::Compress_t, raw::Channel_t

namespace raw {

  class SparseRawDigit {
  public:

    // Default constructor.

    SparseRawDigit();

    // Initializing constructor.

    SparseRawDigit(ChannelID_t channel,
		   geo::View_t view,
		   float pedestal,
		   float sigma,
		   const lar::sparse_vector<short>& adcvec);

    // Initializing move constructor.

    SparseRawDigit(ChannelID_t channel,
		   geo::View_t view,
		   float pedestal,
		   float sigma,
		   lar::sparse_vector<short>&& adcvec);

    // Accessors.

    const lar::sparse_vector<short>& ADCs() const;  // Waveform
    std::vector<short> ADCvector() const;           // Waveform as vector.
    size_t NADC() const;                            // Size of waveform.
    short ADC(int i) const;                         // Waveform element access.
    ChannelID_t Channel() const;                    // Readout channel.
    geo::View_t View() const;                       // View.
    float GetPedestal() const;                      // Pedestal.
    float GetSigma() const;                         // Pedestal sigma.

    // Modifiers.

    void SetPedestal(float ped, float sigma = 1.);  // Set pedestal and sigma.

  private:

    // Data members.

    ChannelID_t fChannel;             // Readout channel.
    geo::View_t fView;                // View (from corresponding Wire).
    float fPedestal;                  // Pedestal (from corresonding RawDigit).
    float fSigma;                     // Pedestal sigma (from corresonding RawDigit).
    lar::sparse_vector<short> fADC;   // Waveform.
  };

  // Inlines.

  inline const lar::sparse_vector<short>& SparseRawDigit::ADCs() const {return fADC;}
  inline size_t SparseRawDigit::NADC() const {return fADC.size();}
  inline short SparseRawDigit::ADC(int i) const {return fADC[i];}
  inline ChannelID_t SparseRawDigit::Channel() const {return fChannel;}
  inline geo::View_t SparseRawDigit::View() const {return fView;}
  inline float SparseRawDigit::GetPedestal() const {return fPedestal;}
  inline float SparseRawDigit::GetSigma() const {return fSigma;}

} // namespace raw

#endif // RAW_SPARSE_RAWDIGIT_H
