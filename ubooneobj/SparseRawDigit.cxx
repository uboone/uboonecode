//==============================================================================
//
// Name: SparseRawDigit.cxx
//
// Purpose: Implementation for class SparseRawDigit.
//
// Created: 17-Nov-2017  H. Greenlee
//
//==============================================================================

#include "ubooneobj/SparseRawDigit.h"

namespace raw {

  // Default constructor.

  SparseRawDigit::SparseRawDigit() :
    fChannel(InvalidChannelID),
    fView(geo::kUnknown),
    fPedestal(0.),
    fSigma(0.),
    fADC()
  {}

  // Initializing constructor.

  SparseRawDigit::SparseRawDigit(ChannelID_t channel,
				 geo::View_t view,
				 float pedestal,
				 float sigma,
				 const lar::sparse_vector<short>& adcvec) :
    fChannel(channel),
    fView(view),
    fPedestal(pedestal),
    fSigma(sigma),
    fADC(adcvec)
  {}

  // Initializing move constructor.

  SparseRawDigit::SparseRawDigit(ChannelID_t channel,
				 geo::View_t view,
				 float pedestal,
				 float sigma,
				 lar::sparse_vector<short>&& adcvec) :
    fChannel(channel),
    fView(view),
    fPedestal(pedestal),
    fSigma(sigma),
    fADC(std::move(adcvec))
  {}

  // Waveform as vector.

  std::vector<short> SparseRawDigit::SparseRawDigit::ADCvector() const
  {
    return std::vector<short>(fADC.begin(), fADC.end());
  }

  // Set Pedestal.

  void SparseRawDigit::SetPedestal(float ped, float sigma)
  {
    fPedestal = ped;
    fSigma = sigma;
  }

} // namespace raw.
