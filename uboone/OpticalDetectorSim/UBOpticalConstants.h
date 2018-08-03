#ifndef UBOPTICALCONSTANTS_H
#define UBOPTICALCONSTANTS_H

#include <limits>
#include <climits>

namespace opdet {

  const unsigned short kADC_MAX = 4095;
  const unsigned short kADC_MIN = 0;

  const unsigned short kINVALID_CHANNEL = std::numeric_limits<unsigned short>::max();

  const unsigned short kLogicStartChannel = 40;
  const unsigned short kLogicNChannel     = 8;
  
  enum LogicChannelType_t {
    kFEMChannelBNB  = 46,
    kFEMChannelNuMI = 47
  };

  // These parameters depend on channel number only
  enum ChConfigType_t {
    kPedestalMean=0, // Pedestal mean in ADC count
    kPedestalSpread, // Pedestal standard deviation in ADC count
    kPMTGain,        // p.e/photon
    kQE,             // Quantum efficiency
    kSplitterGain,   // adc/p.e.
    kGainSpread,     // Spread in gain (in fraction)
    kT0,             // T0 in ns
    kT0Spread,       // T0 spread in ns
    kDarkRate,       // Dark Rate in GHz
    kDisc0Threshold, // Discriminator 0 threshold
    kDisc1Threshold, // Discriminator 1 threshold
    kDisc3Threshold, // Discriminator 3 threshold
    
    kChConfigTypeMax
  };

  // These parameters depend on both the channel number and the energy of the photon
  enum ChSpectrumConfigType_t {
    kEnergySpectrum=0, // The energy spectrum (to each value corresponds a QE value below)
    kQESpec,           // Quantum efficiency

    kChSpectrumConfigTypeMax
  };
  
}
#endif


