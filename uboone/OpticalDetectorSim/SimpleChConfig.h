/**
 * \file SimpleChConfig.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class SimpleChConfig
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetector

    @{*/
#ifndef SIMPLECHCONFIG_H
#define SIMPLECHCONFIG_H

#include <map>
#include <vector>
#include <iostream>

#include "UBOpticalException.h"
#include "UBOpticalConstants.h"

namespace opdet {
  /**
     \class SimpleChConfig
     User defined class SimpleChConfig ... these comments are used to generate
     doxygen documentation!
  */
  class SimpleChConfig{
    
  protected:
    
    /// Default constructor
  SimpleChConfig() : fInitialized(false) {};

    /// Default destructor
    virtual ~SimpleChConfig(){};

  public:

    const std::map<unsigned int, float>& GetFloat (const ChConfigType_t type);
    const std::map<unsigned int, int>&   GetInt   (const ChConfigType_t type);

    float GetFloat ( const ChConfigType_t type, const unsigned int ch);
    int   GetInt   ( const ChConfigType_t type, const unsigned int ch);

    float GetFloatAtSpectrumValue ( const ChSpectrumConfigType_t type, const unsigned int ch, const float energy);

  protected:
    bool fInitialized;
    virtual void doInitialization() = 0;
    float Interpolate( std::vector<float> & xData, std::vector<float> & yData, float x, bool extrapolate = true);

    // These are parameters that depend on channel number only
    std::map< ChConfigType_t, std::map< unsigned int, float > > fFloatParams; // [ Parameter Enum, [ Channel Number, Value ] ]
    std::map< ChConfigType_t, std::map< unsigned int, int   > > fIntParams;   // [ Parameter Enum, [ Channel Number, Value ] ]

    // These are parameters that depend on both the channel number and the energy of the photon
    std::vector< float > fSpectrum;                                                                            // The energy spectrum discrete values
    std::map< ChSpectrumConfigType_t, std::map< unsigned int, std::vector< float > > > fSpectrumFloatParams;   // [ Spectrum Parameter Enum, [ Channel Number, Vector of Values ] ]

  };

}

#endif
/** @} */ // end of doxygen group 

