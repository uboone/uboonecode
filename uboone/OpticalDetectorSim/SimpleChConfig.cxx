#ifndef SIMPLECHCONFIG_CXX
#define SIMPLECHCONFIG_CXX
#include "SimpleChConfig.h"

namespace opdet {

  const std::map<unsigned int, float>& SimpleChConfig::GetFloat (const ChConfigType_t type)
  {
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = fFloatParams.find( type );
    if ( par_chdata!=fFloatParams.end() )
      return (*par_chdata).second;
    
    // should not compile if mistake made
    throw UBOpticalException("Invalid parameter type!");
  }

  const std::map<unsigned int, int>& SimpleChConfig::GetInt (const ChConfigType_t type)
  {
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = fIntParams.find( type );
    if ( par_chdata!=fIntParams.end() )
      return (*par_chdata).second;
    
    // should not compile if mistake made
    throw UBOpticalException("Invalid parameter type!");
  }

  float SimpleChConfig::GetFloat( const ChConfigType_t type, const unsigned int ch)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto const& par_chdata = GetFloat( type );
    auto par_value = par_chdata.find( ch );
    if ( par_value!=par_chdata.end() )
      return (*par_value).second;
    
    char err[100];
    sprintf(err, "Invalid channel number=%d provided!",ch );
    throw UBOpticalException(err);
  }

  int SimpleChConfig::GetInt( const ChConfigType_t type, const unsigned int ch)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto const& par_chdata = GetInt( type );
    auto par_value = par_chdata.find( ch );
    if ( par_value!=par_chdata.end() )
      return (*par_value).second;
    
    char err[100];
    sprintf(err, "Invalid channel number=%d provided!",ch );
    throw UBOpticalException(err);
  }

  float SimpleChConfig::GetFloatAtSpectrumValue ( const ChSpectrumConfigType_t type, const unsigned int ch, const float energy)
  {
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }

    // Find the parameter
    auto iter = fSpectrumFloatParams.find(type);
    if (iter == fSpectrumFloatParams.end()) {
      char err[100];
      sprintf(err, "Cannot find parameter %d during the call of GetFloatAtSpectrumValue",type );
      throw UBOpticalException(err);
    }
    std::map< unsigned int, std::vector< float > > ch_to_vector_map = iter->second;

    // Find the channel
    auto iter2 = ch_to_vector_map.find(ch);
    if (iter2 == ch_to_vector_map.end()) {
      char err[100];
      sprintf(err, "Invalid channel number=%d provided!",ch );
      throw UBOpticalException(err);
    }
    std::vector< float > values_per_energy = iter2->second;

    // Go from energy to wavelength
    double wl = 1.2398e3/energy;

    // Now interpolate
    return this->Interpolate(fSpectrum, values_per_energy, wl);

  }

  float SimpleChConfig::Interpolate( std::vector<float> & xData, std::vector<float> & yData, float x, bool extrapolate )
  {
    int size = xData.size();
  
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                                                 // special case: beyond right end
    {
       i = size - 2;
    }
    else
    {
       while ( x > xData[i+1] ) i++;
    }
    float xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
       if ( x < xL ) yR = yL;
       if ( x > xR ) yL = yR;
    }
  
    float dydx = ( yR - yL ) / ( xR - xL );                                    // gradient
 
    //std::cout << "SimpleChConfig::Interpolate :: Wavelength " << x << ", xL " << xL << ", xR " << xR << ", yL " << yL << ", yR " << yR << ", dydx " << dydx << "  =>  " << yL + dydx * ( x - xL ) << std::endl; 

    return yL + dydx * ( x - xL );                                              // linear interpolation
  }
}

#endif

