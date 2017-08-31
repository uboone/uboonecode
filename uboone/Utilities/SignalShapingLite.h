////////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingLite.h
///
/// \brief  Generic class for shaping signals on wires.
///
/// \author Brandon Eberly
///
/// This is a generic class for shaping signals on wires during simulation
/// (convolution) and reconstruction (deconvolution).
///
/// This class is based on the SignalShaping class in lardata, but 
/// optimized to use memory as efficiently as possible.  To this end, 
/// only the convolution and deconvolution kernels are stored , and they 
/// are stored as vectors of floats rather than doubles.
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGLITE_H
#define SIGNALSHAPINGLITE_H

#include <vector>
#include "TComplex.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"

namespace util {

  struct ComplexF {
    float Im;
    float Re;
    
    float Rho() const {return TMath::Sqrt(Re*Re+Im*Im);}
  };
    

  class SignalShapingLite {
  public:

    // Constructor, destructor.
    SignalShapingLite();
    virtual ~SignalShapingLite();

    // Accessors.
    const std::vector<float>& Response() const {return fResponse;}
    const std::vector<ComplexF>& ConvKernel() const {return fConvKernel;}
    const std::vector<ComplexF>& DeconvKernel() const {return fDeconvKernel;}

    // Signal shaping methods.

    // Convolute a time series with convolution kernel.
    template <class T> void Convolute(std::vector<T>& func) const;

    // Convolute a time series with deconvolution kernel.
    template <class T> void Deconvolute(std::vector<T>& func) const;


    // Configuration methods.

    // Reset this class to default-constructed state.
    void Reset();

    void set_normflag(bool flag){fNorm = flag;}

    // Add a time domain response function.
    // Updates overall response function and convolution kernel.
    void AddResponseFunction(const std::vector<float>& resp, bool ResetResponse = false );

    // Shift response function in time.
    // Updates overall response function and convolution kernel.
    void ShiftResponseTime(double ticks);
    void SetPeakResponseTime(double tick);

    //Add DeconvKernel Polarity switch to decide how to normalize
    //deconvoluted signal w.r.t. RawDigits. If +1 then normalize
    //to Max ADC, if -1 to Min ADC
    void SetDeconvKernelPolarity(int pol);

    // Calculate convolution and deconvolution kernels 
    void CalculateConvKernel();
    void CalculateDeconvKernel(const std::vector<TComplex>& filterfunc);

  private:

    // Overall response.
    std::vector<float> fResponse;

    // Convolution kernel (fourier transform of response function).
    std::vector<ComplexF> fConvKernel;

    // Deconvolution kernel (= fFilter / fConvKernel).
    std::vector<ComplexF> fDeconvKernel;

    // Deconvolution Kernel Polarity Flag
    // Set to +1 if deconv signal should be deconv to + ADC count
    // Set to -1 if one wants to normalize to - ADC count
    int fDeconvKernelPolarity;

    // Xin added
    bool fNorm; 
  };
}

//----------------------------------------------------------------------
// Convolute a time series with current response.
template <class T> inline void util::SignalShapingLite::Convolute(std::vector<T>& func) const
{

  // Get FFT service.
  art::ServiceHandle<util::LArFFT> fft;

  // Make sure that time series has the correct size.
  if(func.size() != (size_t)fft->FFTSize())
    throw cet::exception("SignalShaping") << "Bad time series size = " << func.size() << "\n";

  //make a temporary std::vector<TComplex> convolution kernel
  std::vector<TComplex> temp_kernel(fConvKernel.size());
  for (unsigned int i=0; i!=temp_kernel.size(); ++i) {
    temp_kernel[i] = TComplex(fConvKernel[i].Re, fConvKernel[i].Im);
  }

  // Do convolution.
  fft->Convolute(func, const_cast<std::vector<TComplex>&>(temp_kernel));
}

//----------------------------------------------------------------------
// Convolute a time series with deconvolution kernel.
template <class T> inline void util::SignalShapingLite::Deconvolute(std::vector<T>& func) const
{

  // Get FFT service.
  art::ServiceHandle<util::LArFFT> fft;

  // Make sure that time series has the correct size.
  if(func.size() != (size_t)fft->FFTSize())
    throw cet::exception("SignalShaping") << "Bad time series size = " << func.size() << "\n";

  //make a temporary std::vector<TComplex> convolution kernel
  std::vector<TComplex> temp_kernel(fDeconvKernel.size());
  for (unsigned int i=0; i!=temp_kernel.size(); ++i) {
    temp_kernel[i] = TComplex(fDeconvKernel[i].Re, fDeconvKernel[i].Im);
  }

  // Do convolution.
  fft->Convolute(func, temp_kernel);
}




#endif
