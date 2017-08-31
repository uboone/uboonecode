//////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingLite.cxx
///
/// \brief  Generic signal shaping class.
///
/// \author Brandon Eberly 
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "cetlib/exception.h"
#include "SignalShapingLite.h"


//----------------------------------------------------------------------
// Constructor.
//
util::SignalShapingLite::SignalShapingLite() 
  : fNorm (true)
{}


//----------------------------------------------------------------------
// Destructor.
//
util::SignalShapingLite::~SignalShapingLite()
{}


//----------------------------------------------------------------------
// Reset this class to its default-constructed state.
void util::SignalShapingLite::Reset()
{
  fResponse.clear();
  fResponse.shrink_to_fit();
  fConvKernel.clear();
  fConvKernel.shrink_to_fit();
  fDeconvKernel.clear();
  fDeconvKernel.shrink_to_fit();
  
  //Set deconvolution polarity to + as default
  fDeconvKernelPolarity = +1;
}


//----------------------------------------------------------------------
// Add a time domain response function.
void util::SignalShapingLite::AddResponseFunction(const std::vector<float>& resp, bool ResetResponse )
{
  if (!fConvKernel.empty()) {
    std::cout<<"SignalShapingLite: Can't call AddResponseFunction if the convolution kernel is already configured"<<std::endl;
    return;
  }

  // Get FFT service.
  art::ServiceHandle<util::LArFFT> fft;
  int nticks = fft->FFTSize();

  // Copy new response function into fResponse attribute, and pad or
  // truncate to correct size.
  if (fResponse.size() == 0 || ResetResponse) {
    fResponse = resp;
    fResponse.resize(nticks, 0.);
  }
  else { 
    std::vector<float> tmp_resp = resp;   
    if (resp.size() != fResponse.size()) tmp_resp.resize(nticks, 0.);    
    fft->Convolute(fResponse, tmp_resp);
  }
}


//----------------------------------------------------------------------
// Shift the convolution kernel by the specified
// number of ticks.
void util::SignalShapingLite::ShiftResponseTime(double ticks)
{

  if (!fConvKernel.empty()) {

    // Get FFT service.
    art::ServiceHandle<util::LArFFT> fft;
    
    //make temporary std::vector<TComplex>
    std::vector<TComplex> tmp(fConvKernel.size());
    for (unsigned int i=0; i!=tmp.size(); ++i) {
      tmp[i] = TComplex(fConvKernel[i].Re, fConvKernel[i].Im);
    }
    
    // Update convolution kernel.
    fft->ShiftData(tmp, ticks);
    
    fConvKernel.clear();
    fConvKernel.resize(tmp.size());
    for (unsigned int i=0; i!=tmp.size(); ++i) {
      fConvKernel[i].Re = tmp[i].Re();
      fConvKernel[i].Im = tmp[i].Im();
    }
    
  }
}


//----------------------------------------------------------------------
// Set the peak response time to be at the specified tick.
void util::SignalShapingLite::SetPeakResponseTime(double tick)
{

  if (!fConvKernel.empty()) {
    // Get FFT service.
    art::ServiceHandle<util::LArFFT> fft;

    // Construct a delta-function response centered at tick zero.
    std::vector<double> delta(fft->FFTSize(), 0.);
    delta[0] = 1.;
    
    std::vector<double> tmp_response(fft->FFTSize());
    std::vector<TComplex> tmp_convkernel(fConvKernel.size());
    for (unsigned int i=0; i!=tmp_convkernel.size(); ++i) {
      tmp_convkernel[i] = TComplex(fConvKernel[i].Re, fConvKernel[i].Im);
    }  
    fft->DoInvFFT(tmp_convkernel, tmp_response);

    // Figure out peak of current overall response.
    double peak = fft->PeakCorrelation(delta, tmp_response);

    // Shift peak response to desired tick.
    ShiftResponseTime(tick - peak);
  }
}

//----------------------------------------------------------------------
// Add a DeconvKernel Polarity Flag to decide how to normalize
void util::SignalShapingLite::SetDeconvKernelPolarity(int pol)
{

  if ( (pol != 1) and (pol != -1) ) {
    throw cet::exception("SignalShaping") << __func__
      << ": DeconvKernelPolarity should be +1 or -1 (got " << pol << "). Setting to +1\n";
    fDeconvKernelPolarity = +1;
  }

  else
    fDeconvKernelPolarity = pol;

}


//----------------------------------------------------------------------
// Calculate the convolution kernel
void util::SignalShapingLite::CalculateConvKernel()
{
  if (fConvKernel.empty() &&!fResponse.empty()) {

    // Get FFT service.
    art::ServiceHandle<util::LArFFT> fft;
    unsigned int nticks = fft->FFTSize();
    
    if (fResponse.size() != nticks) {
      throw cet::exception("SignalShaping") << __func__ << ": inconsistent kernel size, "
        << fResponse.size() << " vs. " << nticks;
    }

    fConvKernel.resize(nticks/2 + 1);
    std::vector<TComplex> tmp(fConvKernel.size());
    
    fft->DoFFT(fResponse, tmp);
    for (unsigned int i=0; i!=tmp.size(); ++i) {
      fConvKernel[i].Re = tmp[i].Re();
      fConvKernel[i].Im = tmp[i].Im();
    }
    
    //try to save some memory, since the response is no longer needed
    fResponse.clear(); 
    fResponse.shrink_to_fit(); //implementation-dependent: may not actually do a reallocation (tested and works on the gpvms)
  }
}

//----------------------------------------------------------------------
// Calculate the deconvolution kernel as the ratio
// of the filter function and convolution kernel.
void util::SignalShapingLite::CalculateDeconvKernel(const std::vector<TComplex>& filterfunc)
{

  if (!fDeconvKernel.empty() || fConvKernel.empty()) return;

 
  // Get FFT service.
  art::ServiceHandle<util::LArFFT> fft;
  unsigned int n = fft->FFTSize();

  //Resize filter
  std::vector<TComplex> tmp_filter = filterfunc;
  tmp_filter.resize(n/2 + 1);

  // Make sure filter function has the correct size.
  // (Should always be the case if we get here.)
  if (tmp_filter.size() != fConvKernel.size()) {
    throw cet::exception("SignalShaping") << __func__ << ": inconsistent size, "
      << tmp_filter.size() << " vs. " << fConvKernel.size() << "\n";
  }
  
  // Calculate deconvolution kernel as the ratio of the 
  // filter function and the convolution kernel.
  std::vector<TComplex> tmp_deconvkernel = tmp_filter;
  for(unsigned int i=0; i<tmp_deconvkernel.size(); ++i) {
    if(std::abs(fConvKernel[i].Re) <= 0.0001 && std::abs(fConvKernel[i].Im) <= 0.0001) {
      tmp_deconvkernel[i] = 0.; 
    }
    else {
      TComplex val(fConvKernel[i].Re, fConvKernel[i].Im);
      tmp_deconvkernel[i] /= val; 
    }
  }


  // Normalize the deconvolution kernel.
  // Calculate the unnormalized deconvoluted response
  // (inverse FFT of filter function).
  std::vector<double> deconv(n, 0.);
  fft->DoInvFFT(const_cast<std::vector<TComplex>&>(tmp_filter), deconv);
  if (fNorm){
    std::vector<float> tmp_response(n);
    std::vector<TComplex> tmp_convkernel(n/2+1);
    for (unsigned int i=0; i!=tmp_convkernel.size(); ++i) {
      tmp_convkernel[i] = TComplex(fConvKernel[i].Re, fConvKernel[i].Im);
    }  
    fft->DoInvFFT(tmp_convkernel, tmp_response);
  
    // Find the peak value of the response
    // Should normally be at zero, but don't assume that.
    // Use DeconvKernelPolairty to find what to normalize to
    double peak_response = 0;
    if ( fDeconvKernelPolarity == -1 )
      peak_response = 4096;
    for(unsigned int i = 0; i < tmp_response.size(); ++i) {
      if( (tmp_response[i] > peak_response) 
	  and (fDeconvKernelPolarity == 1))
	peak_response = tmp_response[i];
      else if ( (tmp_response[i] < peak_response)
		and ( fDeconvKernelPolarity == -1) )
	peak_response = tmp_response[i];
    }
    if ( fDeconvKernelPolarity == -1 )
      peak_response *= -1;
    if (peak_response <= 0.) {
      throw cet::exception("SignalShaping") << __func__
					    << ": peak should always be positive (got " << peak_response << ")\n";
    }
    
    // Find the peak value of the deconvoluted response
    // Should normally be at zero, but don't assume that.
    
    double peak_deconv = 0.;
    for(unsigned int i = 0; i < deconv.size(); ++i) {
      if(deconv[i] > peak_deconv)
	peak_deconv = deconv[i];
    }
    if (peak_deconv <= 0.) {
      throw cet::exception("SignalShaping") << __func__
					    << ": deconvolution peak should always be positive (got " << peak_deconv << ")\n";
    }
    
    // Multiply the deconvolution kernel by a factor such that
    // (Peak of response) = (Peak of deconvoluted response).
    
    double ratio = peak_response / peak_deconv;
    for(unsigned int i = 0; i < tmp_deconvkernel.size(); ++i)
      tmp_deconvkernel[i] *= ratio;
  }//end if(fNorm)

  //copy over to fDeconvKernel
  fDeconvKernel.resize(tmp_deconvkernel.size());
  for (unsigned int i=0; i != tmp_deconvkernel.size(); ++i) {
    fDeconvKernel[i].Re = (float)tmp_deconvkernel[i].Re();
    fDeconvKernel[i].Im = (float)tmp_deconvkernel[i].Im();
  }

}
