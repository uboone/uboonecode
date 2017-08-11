///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceMicroBooNE.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution)./Users/Lsrea/newSim/SignalShapingServiceMicroBooNE.h
///
/// \author H. Greenlee, major mods by L. Rochester
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response (generated from the histogram).
/// Col3DCorrection - 3D path length correction for collection plane. (not invoked)
/// Ind3DCorrection - 3D path length correction for induction plane.  (not invoked)
/// FieldRespAmpVec - vector of response amplitudes, one for each view
/// ShapeTimeConst  - Time constants for exponential shaping.
/// FilterVec       - vector of filter function parameters, one for each view
/// FilterParamsVec - Vector of filter function parameters.
///
/// \update notes: Leon Rochester (lsrea@slac.stanford.edu, Jan 12, 2015
///                many changes, need to be documented better
///                 1. the three (or n) views are now represented by a vector of views
///                 2. There are separate SignalShaping objects for convolution and
///                    deconvolution
///
///                Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014
///                 1. Read in field responses from input histograms
///                 2. Allow different sampling rates in the input
///                    field response
///                    NOTE: The InputFieldRespSamplingRate parameter has
///                    NOT implemented for the field response input
///                    as a function (UseFunctionFieldShape)
///                 3. Allow different electron drift velocities from
///                    which the input field responses are obtained
///                 4. Convolute the field and electronic responses,
///                    and then sample the convoluted function with
///                    the nominal sampling rate (detinfo::DetectorPropertiesService).
///                    NOTE: Currently this doesn't include the filter 
///                    function and the deconvolution kernel.
///                    We may want to include them later?
///                 5. Disable fColSignalShaping.SetPeakResponseTime(0.),
///                    so that the peak time in the input field response
///                    is preserved.
///                 6. Somebody needs to unify the units of time (microsec
///                    or nanosec); I'm fainting!
///
/// New function:   void SetResponseSampling();
///
/// Modified functions: void init();
///                     void SetFieldResponse();
///
/// New FCL parameters:
/// DefaultDriftVelocity       - The electron drift velocity used to obtain
///                              the input field response waveforms
/// InputFieldRespSamplingRate - The sampling rate in the input field response
/// UseHistogramFieldShape     - Use the field response from an input histogram,
///                              if both UseFunctionFieldShape and 
///                              UseHistogramFieldShape are false, we will
///                              use the toy field responses (a bipolar square
///                              function for induction planes, a ramp function
///                              for collection planes.)
/// FieldResponseFname         - Name of the file containing the input field 
///                              response histograms
/// FieldResponseHistoName     - Name of the field response histograms,
///                              the format in the code will be 
///                              FieldResponseHistoName_U(V,Y)
///update notes: Jyoti Joshi (jjoshi@bnl.gov), Jan 13, 2015 
//               1. Modification to GetShapingTime function to read in different
//                  shaping time for different planes
//               2. Modification to GetASICGain fucntion to read in different gain 
//                  settings for different planes    
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGSERVICEMICROBOONE_H
#define SIGNALSHAPINGSERVICEMICROBOONE_H

#include <vector>
#include <map>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "SignalShapingLite.h"
#include "TF1.h"
#include "TH1D.h"

// LArSoft include
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

using DoubleVec  = std::vector<double>;
using FloatVec   = std::vector<float>;
using DoubleVec2 = std::vector< DoubleVec >;


namespace util {

  // helper class
  class SignalShapingContainer {
  
    public:
    
      SignalShapingContainer() {
        this->ResetAll();
      }
      
      ~SignalShapingContainer() {}
      
      SignalShapingLite& Response(const std::string& response_name) {
        return fResponseMap[response_name];
      }
      
      const SignalShapingLite& Response(const std::string& response_name) const {
        return fResponseMap.at(response_name);
      }
      
      void ResetAll() {
        for (auto& resp : fResponseMap) {
	  resp.second.Reset();
	}
      }
      
    private:
      
      std::map<std::string, SignalShapingLite> fResponseMap;
  };
  
  // little helper class to hold the parameters for charge deposition
  class ResponseParams {
    public:
      ResponseParams(double charge, double y, double z, size_t index) : m_charge(charge), m_y(y), m_z(z), m_index(index) {}
      double getCharge() { return m_charge; }
      double getY() { return m_y; }
      double getZ() { return m_z; }
      size_t getIndex()   { return m_index; }
    private:
      double m_charge;
      double m_y;
      double m_z;
      size_t m_index;
  };

  class SignalShapingServiceMicroBooNE {
 
    public:

      //------------------------------------------------------------
      // Constructor, destructor, and configuration.
      //------------------------------------------------------------
      SignalShapingServiceMicroBooNE(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);   
      ~SignalShapingServiceMicroBooNE();

      void reconfigure(const fhicl::ParameterSet& pset);


      //------------------------------------------------------------
      // Convolution (for simulation) and Deconvolution (for reconstruction)
      //------------------------------------------------------------

      // Do convolution calculation (for simulation).
      template <class T> void Convolute(size_t channel, std::vector<T>& func, double y, double z) const;
      template <class T> void Convolute(size_t channel, std::vector<T>& func, const std::vector<std::unique_ptr<ResponseParams> >& params) const;
      template <class T> void Convolute(size_t channel, std::vector<T>& func, const std::string& response_name) const;

      // Do deconvolution calculation (for reconstruction).
      template <class T> void Deconvolute(size_t channel, std::vector<T>& func, double y, double z);
      template <class T> void Deconvolute(size_t channel, std::vector<T>& func, const std::string& response_name);


      //------------------------------------------------------------
      // Accessors.
      //------------------------------------------------------------

      // Responses and Kernels
      const util::SignalShapingLite& SignalShaping(size_t channel, const std::string& response_name="") const;
      const std::vector<ComplexF>&   GetConvKernel(unsigned int channel, const std::string& response_name ) const;  // M. Mooney 
      void  SetDecon(size_t datasize);
      const std::vector<unsigned int>&   GetNResponses() const { return fNResponses; } 
      
      int FieldResponseTOffset(unsigned int const channel, const std::string& response_name) const;


      // Filter-related    
      double Get2DFilterVal(size_t planeNum, size_t freqDimension, double binFrac) const;  // M. Mooney
      double Get2DFilterNorm(size_t planeNum) const;  // M. Mooney


      // Noise-related
      const DoubleVec2& GetNoiseFactVec() const { return fNoiseFactVec; }
      double     GetDeconNorm() const { return fDeconNorm; }
      double     GetRawNoise(unsigned int const channel) const;
      double     GetDeconNoise(unsigned int const channel) const; 
      
      
      //U- and Y-Shorted Channels and YZ-Regions
      // Determine if YZ coordinate overlaps with the U-shorted YZ region
      static bool IsUShortedYZRegion(double y, double z);
      
      // Determine if Z coordinate overlaps with the Y-shorted Z region
      static bool IsYShortedZRegion(double z);
      
      // Determine if channel's wire overlaps with the U-shorted wires
      static bool IsUShortedOverlapChannel(unsigned int channel);
      
      // Determine if channel's wire overlaps with the Y-shorted wires
      static bool IsYShortedOverlapChannel(unsigned int channel);

    private:

      //------------------------------------------------------------
      // Private configuration methods.
      //------------------------------------------------------------

      // Post-constructor initialization.
      void init() const{const_cast<SignalShapingServiceMicroBooNE*>(this)->init();}
      void init();

      void SetFieldResponse();
      void SetElectResponse();  //changed to read different peaking time for different planes

      // Calculate filter functions.
      void SetFilters();

      // Sample the response function, including a configurable
      // drift velocity of electrons
      void SetResponseSampling();

      // Get the name of the (possibly YZ-dependent) response to use, 
      // as well as a charge_fraction for scaling before convolution/deconvolution
      std::string DetermineResponseName(unsigned int chan, double y, double z, double& charge_fraction) const;
      
      // Determine whether the response indicated by response_name should be stored for channel
      bool StoreThisResponse(const std::string& response_name, unsigned int channel) const;
      


      //------------------------------------------------------------
      // Private Attributes.
      //------------------------------------------------------------     

      bool fInitForConvolution;       ///< True if initialized for convolution
      bool fInitForDeconvolution;     ///< True if initialized for deconvolution
      bool fSetDeconKernelsUponInit;  ///< If true, deconvolution kernels are calculated and set in init().  Otherwise, must use SetDecon()


      //Convolution and Deconvolution
      std::vector<util::SignalShapingContainer> fSignalShapingVec; ///< Stores the convolution and deconvolution kernels   
      unsigned int fConvFFTSize;                                   ///< Stores the FFTSize used for convolution 
      std::vector<int>                          fDeconvPol;        ///< switch for DeconvKernel normalization sign (+ -> max pos ADC, - -> max neg ADC). Entry 0,1,2 = U,V,Y plane settings


      // Electronics Response    
      std::vector<FloatVec>  fElectResponse;            ///< Stores electronics response. Vector elements correspond to channel, then response bins
      unsigned int           fMaxElectResponseBins;     ///< Maximum number of bins to use to store the electronics response
      double                 fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC
      bool                   fIgnoreMisconfigStatus;    ///< If true, use nominal (default) gain and shaping time to calculate electronics response of misconfigured channels.  Useful for deconvolution if noise filtering already performed the correction
      

      // Field Response
      std::vector< std::map<std::string, FloatVec> >  fFieldResponseVec;      ///< Stores adjusted field response. Vector elements corespond to plane, map key is a response name
      std::vector< std::map<std::string, TH1F*> >     fFieldResponseHistVec;  ///< Stores input field response.  Vector elements correspond to planes, map key is a response name
      std::vector<unsigned int>                       fNResponses;            ///< Number of input field responses per view
      DoubleVec                                       fFieldRespAmpVec;       ///< Amplitudes applied to adjusted field response   
      bool                                            fYZdependentResponse;   ///< Using YZ-dependent responses
      bool                                            fdatadrivenResponse;    ///< Using data-driven responses
      std::vector< std::map<std::string, double> >    fFieldResponseTOffset;  ///< Time offset for field response in ns. Vector elements correspond to plane, map key is a response name
      DoubleVec  	      	      	              fCalibResponseTOffset;  ///< calibrated time offset to align U/V/Y Signals												 
 													 
 												 


      // Filter Parameters
      bool                                fGetFilterFromHisto;    ///< Flag that allows to use a filter function from a histogram instead of the functional dependency
      std::vector<TH1D*>                  fFilterHistVec;
      std::vector<TF1*>                   fFilterTF1Vec;          ///< Vector of Parameterized filter functions
      std::vector<std::string>            fFilterFuncVec;
      std::vector<std::vector<TComplex> > fFilterVec;
      DoubleVec2                          fFilterParamsVec;
      DoubleVec                           fFilterWidthCorrectionFactor;  // a knob

      // Filter Parameters - Induced charge deconvolution additions (M. Mooney)
      std::vector<TF1*>        fFilterTF1VecICTime;
      std::vector<std::string> fFilterFuncVecICTime;
      std::vector<TF1*>        fFilterTF1VecICWire;
      std::vector<std::string> fFilterFuncVecICWire;
      DoubleVec                fFilterScaleVecICTime;
      DoubleVec                fFilterScaleVecICWire;
      DoubleVec                fFilterNormVecIC;
      std::vector<double>      fFilterICTimeMaxFreq;
      DoubleVec                fFilterICTimeMaxVal;
      DoubleVec                fFilterICWireMaxFreq;
      DoubleVec                fFilterICWireMaxVal;


      // Noise
      DoubleVec2 fNoiseFactVec;       ///< RMS noise in ADCs for lowest gain setting
      double     fDeconNorm;          ///< Set Decon Noise Scale


      // Diagnostics
      int fDiagnosticChannel;
      std::string fDiagnosticResponse;
      
      bool fPrintResponses;
      bool fHistDoneF[3];   
      TH1D* fHFullResponse[3];
      TH1D* fHSampledResponse[3];
      
      TH1D* fHist_FieldResponseHist;
      TH1D* fHist_FieldResponseVec;
      TH1D* fHist_ElectResponse;
      TH1D* fHist_TotalResponse;
      TH1D* fHist_ResampledConvKernelRe;
      TH1D* fHist_ResampledConvKernelIm;
      
      TH1D* fHist_PreConv;
      TH1D* fHist_PostConv;
      TH1D* fHist_PostOffset;
      TH1D* fHist_PreDeconv;
      TH1D* fHist_PostDeconvOffset;
    
  };
} //namespace util



//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(size_t channel, 
									       std::vector<T>& func, 
									       const std::vector<std::unique_ptr<util::ResponseParams> >& params) const
{

  //loop over params and construct input functions for each of the channel's responses
  std::map<std::string, std::vector<T> > input_map;
  for (const auto& item : params) {
    double charge_fraction = 1.0;
    std::string response_name = this->DetermineResponseName(channel, item->getY(), item->getZ(), charge_fraction);
    
    if (input_map.find(response_name)==input_map.end()) {
      input_map[response_name].resize(func.size(), 0.0);
    }
    
    if (item->getIndex() < 0 || item->getIndex() >= func.size()) continue;
    input_map[response_name].at(item->getIndex()) += item->getCharge()*charge_fraction;
  }
  
  //convolute each input function, and add to output
  std::vector<T> output;
  output.resize(func.size(), 0.0);
  for (auto input : input_map) {
    this->Convolute(channel, input.second, input.first);
    std::transform(input.second.begin(), input.second.end(), output.begin(), output.begin(), std::plus<T>()); 
  }
  
  //copy to func
  func = output;
}
    
      
//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(size_t channel, std::vector<T>& func, double y, double z) const
{

  double charge_fraction = 1.0;
  std::string response_name = this->DetermineResponseName(channel, y, z, charge_fraction);
  if (charge_fraction != 1.0) {
    for (auto& element : func) {
      element *= charge_fraction;
    }
  }

  this->Convolute(channel, func, response_name); 
}

//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(size_t channel, std::vector<T>& func, const std::string& response_name) const
{

  init();
  
  if ( !StoreThisResponse(response_name, channel) ) {
    throw cet::exception(__FUNCTION__) << "You requested to use an invalid response "<<response_name<<" for channel "<<channel<<std::endl;
  }
  
  //make diagnostic histogram
  if ((int)channel==fDiagnosticChannel && response_name==fDiagnosticResponse) {
    for (unsigned int bin=0; bin!=func.size(); ++bin) {
      fHist_PreConv->SetBinContent(bin+1, func.at(bin));
    }
  }
  
  //convolute
  fSignalShapingVec[channel].Response(response_name).Convolute(func);

  //make diagnostic histogram
  if ((int)channel==fDiagnosticChannel && response_name==fDiagnosticResponse) {
    for (unsigned int bin=0; bin!=func.size(); ++bin) {
      fHist_PostConv->SetBinContent(bin+1, func.at(bin));
    }
  }
  
  //Add time offset
  int time_offset = FieldResponseTOffset(channel,response_name);  
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.begin(),func.begin()-time_offset);
    func.erase(func.begin(),func.begin()-time_offset);
    func.insert(func.end(),temp.begin(),temp.end());
  }
  else{
    temp.assign(func.end()-time_offset,func.end());
    func.erase(func.end()-time_offset,func.end());
    func.insert(func.begin(),temp.begin(),temp.end());
  }
  
  //make diagnostic histogram
  if ((int)channel==fDiagnosticChannel && response_name==fDiagnosticResponse) {
    for (unsigned int bin=0; bin!=func.size(); ++bin) {
      fHist_PostOffset->SetBinContent(bin+1, func.at(bin));
    }
  }
  
}


//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Deconvolute(size_t channel, std::vector<T>& func, double y, double z)
{
  
  double charge_fraction = 1.0;
  std::string response_name = this->DetermineResponseName(channel, y, z, charge_fraction);
  if (charge_fraction != 1.0 && charge_fraction > 0.0) {
    for (auto& element : func) {
      element /= charge_fraction;
    }
  }

  this->Deconvolute(channel, func, response_name); 
}  
  

//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Deconvolute(size_t channel, std::vector<T>& func, const std::string& response_name)
{

  init();
  
  if ( !StoreThisResponse(response_name, channel) ) {
    throw cet::exception(__FUNCTION__) << "You requested to use an invalid response "<<response_name<<" for channel "<<channel<<std::endl;
  }
  
  //Initialize deconvolution kernels (necessary if this wasn't done in init())
  if (!fInitForDeconvolution) {
    this->SetDecon(func.size());
  }
  
  //make diagnostic histogram
  if ((int)channel==fDiagnosticChannel && response_name==fDiagnosticResponse) {
    for (unsigned int bin=0; bin!=func.size(); ++bin) {
      fHist_PreDeconv->SetBinContent(bin+1, func.at(bin));
    }
  }
  
  //Do deconvolution
  fSignalShapingVec[channel].Response(response_name).Deconvolute(func);

  //Add Time Offset
  int time_offset = FieldResponseTOffset(channel,response_name); 
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.end()+time_offset,func.end());
    func.erase(func.end()+time_offset,func.end());
    func.insert(func.begin(),temp.begin(),temp.end());
  }
  else{
    temp.assign(func.begin(),func.begin()+time_offset);
    func.erase(func.begin(),func.begin()+time_offset);
    func.insert(func.end(),temp.begin(),temp.end());   
  }
  
  //make diagnostic histogram
  if ((int)channel==fDiagnosticChannel && response_name==fDiagnosticResponse) {
    for (unsigned int bin=0; bin!=func.size(); ++bin) {
      fHist_PostDeconvOffset->SetBinContent(bin+1, func.at(bin));
    }
  }

}

DECLARE_ART_SERVICE(util::SignalShapingServiceMicroBooNE, LEGACY)
#endif
