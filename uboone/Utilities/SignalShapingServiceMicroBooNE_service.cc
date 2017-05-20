
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceMicroBooNE_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "TFile.h"

#include <fstream>

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceMicroBooNE::SignalShapingServiceMicroBooNE(const fhicl::ParameterSet& pset,
                                                                     art::ActivityRegistry& /* reg */)
: fInit(false), fInitConfigMap(false)
{
  for(size_t i=0; i<3; ++i) {
    fHRawResponse[i] = 0;
    fHStretchedResponse[i] = 0;
    fHFullResponse[i] = 0;
    fHistDone[i] = false;
    fHistDoneF[i] = false;
  }
    
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::SignalShapingServiceMicroBooNE::~SignalShapingServiceMicroBooNE()
{}


//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceMicroBooNE::reconfigure(const fhicl::ParameterSet& pset)
{
  //Services
  art::ServiceHandle<geo::Geometry> geo;
 
  // Reset initialization flag.
  fInit = false;
  fInitConfigMap = false;

  //two configs - misconfigured and not
  fNConfigs = 2;
  
  fNViews = geo->Nplanes();
  
  //fetch fcl parameters
  fViewForNormalization = pset.get<size_t>("ViewForNormalization");
  fNResponses       = pset.get<std::vector<std::vector<size_t> > >("NResponses");
  fNYZResponses     = pset.get<std::vector<std::vector<size_t> > >("NYZResponses");
  fNdatadrivenResponses = pset.get<std::vector<std::vector<size_t> > >("NdatadrivenResponses");
  fNActiveResponses = pset.get<std::vector<std::vector<size_t> > >("NActiveResponses");
  fNYZActiveResponses = pset.get<std::vector<std::vector<size_t> > >("NYZActiveResponses");
  fNdatadrivenActiveResponses = pset.get<std::vector<std::vector<size_t> > >("NdatadrivenActiveResponses");

  fPrintResponses   = pset.get<bool>("PrintResponses");
  fYZdependentResponse = pset.get<bool>("YZdependentResponse");
  fdatadrivenResponse = pset.get<bool>("datadrivenResponse");
  fIncludeMisconfiguredU = pset.get<bool>("IncludeMisconfiguredU");

  fDeconNorm = pset.get<double>("DeconNorm");
  fASICGainInMVPerFC    = pset.get< DoubleVec2 >("ASICGainInMVPerFC");
  fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");

  fNoiseFactVec       = pset.get<DoubleVec2>("NoiseFactVec");
  f3DCorrectionVec    = pset.get<DoubleVec>("Drift3DCorrVec");
  fFieldRespAmpVec    = pset.get<DoubleVec>("FieldRespAmpVec");
  fShapeTimeConst     = pset.get<DoubleVec2 >("ShapeTimeConst");
  fDeconvPol          = pset.get<std::vector<int> >("DeconvPol");
  
  fDefaultEField                 = pset.get<double>("DefaultEField");
  fDefaultTemperature            = pset.get<double>("DefaultTemperature");
  fTimeScaleParams               = pset.get<DoubleVec>("TimeScaleParams");
  fStretchFullResponse           = pset.get<bool>("StretchFullResponse");
  
  fCalibResponseTOffset = pset.get< DoubleVec >("CalibResponseTOffset");
  
  // Construct parameterized collection filter function.
  fGetFilterFromHisto = pset.get<bool>("GetFilterFromHisto");  
  fFilterWidthCorrectionFactor = pset.get<DoubleVec>("FilterWidthCorrectionFactor", DoubleVec() = {1.0, 1.0, 1.0});   
  if(!fGetFilterFromHisto) {

    fFilterFuncVec.resize(fNViews);
    mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting Filters from .fcl file" ;

    fFilterParamsVec = pset.get< DoubleVec2 >("FilterParamsVec");
    fFilterFuncVec = pset.get<std::vector<std::string> > ("FilterFuncVec");

    fFilterTF1Vec.resize(fNViews);
    for(unsigned int _vw=0;_vw<fNViews; ++_vw) {
      std::string name = Form("Filter_vw%02i_wr%02i", (int)_vw, 0);
      fFilterTF1Vec[_vw] = new TF1(name.c_str(), fFilterFuncVec[_vw].c_str() );
      for(unsigned int _ind=0; _ind<fFilterParamsVec[_vw].size(); ++_ind) {
        fFilterTF1Vec[_vw]->SetParameter(_ind, fFilterParamsVec[_vw][_ind]);
      }
    }
  } else {

    std::string histoname = pset.get<std::string>("FilterHistoName");
    mf::LogInfo("SignalShapingServiceMicroBooNE") << " using filter from .root file " ;

    // constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);

    TFile * in=new TFile(fname.c_str(),"READ");
    for(unsigned int _vw=0;_vw<fNViews;_vw++){
      std::string name = Form("%s_vw%02i", histoname.c_str(), (int)_vw);
      fFilterHistVec[_vw] = (TH1D *)in->Get(name.c_str());
    }

    in->Close();
    delete in;
  }
  
  
  // Load 2D filters for induced charge deconvolution (M. Mooney)
  fFilterFuncVecICTime.resize(fNViews);
  fFilterFuncVecICWire.resize(fNViews);
  mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting 2D Filters from .fcl file" ;

  DoubleVec2 paramsICTime = pset.get< DoubleVec2 >("FilterParamsVecICTime");
  fFilterFuncVecICTime = pset.get<std::vector<std::string> > ("FilterFuncVecICTime");
  
  DoubleVec2 paramsICWire = pset.get< DoubleVec2 >("FilterParamsVecICWire");
  fFilterFuncVecICWire = pset.get<std::vector<std::string> > ("FilterFuncVecICWire");

  fFilterTF1VecICTime.resize(fNViews);
  fFilterICTimeMaxFreq.resize(fNViews);
  fFilterICTimeMaxVal.resize(fNViews);
  
  fFilterTF1VecICWire.resize(fNViews);
  fFilterICWireMaxFreq.resize(fNViews);
  fFilterICWireMaxVal.resize(fNViews);
  for(unsigned int _vw=0;_vw<fNViews; ++_vw) {
    std::string name = Form("FilterICTime_vw%02i_wr%02i", (int)_vw, 0);
    fFilterTF1VecICTime[_vw] = new TF1(name.c_str(), fFilterFuncVecICTime[_vw].c_str() );
    
    name = Form("FilterICWire_vw%02i_wr%02i", (int)_vw, 0);
    fFilterTF1VecICWire[_vw] = new TF1(name.c_str(), fFilterFuncVecICWire[_vw].c_str() );
    for(unsigned int _ind=0; _ind<paramsICTime[_vw].size(); ++_ind) {
      fFilterTF1VecICTime[_vw]->SetParameter(_ind, paramsICTime[_vw][_ind]);
      fFilterICTimeMaxFreq[_vw] = fFilterTF1VecICTime[_vw]->GetMaximumX();
      fFilterICTimeMaxVal[_vw] = fFilterTF1VecICTime[_vw]->GetMaximum();
    }
    for(unsigned int _ind=0; _ind<paramsICWire[_vw].size(); ++_ind) {
      fFilterTF1VecICWire[_vw]->SetParameter(_ind, paramsICWire[_vw][_ind]);
      fFilterICWireMaxFreq[_vw] = fFilterTF1VecICWire[_vw]->GetMaximumX();
      fFilterICWireMaxVal[_vw] = fFilterTF1VecICWire[_vw]->GetMaximum();
    }
  }
  
    
  // get the scale factor between the bulk drift velocity used to generate the field response
  //   and that used for this simulation.
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  double defaultVelocity = detprop->DriftVelocity(fDefaultEField, fDefaultTemperature);
  double thisVelocity    = detprop->DriftVelocity( detprop->Efield(0), detprop->Temperature() );
  double vRatio = defaultVelocity/thisVelocity;
  double vDiff = vRatio -1.0;

  fTimeScaleFactor = 0.0;
  double term = 1.0;

  // the time scale params are from a fit to Garfield simulations at different E Fields
  for(size_t i = 0;i<fTimeScaleParams.size(); ++i) {
    fTimeScaleFactor += fTimeScaleParams[i]*term;
    term *= vDiff;
  }

  std::cout << "Current E field = " << detprop->Efield(0) << " KV/cm, Ratio of drift velocities = " << vRatio << ", timeScaleFactor = " << fTimeScaleFactor << std::endl;

  //do some checks on fcl-set parameters
  for(size_t ktype=0; ktype<2; ++ktype) {
    bool badN = false;

    for(size_t ktype=0; ktype<2; ++ktype) {
      for(unsigned int _vw=0; _vw<fNViews; ++_vw) {
	if(!fYZdependentResponse) {
	  if(fNResponses[ktype][_vw] < fNActiveResponses[ktype][_vw])  {
	    badN = true;
	  }
	}	
	else {
	  if(fNResponses[ktype][_vw] > 1){
	    throw art::Exception( art::errors::Configuration ) << "Don't use DIC with YZ dependent responses enabled; make all elements of NResponses == 1" << std::endl;
	  }
	  if(fdatadrivenResponse){
	    if(fNdatadrivenResponses[ktype][_vw]<fNdatadrivenActiveResponses[ktype][_vw]){
	      badN = true;
	    }
	  }
	  else{
	    if(fNYZResponses[ktype][_vw] < fNYZActiveResponses[ktype][_vw]) {
	      badN = true;
	    }
	  }
	}
      }
    }
    if(badN) {
      throw art::Exception( art::errors::InvalidNumber )
	<< "check NResponses/NActiveRespones or NYZResponses/NYZActiveResponses or NdatadrivenResponses/NdatadrivenActiveResponses" << std::endl;
    }
  }


  // Reset kernels.
  fSignalShapingVec.resize(fNConfigs);
  fFieldResponseVec.resize(fNConfigs);
  for (unsigned int config=0; config!= fNConfigs; ++config) {
    fSignalShapingVec[config].resize(2);
    fFieldResponseVec[config].resize(2);
    for (unsigned int kset=0; kset!=2; ++kset) {
      fSignalShapingVec[config][kset].resize(fNViews);
      fFieldResponseVec[config][kset].resize(fNViews);
      size_t nWires = 0;
      for (unsigned int plane=0; plane!=fNViews; ++plane) {
        if(!fYZdependentResponse) {
	  nWires = fNResponses[kset][plane];
	}
	else {
	  if(fdatadrivenResponse){
            nWires = fNdatadrivenResponses[kset][plane]; 
	  }
	  else{
	    nWires = fNYZResponses[kset][plane];
	  }
	}
        fSignalShapingVec[config][kset][plane].resize(nWires);
	fFieldResponseVec[config][kset][plane].resize(nWires);
	for ( auto& ss : fSignalShapingVec[config][kset][plane] ) {
	  ss.Reset();
	}
      }//loop over planes
    }//loop over kset
  }//loop over config
  
  fFieldResponseTOffset.resize(2);
  for(size_t ktype=0;ktype<2;++ktype) {
    fFieldResponseTOffset[ktype].resize(fNViews);
  }
  
  
  // constructor decides if initialized value is a path or an environment variable
  mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the field response provided from a .root file " ;
  std::string fileNameBase = pset.get<std::string>("FieldResponseFNameBase");
  std::vector<std::string> version      = pset.get<std::vector<std::string> >("FieldResponseFVersion"); 
  std::string histNameBase = pset.get<std::string>("FieldResponseHNameBase");
  cet::search_path sp("FW_SEARCH_PATH");
  DoubleVec tOffset(fNViews, 0.0);
  
  fFieldResponseHistVec.resize(fNConfigs);
  for(size_t config=0;config<fNConfigs; ++config) {
    fFieldResponseHistVec[config].resize(2);
    for(size_t ktype=0;ktype<2;++ktype) {
      fFieldResponseHistVec[config][ktype].resize(fNViews);
      unsigned int _vw = 0;
      //std::cout << "config " << config << " in " << fNConfigs << std::endl;
      for(auto& plane : fFieldResponseHistVec[config][ktype]) {
	std::string fname0, fname;

	if (!fdatadrivenResponse) {
	  plane.resize(fNResponses[ktype][_vw]);
	}
	else {
	  plane.resize(fNdatadrivenResponses[ktype][_vw]);
	}


	// load up the response functions
	unsigned int _wr = 0;
	int ctr = 0;
	for(auto& resp : plane) {

          //get response function file name
	  if (fdatadrivenResponse) {
	    if( (ctr == 0) || (_vw == 0 && ctr == 2) ){ 
	      fname0 = Form("%s_vw%02i_nominal_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
	    }
	    if( (_vw == 0 && ctr == 1) || (_vw == 1 && ctr == 2) ){
	      fname0 = Form("%s_vw%02i_shortedY_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
	    }
	    if( (_vw == 1 && ctr == 1) || (_vw == 2 && ctr == 1) ){
	      fname0 = Form("%s_vw%02i_shortedU_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
	    }
	  }
	  else {
	    fname0 = Form("%s_vw%02i_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
	  }	  
          sp.find_file(fname0, fname);
	  std::unique_ptr<TFile> fin(new TFile(fname.c_str(), "READ"));
	  
	  TString histName = Form("%s_vw%02i_wr%02i", histNameBase.c_str(), (int)_vw, (int)_wr);
	  resp = (TH1F*)fin->Get(histName);
	  auto Xaxis = resp->GetXaxis();
	  fNFieldBins[ktype] = Xaxis->GetNbins();
	  fFieldLowEdge[ktype] = resp->GetBinCenter(1) - 0.5*resp->GetBinWidth(1);
	  fFieldBin1Center[ktype] = resp->GetBinCenter(1);
	  // internal time is in nsec
	  fFieldBinWidth[ktype] = resp->GetBinWidth(1)*1000.;
	  
	  // get the offsets for each plane... use wire 0 and either peak or zero-crossing
	  double tOffset = 0.0;

	  // this is for the standard response
	  if(_wr==0 && _vw==fViewForNormalization) {
	    // for the collection plane, find the peak
	    int binMax = resp->GetMaximumBin();
	    tOffset = (resp->GetXaxis()->GetBinCenter(binMax) - resp->GetXaxis()->GetBinCenter(1));
	    // for later, to be a bit cleverer, take weighted average of 3 bins
	    //          for(int bin=binMax-1; bin<=binMax+1; ++ bin) {
	    //            content = resp->GetBinContent(bin);
	    //            binVal = resp->GetXaxis()->GetBinCenter(bin);
	    //            numer += content*binVal;
	    //            denom += content;
	    //          }
	    //          tOffset[_vw] = numer/denom*delta - resp->GetXaxis()->GetBinCenter(1);
	  } else {
	    // for the other planes, find the zero-crossing
	    // lets find the minimum, and work backwards!

	    int binMin = resp->GetMinimumBin();
	    for(int bin=binMin;bin>0; --bin) {
	      double content = resp->GetBinContent(bin);
	      bool found = false;
	      if(content>0) {
        	double binVal = resp->GetXaxis()->GetBinCenter(bin);
        	tOffset = binVal - resp->GetXaxis()->GetBinCenter(1);
        	found = true;
        	// for later
        	//            } else if (content>0) {
        	//              // If it's already gone through zero, split the difference
        	//              std::cout << resp->GetBinContent(bin) << " " << resp->GetBinContent(bin+1) << " " << bin << std::endl;
        	//              binVal = resp->GetXaxis()->GetBinCenter(bin);
        	//              numer = resp->GetBinContent(bin)*binVal - resp->GetBinContent(bin+1)*(binVal+delta);
        	//              denom = resp->GetBinContent(bin) - resp->GetBinContent(bin+1);
        	//              found = true;
	      }
	      if(found) break;
	    }
	  }

	  tOffset *= f3DCorrectionVec[_vw]*fTimeScaleFactor;
	  fFieldResponseTOffset[ktype].at(_vw) = (-tOffset + fCalibResponseTOffset[_vw])*1000.;
	  
	  fin->Close();
	  
	  if (fdatadrivenResponse) ctr++;
	  else _wr++;
	}
	_vw++;

      }//end loop over plane
    }//end loop over ktype
  }//end loop over config

}


//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceMicroBooNE::init()
{

  //do nothing if already initialized
  if (fInit) return;

  fInit = true;

  //first set misconfigured asic information
  //Due to structure of this code, all events processed in a job will use the 
  //misconfigured channel status of the first event (because init() is only called 
  //once, presumably for the first event)
  art::ServiceHandle<geo::Geometry> geo;
  if (!fInitConfigMap) {
    /*const size_t N_CHANNELS = geo->Nchannels();
    const lariov::ElectronicsCalibProvider& elec_provider = art::ServiceHandle<lariov::ElectronicsCalibService>()->GetProvider();
    for (unsigned int channel=0; channel!= N_CHANNELS; ++channel) {
      if (elec_provider.ExtraInfo(channel).GetBoolData("is_misconfigured")) {
	fConfigMap[channel] = 1;
      }
    }*/
    fInitConfigMap = true;
  }

  // Do microboone-specific configuration of SignalShaping by providing
  // microboone response and filter functions.

  // re-initialize the FFT service for the request size
  art::ServiceHandle<util::LArFFT> fFFT;
  int fftsize = (int) fFFT->FFTSize();

  // Calculate field and electronics response functions.

  std::string kset[2] = { "Convolution ", "Deconvolution "};

  for(size_t config=0;config<fNConfigs;++config) {
    for(size_t ktype=0;ktype<2;++ktype) {
      if (fNFieldBins[ktype]*4>fftsize)
        fFFT->ReinitializeFFT( (size_t)fNFieldBins[ktype]*4, fFFT->FFTOptions(), fFFT->FFTFitBins());

      // call this first, so that the binning will be known to SetElectResponse
      SetFieldResponse(ktype);

      for(unsigned int _vw=0;_vw<fNViews; ++_vw) {
        
	SetElectResponse(ktype,fShapeTimeConst[config].at(_vw),fASICGainInMVPerFC[config].at(_vw));
        if(fPrintResponses) {
          for(size_t i = 0; i<100; ++i) {
            std::cout << "Electronic Response " << fElectResponse[ktype][i] << " " ;
            if((i+1)%10==0) std::cout << std::endl;
          }
          std::cout << std::endl;
        }
	
	if (fYZdependentResponse || !fdatadrivenResponse) {		
	  unsigned int NResponses = fNdatadrivenResponses[ktype][_vw];
	  if (!fdatadrivenResponse) {
	    NResponses = fYZdependentResponse ? fNYZResponses[ktype][_vw] : fNResponses[ktype][_vw];
	  }

	  bool YZflag = true;
	  for (unsigned int _wr=0; _wr < NResponses; ++_wr) {

	    if (YZflag == false) {
	      _wr = 0;
	      _vw = 2;
	    }

	    if(fPrintResponses) {          std::cout << "Input field response for view " << _vw << " wire " << _wr
						     << ", " << (fFieldResponseVec[config][ktype][_vw][_wr]).size() << " bins" << std::endl;
	      for(size_t i = 0; i<(fFieldResponseVec[config][ktype][_vw][_wr]).size(); ++i) {
		std::cout << fFieldResponseVec[config][ktype][_vw][_wr][i] << " " ;
		if((i+1)%10==0) std::cout << std::endl;
	      }
	      std::cout << std::endl;
	    }

	    (fSignalShapingVec[config][ktype][_vw][_wr]).AddResponseFunction(fFieldResponseVec[config][ktype][_vw][_wr]);
	    (fSignalShapingVec[config][ktype][_vw][_wr]).AddResponseFunction(fElectResponse[ktype]);
	    (fSignalShapingVec[config][ktype][_vw][_wr]).save_response();
	    (fSignalShapingVec[config][ktype][_vw][_wr]).set_normflag(false);

	    if (fYZdependentResponse && !fdatadrivenResponse) {
	      if (YZflag == false) {
		_wr = 2;
		_vw = 1;
              }
	      YZflag = false;
	    }
	  }//end loop over _wr
	}
	
      }//end loop over _vw
      
    
      // see if we get the same toffsets
      SetResponseSampling(ktype, config);

      // Currently we only have fine binning "fFieldBinWidth"
      // for the field and electronic responses.
      // Now we are sampling the convoluted field-electronic response
      // with the nominal sampling.
      // We may consider to do the same for the filters as well.
      if ((int)fftsize!=fFFT->FFTSize()){
        fFFT->ReinitializeFFT( (size_t)fftsize, fFFT->FFTOptions(), fFFT->FFTFitBins());
      }


      // Calculate filter functions.
      if(config==0) SetFilters();

      // Configure deconvolution kernels.
      for(unsigned int _vw=0;_vw<fNViews; ++_vw) {
        // std::cout << "filtervec size" << fFilterVec[_vw].size() << std::endl;
	//if(!fYZdependentResponse){
	  for(unsigned int _wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
	    (fSignalShapingVec[config][ktype][_vw][_wr]).AddFilterFunction(fFilterVec[_vw]);
	    (fSignalShapingVec[config][ktype][_vw][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(_vw));
	    (fSignalShapingVec[config][ktype][_vw][_wr]).CalculateDeconvKernel();
	  }
	  //}
	  // YZdependent devonolution disabled --> want to deconvolve with nominal run 1 responses
	  /*	  if(fYZdependentResponse){
	  bool YZflag = true;
          for(_wr=0; _wr<fNYZResponses[ktype][_vw]; ++_wr) {
	    if(YZflag == false){ 
	      _wr = 0; 
	      _vw = 2;
	    }
	    (fSignalShapingVec[config][ktype][_vw][_wr]).AddFilterFunction(fFilterVec[_vw]);
	    (fSignalShapingVec[config][ktype][_vw][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(_vw));
	    if(YZflag == true && _vw == 2 && _wr == 0){
	      fSignalShapingVec[config][ktype][_vw][_wr].AddResponseFunction(fSignalShapingVec[config][ktype][_vw][_wr].Response_save(), true);
	    }

	    (fSignalShapingVec[config][ktype][_vw][_wr]).CalculateDeconvKernel();
	    if(YZflag == false){ 
	      fSignalShapingVec[config][ktype][_vw][_wr].save_response();
	      fSignalShapingVec[config][ktype][_vw][_wr].Reset();
	      _wr = 2; 
	      _vw = 1;

	    }
            YZflag = false;
	  }
	  */
      }//end loop over _vw
    }//end loop over ktype
  }//end loop over config
}


//----------------------------------------------------------------------
// SetDeconvolution kernels on a channel-by-channel basis
// Does the following after calling init():
// 1) Reset all ktype=1 entries in fSignalShapingVec
// 2) Call SetResponseSampling for the given channel for each config, with mode=0 and ktype=1
// 3) Call SetFilters if config=0
// 4) For each config, view, wire and ktype=1 of fSignalShapingVec, do:
//   4a) AddFilterFunction
//   4b) SetDeconvKernelPolarity
//   4c) CalculateDeconvKernel
//
//To-do: can this be merged with init()?
void util::SignalShapingServiceMicroBooNE::SetDecon(size_t fftsize, size_t channel)
{
  init();

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::LArFFT> fFFT;
  
  // streamline this method:
  // if the deconvolution kernel is already appropriate for the datasize (aka fftsize) do nothing
  // otherwise, set it to the appropriate size
  // do this test for *every* ss
  // But it will in general only happen once per run! 
  bool setDecon = false; 
  size_t FFTSize = fFFT->FFTSize();
  if (fftsize>FFTSize||fftsize<=FFTSize/2){
    fFFT->ReinitializeFFT( (size_t)fftsize, fFFT->FFTOptions(), fFFT->FFTFitBins() );
    setDecon = true;
  }  
  if(!setDecon) return;
  
  size_t ktype = 1;
 
  for(size_t config=0; config<fNConfigs; ++config) {
    
    //reset fSignalShapingVec for ktype=1
    for (size_t view=0;view<fNViews; ++view) {
      for(unsigned int _wr=0; _wr<fNResponses[ktype][view]; ++_wr) {
        (fSignalShapingVec[config][ktype][view][_wr]).Reset();
      }
    }      

    int mode = 0;
    SetResponseSampling(ktype, config, mode, channel);

    if(config==0) SetFilters();

    // Configure deconvolution kernels.
    for(size_t view=0; view < fNViews; ++view) {
      for(unsigned int _wr=0; _wr<fNResponses[ktype][view]; ++_wr) {
        (fSignalShapingVec[config][ktype][view][_wr]).AddFilterFunction(fFilterVec[view]);
        (fSignalShapingVec[config][ktype][view][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(view));
        (fSignalShapingVec[config][ktype][view][_wr]).CalculateDeconvKernel();
      }
    }
  }//end loop over config
}








//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetFieldResponse(size_t ktype)
{

  art::ServiceHandle<art::TFileService> tfs; 
  
  char buff0[80];

  // Ticks in nanosecond
  // Calculate the normalization of the collection plane
  for(size_t config=0; config<fNConfigs; ++config) {
    double integral = fFieldResponseHistVec[config][ktype][fViewForNormalization][0]->Integral();
    double weight = 1./integral;
    
    // we adjust the size of the fieldresponse vector to account for the stretch
    // and interpolate the histogram to fill the vector with the stretched response        
    for(unsigned int _vw=0; _vw<fNViews; ++_vw) {
      double timeFactor = f3DCorrectionVec[_vw];
      if(!fStretchFullResponse) timeFactor *= fTimeScaleFactor;

      unsigned int NResponses = fNdatadrivenResponses[ktype][_vw];
      if (!fdatadrivenResponse) {
	NResponses = fYZdependentResponse ? fNYZResponses[ktype][_vw] : fNResponses[ktype][_vw];
      }
      
      bool YZflag = true;
      for(unsigned int _wr=0; _wr<NResponses; ++_wr) {
	
	if (YZflag == false) {
	  _wr = 0;
	  _vw = 2;
	}
	
	DoubleVec* responsePtr = &fFieldResponseVec[config][ktype][_vw][_wr];
	TH1F*      histPtr     = fFieldResponseHistVec[config][ktype][_vw][_wr];

	size_t nBins = histPtr->GetNbinsX();
	size_t nResponseBins = nBins*timeFactor;
	responsePtr->resize(nResponseBins);
	double x0 = histPtr->GetBinCenter(1);
	double xf = histPtr->GetBinCenter(nBins);
	double deltaX = (xf - x0)/(nBins-1);

	for(unsigned int _bn=1; _bn<=nResponseBins; ++_bn) {
	  double xVal = x0 + deltaX*(_bn-1)/timeFactor;
	  double yVal = histPtr->Interpolate(xVal);
	  responsePtr->at(_bn-1) = yVal; 
	  responsePtr->at(_bn-1) *= fFieldRespAmpVec[_vw]*weight;
	}

	// fill some histos
	if(_wr==0 && config==0 && ktype==0 && !fHistDone[_vw]) {
	  sprintf(buff0, "hRawResp%i", (int)_vw);
	  fHRawResponse[_vw] = tfs->make<TH1D>(buff0, buff0, nBins, x0-0.5*deltaX, xf+0.5*deltaX);
	  sprintf(buff0, "hStretchedResp%i", (int)_vw);
	  double x0S = timeFactor*x0 - 0.5*deltaX/timeFactor;
	  double xfS = timeFactor*xf + 0.5*deltaX/timeFactor;
	  //std::cout << "title " << buff0 << std::endl;
	  fHStretchedResponse[_vw] = tfs->make<TH1D>(buff0, buff0, nResponseBins, x0S, xfS);
	  for(size_t i=0;i<nBins; ++i) {
	    fHRawResponse[_vw]->SetBinContent(i, histPtr->GetBinContent(i));
	  }
	  for(size_t i=0;i<nResponseBins; ++i) {
	    fHStretchedResponse[_vw]->SetBinContent(i+1, responsePtr->at(i));
	  }
	  fHistDone[_vw] = true;
	}
	
	if (fYZdependentResponse && !fdatadrivenResponse) {
	  if (YZflag == false) {
	    _wr = 2;
	    _vw = 1;
	  }
	  YZflag = false;
	}
      }//end loop over NResponses     
    }//end loop over views
  }//end loop over configs

  return;
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetElectResponse(size_t ktype,double shapingtime, double gain)
{
 
  LOG_DEBUG("SignalShapingMicroBooNE") << "Setting MicroBooNE electronics response function...";

  art::ServiceHandle<util::LArFFT> fft;

  size_t nticks = fft->FFTSize();
  DoubleVec time(nticks,0.);

  fElectResponse.resize(2);
  for(auto& resp : fElectResponse) {
    resp.resize(nticks, 0.);
  }

  //Gain and shaping time variables from fcl file:
  double Ao = 1.0;//Gain
  double To = shapingtime;  //peaking time

  // this is actually sampling time, in ns
  // mf::LogInfo("SignalShapingMicroBooNE") << "Check sampling intervals: "
  //                                  << fSampleRate << " ns"
  //                                  << "Check number of samples: " << fNTicks;

  // The following sets the microboone electronics response function in
  // time-space. Function comes from BNL SPICE simulation of MicroBooNE
  // electronics. SPICE gives the electronics transfer function in
  // frequency-space. The inverse laplace transform of that function
  // (in time-space) was calculated in Mathematica and is what is being
  // used below. Parameters Ao and To are cumulative gain/timing parameters
  // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain.
  // They have been adjusted to make the SPICE simulation to match the
  // actual electronics response. Default params are Ao=1.4, To=0.5us.


  // For the cold electronics,  the gain (i.e. 4.7 mV/fC) represents the peak 
  // height. The shaping time will not affect the peak height, but make the 
  // peak broader  

  for(size_t i=0; i<nticks;++i) {
    time[i] = (1.*i) * fFieldBinWidth[ktype]*1.e-3;
  }
  
  double max = 0;
  for (unsigned int i=0; i < fElectResponse[ktype].size(); ++i) {
    fElectResponse[ktype][i] =
    4.31054*exp(-2.94809*time[i]/To)*Ao - 2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*Ao
    -2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*cos(2.38722*time[i]/To)*Ao
    +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*Ao
    +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*cos(5.18561*time[i]/To)*Ao
    +0.762456*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*Ao
    -0.762456*exp(-2.82833*time[i]/To)*cos(2.38722*time[i]/To)*sin(1.19361*time[i]/To)*Ao
    +0.762456*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao
    -2.6202*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao
    -0.327684*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*Ao +
    +0.327684*exp(-2.40318*time[i]/To)*cos(5.18561*time[i]/To)*sin(2.5928*time[i]/To)*Ao
    -0.327684*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao
    +0.464924*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao;
    
    if (fElectResponse[ktype][i] > max) max = fElectResponse[ktype][i];
  }// end loop over time buckets

  LOG_DEBUG("SignalShapingMicroBooNE") << " Done.";

  // normalize fElectResponse[i], before the convolution
  // Put in overall normalization in a pedantic way:
  // first put in the pulse area per eleectron at the lowest gain setting,
  // then normalize by the actual ASIC gain setting used.
  // This code is executed only during initialization of service,
  // so don't worry about code inefficiencies here.

  //Normalization are the following
  // Peak is firstly normalized to 1
  // thus we expect peak to be 1 * 9390 (fADCPerPCtAtLowestAsicGain) * 1.602e-7 * (1 fC) = 9.39 ADC
  // At 4.7 mV/fC, the ADC value should be 4.7 (mV/fC) * 2 (ADC/mV) ~ 9.4 ADC/fC
  // so the normalization are consistent
  
  double last_integral=0;
  double last_max=0;
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  for(auto& element : fElectResponse[ktype]){
    element /= max;
    element *= fADCPerPCAtLowestASICGain * 1.60217657e-7;
    element *= gain / 4.7;

    if(element > last_max) last_max = element;
    last_integral += element * fFieldBinWidth[ktype] / detprop->SamplingRate();
  }
  return;
}


//----------------------------------------------------------------------
// Calculate microboone filter functions.
void util::SignalShapingServiceMicroBooNE::SetFilters()
{

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<util::LArFFT> fft;

  double ts = detprop->SamplingRate();
  size_t nFFT2 = fft->FFTSize() / 2;

  // Calculate collection filter.
  fFilterVec.resize(fNViews);
  for(auto& filter : fFilterVec) {
    filter.resize(nFFT2+1);
  }

  if(!fGetFilterFromHisto) {
    unsigned int _vw = 0;
    for(auto& func : fFilterTF1Vec) {
      func->SetRange(0, double(nFFT2));
      size_t count = 0;
      
      // now to scale the filter function!
      // only scale params 1,2 &3
      
      double timeFactor = fTimeScaleFactor*f3DCorrectionVec[_vw]*fFilterWidthCorrectionFactor[_vw];
      for(size_t i=1;i<4;++i) {
        func->SetParameter(i, fFilterParamsVec[_vw][i]/timeFactor);
      }
      
      for(unsigned int _bn=0; _bn<=nFFT2; ++_bn) {
        double freq = 500.*_bn/(ts*nFFT2);
        double f = func->Eval(freq);
        if(f!=0.0) count++;
        fFilterVec[_vw][_bn] = TComplex(f, 0.);
      }
      _vw++;
    }
  } 
  else {
    unsigned int _vw = 0;
    for(auto hist : fFilterHistVec) {
      for(unsigned int _bn=1; _bn<=nFFT2+1; ++_bn) {
        double f = hist->GetBinContent(_bn);
        fFilterVec[_vw][_bn-1] = TComplex(f, 0.);
      }
      _vw++;
    }
  }

}


//----------------------------------------------------------------------
// Sample microboone response (the convoluted field and electronic
// response), will probably add the filter later
void util::SignalShapingServiceMicroBooNE::SetResponseSampling(size_t ktype, size_t config, int mode /*=0*/, size_t channel /*=0*/)
{
  // Get services
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* geo = lar::providerFrom<geo::Geometry>();
  art::ServiceHandle<util::LArFFT> fft;
  art::ServiceHandle<art::TFileService> tfs; 
  
  char buff0[80], buff1[80];


  /* This could be a warning, but in principle, there's no reason to restrict the binning

  // Operation permitted only if output of rebinning has a larger bin size
  if( fFieldBinWidth > samplingRate )
  throw cet::exception(__FUNCTION__) << "\033[93m"
  << "Invalid operation: cannot rebin to a more finely binned vector!"
  << "\033[00m" << std::endl;

  */

  size_t view0, view1;
  if(mode==0) {
    view0 = 0;
    view1 = fNViews;
  } else {
    geo::View_t view = geo->View(channel);
    view0 = view;
    view1 = std::min(fNViews,(size_t)view+1);
  }

  size_t nticks = fft->FFTSize();
  DoubleVec SamplingTime( nticks, 0. );
  double deltaInputTime = fFieldBinWidth[ktype];
  for ( size_t itime = 0; itime < nticks; itime++ ) {
    SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
  }
  // Sampling


  // we want to implement new scheme (fStretchFullResponse==false) while retaining the old
  for(size_t view=view0; view<view1; ++view) {
    
    double timeFactor = f3DCorrectionVec[view];
    double plotTimeFactor = 1.0;
    if(fStretchFullResponse) timeFactor *= fTimeScaleFactor;
    else plotTimeFactor = f3DCorrectionVec[view]*fTimeScaleFactor;

    unsigned int NResponses = fNdatadrivenResponses[ktype][view];
    if (!fdatadrivenResponse) {
      NResponses = fYZdependentResponse ? fNYZResponses[ktype][view] : fNResponses[ktype][view];
    }

    bool YZflag = true;
    for(unsigned int _wr=0; _wr<NResponses; ++_wr) {

      if (YZflag == false) {
	_wr = 0;
	view = 2;
      }
	
      const DoubleVec* pResp = &((fSignalShapingVec[config][ktype][view][_wr]).Response_save());
      if(!fHistDoneF[view] &&config==0 && ktype==0 && _wr==0) {
	double xLowF = fFieldLowEdge[ktype]*plotTimeFactor;
	double xHighF = xLowF + 0.001*(fNFieldBins[ktype]+1)*fFieldBinWidth[ktype]*plotTimeFactor;
	double nBins = fNFieldBins[ktype]*plotTimeFactor;
	sprintf(buff0, "FullResponse%i", (int)view);
	fHFullResponse[view] = tfs->make<TH1D>(buff0, buff0, nBins, xLowF, xHighF);                
	for (size_t i=0; i<nBins; ++i) {
	  fHFullResponse[view]->SetBinContent(i+1, pResp->at(i));
	}
      }

      size_t nticks_input = pResp->size();
      DoubleVec InputTime(nticks_input, 0. );
      for (size_t itime = 0; itime < nticks_input; itime++ ) {
	InputTime[itime] = (1.*itime) * deltaInputTime*timeFactor;
      }
      DoubleVec SamplingResp(nticks, 0. );
      size_t SamplingCount = 0;
      size_t startJ = 1;
      SamplingResp[0] = (*pResp)[0];
      for ( size_t itime = 1; itime < nticks; itime++ ) {
	size_t low, high;
	for ( size_t jtime = startJ; jtime < nticks_input; jtime++ ) {
	  if ( InputTime[jtime] >= SamplingTime[itime] ) {
	    low  = jtime - 1;
	    high = jtime;
	    double interpolationFactor = ((*pResp)[high]-(*pResp)[low])/deltaInputTime;
	    SamplingResp[itime] = ((*pResp)[low] + ( SamplingTime[itime] - InputTime[low] ) * interpolationFactor);
	    SamplingResp[itime] /= timeFactor;
	    SamplingCount++;
	    startJ = jtime;
	    break;
	  }
	} // for (  jtime = 0; jtime < nticks; jtime++ )
      } // for (  itime = 0; itime < nticks; itime++ )

      if(!fHistDoneF[view] &&config==0 && ktype==0 && _wr==0) {
	double plotTimeFactor = f3DCorrectionVec[view]*fTimeScaleFactor;
	double xLowF = fFieldLowEdge[ktype]*plotTimeFactor;
	double xHighF = xLowF + 0.001*(fNFieldBins[ktype])*fFieldBinWidth[ktype]*plotTimeFactor;
	double binWidth = 0.5;
	size_t nBins = (xHighF-xLowF+1)/binWidth;
	sprintf(buff1, "SampledResponse%i", (int)view);
	fHSampledResponse[view] = tfs->make<TH1D>(buff1, buff1, nBins, xLowF, xHighF);                
	for (size_t i=0; i<nBins; ++i) {
	  fHSampledResponse[view]->SetBinContent(i+1, SamplingResp[i]);
	}
	fHistDoneF[view] = true;
      }

      if(fPrintResponses) {
	size_t printCount = 0;
	int inc = 1;
	for(size_t i = 0; i<nticks; i+=inc) {
	  if((printCount+1)%10==0) std::cout << std::endl;
	  printCount++;
	  if (printCount>=100) {inc = 100;}
	}
      }

      (fSignalShapingVec[config][ktype][view][_wr]).AddResponseFunction( SamplingResp, true);

      if (fYZdependentResponse && !fdatadrivenResponse) {
        if(YZflag == false){ 
	  _wr = 2; 
	  view = 1;
	}
	YZflag = false;
      }

    }  //  loop over wires
  } // loop over views
  
  return;
}


//-----Give Gain Settings to SimWire-----//jyoti
double util::SignalShapingServiceMicroBooNE::GetASICGain(unsigned int const channel) const
{
  init();

  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);
  
  size_t config = GetConfig(channel);
  
  double gain = 0;
  switch(view){
    case geo::kU:
      gain = fASICGainInMVPerFC[config].at(0);
      break;
    case geo::kV:
      gain = fASICGainInMVPerFC[config].at(1);
      break;
    case geo::kZ:
      gain = fASICGainInMVPerFC[config].at(2);
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  return gain;
}


//-----Give Shaping time to SimWire-----//jyoti
double util::SignalShapingServiceMicroBooNE::GetShapingTime(unsigned int const channel) const
{
  init();

  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  size_t config = GetConfig(channel);

  double shaping_time = 0;
  switch(view){
    case geo::kU:
      shaping_time = fShapeTimeConst[config].at(0);
      break;
    case geo::kV:
      shaping_time = fShapeTimeConst[config].at(1);
      break;
    case geo::kZ:
      shaping_time = fShapeTimeConst[config].at(2);
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  return shaping_time;
}


double util::SignalShapingServiceMicroBooNE::GetRawNoise(unsigned int const channel) const
{
  init();

  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  size_t config = GetConfig(channel);
  
  unsigned int plane;
  switch(view){
    case geo::kU:
      plane = 0;
      break;
    case geo::kV:
      plane = 1;
      break;
    case geo::kZ:
      plane = 2;
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }

  double shapingtime = fShapeTimeConst[config].at(plane);
  double gain = fASICGainInMVPerFC[config].at(plane);
  int temp;
  if (std::abs(shapingtime - 0.5)<1e-6){
    temp = 0;
  }else if (std::abs(shapingtime - 1.0)<1e-6){
    temp = 1;
  }else if (std::abs(shapingtime - 2.0)<1e-6){
    temp = 2;
  }else{
    temp = 3;
  }

  auto tempNoise = fNoiseFactVec.at(plane);
  double rawNoise = tempNoise.at(temp);

  rawNoise *= gain/4.7;
  return rawNoise;
}


double util::SignalShapingServiceMicroBooNE::GetDeconNoise(unsigned int const channel) const
{ 

  double deconNoise = this->GetRawNoise(channel);
  deconNoise *= 4.7/this->GetASICGain(channel);
  deconNoise = deconNoise /4096.*2000./4.7 *6.241*1000/fDeconNorm;
  
  return deconNoise;
}


int util::SignalShapingServiceMicroBooNE::FieldResponseTOffset(unsigned int const channel, size_t ktype=0) const
{
  //no init needed - fFieldResponseTOffset is initialized in reconfigure()
  
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  double time_offset = 0;
  switch(view){
    case geo::kU:
      time_offset = fFieldResponseTOffset[ktype].at(0);
      break;
    case geo::kV:
      time_offset = fFieldResponseTOffset[ktype].at(1);
      break;
    case geo::kZ: 
      time_offset = fFieldResponseTOffset[ktype].at(2); 
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  auto tpc_clock = lar::providerFrom<detinfo::DetectorClocksService>()->TPCClock();
  return tpc_clock.Ticks(time_offset/1.e3);
}


//----------------------------------------------------------------------
// Return the correct configuration for this channel
size_t util::SignalShapingServiceMicroBooNE::GetConfig(size_t channel) const
{
  if(fNConfigs<=1 || fConfigMap.size()==0) return 0;

  // for a test with the special sim event, set n to 1000

  if(fConfigMap.find(channel)==fConfigMap.end()) return 0;

  return fConfigMap[channel];
}


//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceMicroBooNE::SignalShaping(size_t channel, size_t wire, size_t ktype) const
{
  init();

  // Figure out plane type.
  art::ServiceHandle<geo::Geometry> geom;
  auto view = (size_t)geom->View(channel);
  if(view<0||view>fNViews-1) {
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine view\n";
  }

  //modify view/wire if needed
  if(fYZdependentResponse == true && view == 1 && wire == 1 && fdatadrivenResponse == false){
    view = 2;
    wire = 0;
  }

  size_t config = GetConfig(channel);

  // Return appropriate shaper.
  return fSignalShapingVec[config][ktype][view][wire];
}


//----------------------------------------------------------------------
// Get convolution kernel from SignalShaping service for use in CalWire's
// DeconvoluteInducedCharge() - added by M. Mooney
const std::vector<TComplex>& util::SignalShapingServiceMicroBooNE::GetConvKernel(unsigned int channel, unsigned int wire) const
{
  init();

  art::ServiceHandle<geo::Geometry> geom;
  auto view = (size_t)geom->View(channel);
  if(view<0||view>fNViews-1) {
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
    << " View\n";
  }

  size_t config = GetConfig(channel);
  
  // Return appropriate shaper.
  return fSignalShapingVec[config][0][view][wire].ConvKernel();
}


//----------------------------------------------------------------------
// Evaluate 2D filter used in induced charge deconvolution (M. Mooney)
double util::SignalShapingServiceMicroBooNE::Get2DFilterVal(size_t planeNum, size_t freqDimension, double binFrac) const
{
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double ts = detprop->SamplingRate();

  double freq;
  double filtVal;
  double val;
  if(freqDimension == 1) {
    freq = fFilterScaleVecICTime.at(planeNum)*(500.0/ts)*2.0*(0.5-fabs(binFrac-0.5));
    filtVal = fFilterTF1VecICTime.at(planeNum)->Eval(freq);

    if(freq < fFilterICTimeMaxFreq.at(planeNum))
      val = fFilterICTimeMaxVal.at(planeNum);
    else
      val = filtVal;

    return val;
  }
  else if(freqDimension == 2) {
    freq = fFilterScaleVecICWire.at(planeNum)*(0.5-fabs(binFrac-0.5));
    filtVal = fFilterTF1VecICWire.at(planeNum)->Eval(freq);

    if(freq < fFilterICWireMaxFreq.at(planeNum))
      val = fFilterICWireMaxVal.at(planeNum);
    else
      val = filtVal;

    return val;
  }
  else {
    return 0.0;
  }
}


//----------------------------------------------------------------------
// Return 2D filter normalization, used in induced charge deconvolution (M. Mooney)
double util::SignalShapingServiceMicroBooNE::Get2DFilterNorm(size_t planeNum) const
{
  return fFilterNormVecIC.at(planeNum);
}

namespace util {
  
  DEFINE_ART_SERVICE(SignalShapingServiceMicroBooNE)
  
}
