////////////////////////////////////////////////////////////////////////
// $Id: SimWireMicroBooNE.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireMicroBooNE class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard library
#include <stdexcept> // std::range_error
#include <vector>
#include <string>
#include <algorithm> // std::fill()
#include <functional>
#include <random>
#include <chrono>

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// ROOT libraries
#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

// art library and utilities
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this is not portable
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardataobj/Simulation/sim.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires.
  class SimWireMicroBooNE : public art::EDProducer {

  public:

    explicit SimWireMicroBooNE(fhicl::ParameterSet const& pset);
    virtual ~SimWireMicroBooNE();

    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void GenNoiseInTime(std::vector<float> &noise, double noise_factor) const;
    void GenNoiseInFreq(std::vector<float> &noise, double noise_factor) const;
    void GenNoisePostFilter(std::vector<float> &noise, double noise_factor, size_t view, int chan);
    void MakeADCVec(std::vector<short>& adc, std::vector<float> const& noise, 
                    std::vector<double> const& charge, float ped_mean) const;

    std::string             fDriftEModuleLabel; ///< module making the ionization electrons
    raw::Compress_t         fCompression;       ///< compression type to use

    double                  fNoiseWidth;        ///< exponential noise width (kHz)
    double                  fNoiseRand;         ///< fraction of random "wiggle" in noise in freq. spectrum
    double                  fLowCutoff;         ///< low frequency filter cutoff (kHz)
    double                  fSampleRate;        ///< sampling rate in ns

    size_t                  fNTicks;	        ///< number of ticks of the clock    
    unsigned int            fNTimeSamples;      ///< number of ADC readout samples in all readout frames (per event)	 	   

    std::vector<TH1D*>      fNoiseDist;     ///< distribution of noise counts, one per plane
    bool                    fGetNoiseFromHisto; ///< if True -> Noise from Histogram of Freq. spectrum
    unsigned short          fGenNoise;          ///< 0 -> no noise, 1: time domain, 2: freq domain, 3: postfilter
    std::string             fNoiseFileFname;
    std::string             fNoiseHistoName;
    TH1D*                   fNoiseHist;         ///< distribution of noise counts

    std::map< double, int > fShapingTimeOrder;
    std::string             fTrigModName;       ///< Trigger data product producer name
    
    bool                    fSimDeadChannels;   ///< if True, simulate dead channels using the ChannelStatus service.  If false, do not simulate dead channels

    bool fMakeNoiseDists;

    bool        fTest; // for forcing a test case
    std::vector<sim::SimChannel> fTestSimChannel_v;
    size_t      fTestWire;
    std::vector<size_t> fTestIndex;
    std::vector<double> fTestCharge;

    int         fSample; // for histograms, -1 means no histos

    //std::vector<std::vector<std::vector<int> > > fYZwireOverlap; //channel ranges for shorted wires and corresponding channel ranges for wires effected on other planes

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;

    ::detinfo::ElecClock fClock; ///< TPC electronics clock

    //
    // Needed for post-filter noise (pfn) generator
    //
    TF1  *_pfn_f1;
    TF1  *_pfn_MyPoisson;
    TVirtualFFT *_pfn_ifft;
    std::vector<double> _pfn_rho_v;
    std::vector<double> _pfn_value_re;
    std::vector<double> _pfn_value_im;
    float gammaRand;

  }; // class SimWireMicroBooNE

  /*  namespace {
    size_t _ch = 0;
    size_t _wr = 0;
  }
  */
  DEFINE_ART_MODULE(SimWireMicroBooNE)

  //-------------------------------------------------
  SimWireMicroBooNE::SimWireMicroBooNE(fhicl::ParameterSet const& pset)
  : fNoiseHist(0)
    , _pfn_f1(nullptr)
    , _pfn_MyPoisson(nullptr)
    , _pfn_ifft(nullptr)
    , _pfn_rho_v()
    , _pfn_value_re()
    , _pfn_value_im()
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;

    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed" and "SeedPedestal"
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    Seeds->createEngine(*this, "HepJamesRandom", "noise", pset, "Seed");
    Seeds->createEngine(*this, "HepJamesRandom", "pedestal", pset, "SeedPedestal");

  }

  //-------------------------------------------------
  SimWireMicroBooNE::~SimWireMicroBooNE()
  {
    delete fNoiseHist;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::reconfigure(fhicl::ParameterSet const& p)
  {
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fNoiseWidth       = p.get< double              >("NoiseWidth");
    fNoiseRand        = p.get< double              >("NoiseRand");
    fLowCutoff        = p.get< double              >("LowCutoff");
    fGetNoiseFromHisto= p.get< bool                >("GetNoiseFromHisto");
    fGenNoise         = p.get< unsigned short      >("GenNoise");
    fSimDeadChannels  = p.get< bool                >("SimDeadChannels");

    fMakeNoiseDists   = p.get< bool                >("MakeNoiseDists", false);

    fTrigModName      = p.get< std::string         >("TrigModName");
    fTest             = p.get<bool                 >("Test");
    fTestWire         = p.get< size_t              >("TestWire");
    fTestIndex        = p.get< std::vector<size_t> >("TestIndex");
    fTestCharge       = p.get< std::vector<double> >("TestCharge");
    if(fTestIndex.size() != fTestCharge.size())
      throw cet::exception(__FUNCTION__)<<"# test pulse mismatched: check TestIndex and TestCharge fcl parameters...";
    fSample           = p.get<int                  >("Sample");

    //fYZwireOverlap    = p.get<std::vector<std::vector<std::vector<int> > > >("YZwireOverlap");

    //Map the Shaping Times to the entry position for the noise ADC
    //level in fNoiseFactInd and fNoiseFactColl
    fShapingTimeOrder = { {0.5, 0}, {1.0, 1}, {2.0, 2}, {3.0, 3} };

    if(fGetNoiseFromHisto)
    {
      fNoiseHistoName= p.get< std::string         >("NoiseHistoName");

      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(p.get<std::string>("NoiseFileFname"), fNoiseFileFname);

      /*TFile * in=new TFile(fNoiseFileFname.c_str(),"READ");
      TH1D * temp=(TH1D *)in->Get(fNoiseHistoName.c_str());

      if(temp!=NULL)
      {
        fNoiseHist=new TH1D(fNoiseHistoName.c_str(),fNoiseHistoName.c_str(),temp->GetNbinsX(),0,temp->GetNbinsX());
        temp->Copy(*fNoiseHist);
      }
      else
        throw cet::exception("SimWireMicroBooNE") << "Could not find noise histogram in Root file\n";
      in->Close();
*/
    }
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSampleRate    = detprop->SamplingRate();
    fNTimeSamples  = detprop->NumberTimeSamples();

    // make the histos if not already made
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    return;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::beginJob()
  {

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    char buff0[80], buff1[80];

    if(fMakeNoiseDists) {
      fNoiseDist.resize(3,0);
      for(int view = 0; view<3; ++view) {
        sprintf(buff0, "Noise%i", view);
        sprintf(buff1, ";Noise on Plane %i(ADC);", view);
        fNoiseDist[view]  = tfs->make<TH1D>(buff0, buff1, 1000,   -30., 30.);
      }
    }

    if(fTest){
      art::ServiceHandle<geo::Geometry> geo;  
      if(geo->Nchannels()<=fTestWire)
        throw cet::exception(__FUNCTION__)<<"Invalid test wire channel: "<<fTestWire;

      std::vector<unsigned int> channels;
      for(auto const& plane_id : geo->IteratePlaneIDs())
        channels.push_back(geo->PlaneWireToChannel(plane_id.Plane,fTestWire));

      double xyz[3] = { std::numeric_limits<double>::max() };
      for(auto const& ch : channels) {

        fTestSimChannel_v.push_back(sim::SimChannel(ch));

        for(size_t i=0; i<fTestIndex.size(); ++i){

          fTestSimChannel_v.back().AddIonizationElectrons(-1,
                                                          fTestIndex[i],
                                                          fTestCharge[i],
                                                          xyz,
                                                          std::numeric_limits<double>::max());
        }
      }
    }

    return;

  }

  //-------------------------------------------------
  void SimWireMicroBooNE::endJob()
  {}

  void SimWireMicroBooNE::produce(art::Event& evt)
  {

    //--------------------------------------------------------------------
    //
    // Get all of the services we will be using
    //
    //--------------------------------------------------------------------
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg 
       = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
    //electronics conditions
    const lariov::ElectronicsCalibProvider& elec_provider 
       = art::ServiceHandle<lariov::ElectronicsCalibService>()->GetProvider();
    
    //get rng for pedestals
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("pedestal");   
    
    //channel status for simulating dead channels
    const lariov::ChannelStatusProvider& ChannelStatusProvider
       = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

    //get the FFT
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->ReinitializeFFT(fNTimeSamples,fFFT->FFTOptions(),fFFT->FFTFitBins());
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 )
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize " << fNTicks << " not a power of 2. "
      << "May cause issues in (de)convolution.\n";

    if ( fNTimeSamples > fNTicks )
      mf::LogError("SimWireMicroBooNE") << "Cannot have number of readout samples "
      << fNTimeSamples << " greater than FFTSize " << fNTicks << "!";

    // TFileService
    art::ServiceHandle<art::TFileService> tfs;
    
    //TimeService
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard> tss;
    // In case trigger simulation is run in the same job...
    tss->preProcessEvent(evt);
    auto const* ts = tss->provider();
    
    // Check if trigger data product exists or not. If not, throw a warning
    art::Handle< std::vector<raw::Trigger> > trig_array;
    evt.getByLabel(fTrigModName, trig_array);
    if(!trig_array.isValid())

      std::cout << std::endl << "  "
      << "\033[95m" << "<<" << __PRETTY_FUNCTION__ << ">>" << "\033[00m"
      << std::endl << "  "
      << "\033[93m"
      << " No trigger data exists => will use the default trigger time set in TimeService..."
      << "\033[00m"
      << std::endl;

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    const size_t N_CHANNELS = geo->Nchannels();
    const size_t N_VIEWS = geo->Nplanes();

    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;

   
    //--------------------------------------------------------------------
    //
    // Get the SimChannels, which we will use to produce RawDigits
    //
    //--------------------------------------------------------------------

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(N_CHANNELS,nullptr);
    if(!fTest){
      std::vector<const sim::SimChannel*> chanHandle;
      evt.getView(fDriftEModuleLabel,chanHandle);
      for(size_t c = 0; c < chanHandle.size(); ++c){
        channels.at(chanHandle.at(c)->Channel()) = chanHandle.at(c);
      }
    }else{
      for(size_t c = 0; c<fTestSimChannel_v.size(); ++c)
        channels.at(fTestSimChannel_v[c].Channel()) = &(fTestSimChannel_v[c]);
    }

    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(N_CHANNELS);

    //first vector index is channel
    std::vector< std::vector<std::unique_ptr<util::ResponseParams> > > responseParamsVec;
    responseParamsVec.resize(N_CHANNELS);
    
   
    //--------------------------------------------------------------------
    //
    // Store energy deposits in ResponseParams objects
    //
    //--------------------------------------------------------------------
    std::vector<int> first_channel_in_view(N_VIEWS,-1);
    for(unsigned int chan = 0; chan < N_CHANNELS; ++chan) {	
      size_t view = (size_t)geo->View(chan);

      if (first_channel_in_view[view] == -1) {
	first_channel_in_view[view] = chan;
      }

      const sim::SimChannel* sc = channels.at(chan);
      if( !sc ) continue;

      auto const& timeSlices = sc->TDCIDEMap();

      // remove the time offset
      int time_offset = 0;//sss->FieldResponseTOffset(chan);

      for(auto timeSlice : timeSlices) {
	auto tdc = timeSlice.first;
	if( tdc < 0 ) continue;
	auto t = ts->TPCTDC2Tick(tdc)+1; // +1 added because nominal detsim rounds up (B. Russell)

	int raw_digit_index = (int)( (t + time_offset) >= 0 ? t+time_offset : (fNTicks + (t+time_offset)) );
	if(raw_digit_index <= 0 || raw_digit_index >= (int)fNTicks) continue;

	auto const& energyDeposits = timeSlice.second;
	for(auto energyDeposit : energyDeposits) {
	  double charge = (double)energyDeposit.numElectrons;
	  double y = (double)energyDeposit.y;
	  double z = (double)energyDeposit.z;
	  if(charge == 0) continue;
	  responseParamsVec[chan].emplace_back(new util::ResponseParams(charge, y, z, raw_digit_index));
	}
      }
    } // channels

    
    //--------------------------------------------------------------------
    //
    // Loop over channels, generating pedestal and noise
    // Then loop over channel's energy deposits and convolute appropriate response
    //
    //--------------------------------------------------------------------
       
    // vectors for working in the following for loop
    std::vector<short>    adcvec(fNTimeSamples, 0);
    std::vector<double>   chargeWork(fNTicks,0.);
    std::vector<double>   tempWork(fNTicks,0.);
    std::vector<float>    noisetmp(fNTicks,0.);
    
    int step = 0;
    for(unsigned int chan = 0; chan < N_CHANNELS; chan++) {

      //clean up working vectors from previous iteration of loop
      adcvec.resize(fNTimeSamples); //compression may have changed the size of this vector
      noisetmp.resize(fNTicks); //just in case
      std::fill(adcvec.begin(),     adcvec.end(),     0);
      std::fill(chargeWork.begin(), chargeWork.end(), 0.);
      std::fill(tempWork.begin(),   tempWork.end(),   0.);
      std::fill(noisetmp.begin(),   noisetmp.end(),   0.);
      
      // make sure chargeWork is correct size
      if (chargeWork.size() < fNTimeSamples)
        throw std::range_error("SimWireMicroBooNE: chargeWork vector too small");
	
      //use channel number to set some useful numbers
      size_t view = (size_t)geo->View(chan);
           
      //Get pedestal with random gaussian variation
      CLHEP::RandGaussQ rGaussPed(engine, 0.0, pedestalRetrievalAlg.PedRms(chan));
      float ped_mean = pedestalRetrievalAlg.PedMean(chan) + rGaussPed.fire();
        
      //Generate Noise
      if (fGenNoise) {
	double noise_factor = 0.0;
	auto tempNoiseVec = sss->GetNoiseFactVec();
	double shapingTime = elec_provider.ShapingTime(chan);
	double asicGain    = elec_provider.Gain(chan);
	
	double st_diff = 99999999.9;
	for (auto iST = fShapingTimeOrder.begin(); iST != fShapingTimeOrder.end(); ++iST) {
	  if ( fabs(iST->first - shapingTime) < st_diff ) {
	    st_diff = fabs(iST->first - shapingTime);
	    noise_factor = tempNoiseVec[view][iST->second];
	  }
	}
	noise_factor *= asicGain/4.7;

        if (fGenNoise==1)
          GenNoiseInTime(noisetmp, noise_factor);
        else if(fGenNoise==2)
          GenNoiseInFreq(noisetmp, noise_factor);
	else if(fGenNoise==3)
	  GenNoisePostFilter(noisetmp, noise_factor, view, chan);
	    
   
	//Add Noise to NoiseDist Histogram
	if(fMakeNoiseDists) {
          for (size_t i=step; i < fNTimeSamples; i+=1000) {
            fNoiseDist[view]->Fill(noisetmp[i]);
          }
	  ++step;
	}	
      }//end Generate Noise


      //If the channel is bad, we can stop here
      //if you are using the UbooneChannelStatusService, then this removes disconnected, "dead", and "low noise" channels
      if (fSimDeadChannels && (ChannelStatusProvider.IsBad(chan) || !ChannelStatusProvider.IsPresent(chan)) ) {
        MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
        raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
        rd.SetPedestal(ped_mean);
        digcol->push_back(std::move(rd));
        continue; //on to next channel
      }
      
      
      //Channel is good, so convolute response onto all charges and fill the chargeWork vector
      sss->Convolute(chan, tempWork, responseParamsVec[chan]);
      for(size_t bin = 0; bin < fNTicks; ++bin) {
        chargeWork[bin] += tempWork[bin];
      }
      
      
      /*for (auto& item : responseParamsVec[chan]) {
        std::fill(tempWork.begin(), tempWork.end(), 0.);
	auto charge = item->getCharge();
        if(charge==0) continue;
	auto raw_digit_index = item->getTime();
	if(raw_digit_index > 0 && raw_digit_index < fNTicks) {
          tempWork.at(raw_digit_index) += charge;
	}
	
	sss->Convolute(chan, tempWork, item->getY(), item->getZ());
        
        // now add the result into the "charge" vector
        for(size_t bin = 0; bin < fNTicks; ++bin) {
          chargeWork[bin] += tempWork[bin];
        }
      }*/


      // add this digit to the collection;
      // adcvec is copied, not moved: in case of compression, adcvec will show
      // less data: e.g. if the uncompressed adcvec has 9600 items, after
      // compression it will have maybe 5000, but the memory of the other 4600
      // is still there, although unused; a copy of adcvec will instead have
      // only 5000 items. All 9600 items of adcvec will be recovered for free
      // and used on the next loop.
      MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
      raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
      rd.SetPedestal(ped_mean);
      digcol->push_back(std::move(rd)); // we do move the raw digit copy, though
      
    }// end of loop over channels

 
    evt.put(std::move(digcol));
    return;
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::MakeADCVec(std::vector<short>& adcvec, std::vector<float> const& noisevec, 
                                     std::vector<double> const& chargevec, float ped_mean) const {


    for(unsigned int i = 0; i < fNTimeSamples; ++i) {

       float adcval = noisevec[i] + chargevec[i] + ped_mean;

      //allow for ADC saturation
      if ( adcval > adcsaturation )
	adcval = adcsaturation;
      //don't allow for "negative" saturation
      if ( adcval < 0 )
	adcval = 0;

      adcvec[i] = (unsigned short)TMath::Nint(adcval);

    }// end loop over signal size

    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInTime(std::vector<float> &noise, double noise_factor) const
  {
    //ART random number service
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
    CLHEP::RandGaussQ rGauss(engine, 0.0, noise_factor);

    //In this case noise_factor is a value in ADC counts
    //It is going to be the Noise RMS
    //loop over all bins in "noise" vector
    //and insert random noise value
    for (unsigned int i=0; i<noise.size(); i++)
      noise.at(i) = rGauss.fire();
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInFreq(std::vector<float> &noise, double noise_factor) const
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
    CLHEP::RandFlat flat(engine,-1,1);

    if(noise.size() != fNTicks)
      throw cet::exception("SimWireMicroBooNE")
      << "\033[93m"
      << "Frequency noise vector length must match fNTicks (FFT size)"
      << " ... " << noise.size() << " != " << fNTicks
      << "\033[00m"
      << std::endl;

    // noise in frequency space
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);

    double pval = 0.;
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(size_t i=0; i< fNTicks/2+1; ++i){
      // exponential noise spectrum
      flat.fireArray(2,rnd,0,1);
      //if not from histo or in time --> then hardcoded freq. spectrum
      if( !fGetNoiseFromHisto )
      {
        pval = noise_factor*exp(-(double)i*binWidth/fNoiseWidth);
        // low frequency cutoff
        lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
        // randomize 10%
        
        pval *= lofilter*((1-fNoiseRand)+2*fNoiseRand*rnd[0]);
      }
      
      
      else
      {
        // histogram starts in bin 1!
        pval = fNoiseHist->GetBinContent(i+1)*((1-fNoiseRand)+2*fNoiseRand*rnd[0])*noise_factor;
        //mf::LogInfo("SimWireMicroBooNE")  << " pval: " << pval;
      }
      phase = rnd[1]*2.*TMath::Pi();
      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency.at(i) += tc;
    }
    
    
    // mf::LogInfo("SimWireMicroBooNE") << "filled noise freq";
    
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();
    
    // multiply each noise value by fNTicks as the InvFFT
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    //Also need to scale so that noise RMS matches that asked
    //in fhicl parameter (somewhat arbitrary scaling otherwise)
    //harcode this scaling factor (~20) for now
    for(unsigned int i = 0; i < noise.size(); ++i) noise.at(i) *= 1.*(fNTicks/20.);
    
  }
  
  //---------------------------------------------------------

  void SimWireMicroBooNE::GenNoisePostFilter(std::vector<float> &noise, double noise_factor, size_t view, int chan)
  {
    
     //electronics conditions
    const lariov::ElectronicsCalibProvider& elec_provider 
       = art::ServiceHandle<lariov::ElectronicsCalibService>()->GetProvider();
    
    // noise is a vector of size fNTicks, which is the number of ticks
    const size_t waveform_size = noise.size();
       

    Double_t ShapingTime = elec_provider.ShapingTime(chan);

      if(!_pfn_f1) _pfn_f1 = new TF1("_pfn_f1", "([0]*exp(-0.5*(((x*9592/2)-[1])/[2])**2)*exp(-0.5*pow(x*9592/(2*[3]),[4]))+[5])", 0.0, (double)waveform_size/2);

    if(_pfn_rho_v.empty()) _pfn_rho_v.resize(waveform_size);
    if(_pfn_value_re.empty()) _pfn_value_re.resize(waveform_size);
    if(_pfn_value_im.empty()) _pfn_value_im.resize(waveform_size);

    //**Setting lambda/
    Double_t params[1] = {0.};
    Double_t fitpar[6] = {0.};

    if(ShapingTime > 1.5 && ShapingTime <= 2.5) {
      params[0] = 3.3708; //2us

      // wiener-like
      fitpar[0] = 8.49571e+02;
      fitpar[1] = 6.60496e+02;
      fitpar[2] = 5.68387e+02;
      fitpar[3] = 1.02403e+00;
      fitpar[4] = 1.57143e-01;
      fitpar[5] = 4.79649e+01;


    }
    else if(ShapingTime > 0.75 && ShapingTime <= 1.5) {
      params[0] = 3.5125; //1us
      fitpar[0] = 14.4;
      fitpar[1] = 35.1;
      fitpar[2] = 0.049;
      fitpar[3] = 6.0e-9;
      fitpar[4] = 2.4;
    }else
      throw cet::exception("SimWireMicroBooNE") << "<<" << __FUNCTION__ << ">> not supported shaping time " << ShapingTime << std::endl;
    
    _pfn_f1->SetParameters(fitpar);

    Int_t n = waveform_size;
    
    TVirtualFFT::SetTransform(0);
    
    // seed gamma-distibuted random number with mean params[0]
    // replacing continuous Poisson distribution from ROOT
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::gamma_distribution<double> distribution(params[0]);
    
    // For every tick...
    for(size_t i=0; i<waveform_size; i++){
      
      Double_t freq;
      if (i < n/2.){
        freq = (i)*2./n; //2 MHz digitization
      }else{
        freq = (n-i)*2./n;
      }

      // Draw gamma-dstributed random number
      gammaRand = distribution(generator); 
      
      // Define FFT parameters
      _pfn_rho_v[i] = _pfn_f1->Eval(freq) * gammaRand/params[0];
      Double_t rho = _pfn_rho_v[i];
      Double_t phi = gRandom->Uniform(0,1) * 2. * TMath::Pi();
      
      _pfn_value_re[i] = rho*cos(phi)/((double)waveform_size);
      _pfn_value_im[i] = rho*sin(phi)/((double)waveform_size);
    }
    
    // Inverse FFT
    if(!_pfn_ifft) _pfn_ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    _pfn_ifft->SetPointsComplex(&_pfn_value_re[0],&_pfn_value_im[0]);
    _pfn_ifft->Transform();
    
    // Produce fit histogram from the FFT fo the real part
    TH1 *fb = 0;
    fb = TH1::TransformHisto(_pfn_ifft,fb,"Re");
    
    // Get wire length 
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    std::vector<geo::WireID> wireIDs = geom->ChannelToWire(chan);
    geo::WireGeo const& wire = geom->Wire(wireIDs.front());
    double wirelength = wire.HalfL() * 2;
  
    // Calculate RMS -----------------------------------------------------
    // Calculating using the 16th, 50th, and 84th percentiles.
    // Because the signal is expected to be above the 84th percentile, this 
    // effectively vetos the signal.
    
    Double_t min = fb->GetMinimum();
    Double_t max = fb->GetMaximum();	 
    TH1F* h_rms = new TH1F("h_rms", "h_rms", Int_t(10*(max-min+1)), min, max+1);
    
    
    for(size_t i=0; i < waveform_size; ++i){
      h_rms->Fill(fb->GetBinContent(i+1));
    }
    
    double par[3];
    double rms_quantilemethod = 0.0;
    if (h_rms->GetSum()>0){
      double xq = 0.5-0.34;
      h_rms->GetQuantiles(1, &par[0], &xq);
      
      xq = 0.5;
      h_rms->GetQuantiles(1, &par[1], &xq);
      
      xq = 0.5+0.34;
      h_rms->GetQuantiles(1, &par[2], &xq);
      
      rms_quantilemethod = sqrt((pow(par[1]-par[0],2)+pow(par[2]-par[1],2))/2.);
    
    }
    
    // Scaling noise RMS with wire length dependance
    double baseline = 1.17764;

    double para = 0.4616;
    double parb = 0.19;
    double parc = 1.07;

    // 0.77314 scale factor accounts for fact that original DDN designed based
    // on the Y plane, updated fit takes average of wires on 2400 on each plane
    double scalefactor = 0.83 * (rms_quantilemethod/baseline) * sqrt(para*para + pow(parb*wirelength/100 + parc, 2));
    for(size_t i=0; i<waveform_size; ++i) {
      noise[i] = fb->GetBinContent(i+1)*scalefactor;
    }
    /*
      double average=0;
      for(auto const& v : noise) average += v;
      average /= ((double)waveform_size);
      std::cout<<"\033[93m Average ADC: \033[00m" << average << std::endl;
    */
    delete fb;
  }
  
}
