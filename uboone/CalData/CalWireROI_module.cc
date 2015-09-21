////////////////////////////////////////////////////////////////////////
//
// CalWireROI class - variant of CalWire that deconvolves in 
// Regions Of Interest
//
// baller@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>

// ROOT libraries
#include "TComplex.h"
#include "TH1D.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "RecoBaseArt/WireCreator.h"
#include "Utilities/LArFFT.h"
#include "Utilities/AssociationUtil.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "WaveformPropertiesAlg.h"

/* unused function
namespace {

	void DumpWire(const recob::Wire& wire) {
		
		const size_t pagesize = 10;
		mf::LogDebug log("DumpWire");
		
		const recob::Wire::RegionsOfInterest_t& wireSignalRoI = wire.SignalROI();
		
		log << "\nDumpWire: wire on view " << ((int) wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignalRoI.n_ranges() << " regions of interest:";
		size_t iRoI = 0;
		auto RoI = wireSignalRoI.begin_range(), rend = wireSignalRoI.end_range();
		while (RoI != rend) {
			++iRoI;
			log << "\nDumpWire: [RoI " << iRoI << "] starts at " << RoI->begin_index()
				<< ", " << RoI->size() << " samples:";
			size_t iSample = 0;
			for (auto sample: *RoI) {
				if (iSample % pagesize == 0)
					log << "\nDumpWire: [RoI " << iRoI << "/" << iSample << "]";
				log << '\t' << sample;
				++iSample;
			} // for sample
			++RoI;
		} // for RoI
		
		
		const recob::Wire::RegionsOfInterest_t& wireSignal = wireSignalRoI;
		std::vector<float> buffer(pagesize), prev_buffer;
		log << "\nDumpWire: wire on view " << ((int)wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignal.size() << " samples:";
		size_t i = 0, nSame = 0;
		auto iSample = wireSignal.begin(), send = wireSignal.end();
		while (iSample < send) {
			i += prev_buffer.size();
			buffer.assign(iSample, std::min(iSample + pagesize, send));
			iSample += buffer.size();
			if (buffer == prev_buffer) {
				++nSame;
			}
			else {
				if (nSame > 0) {
					log << "\nDumpWire: [" << i << "]  ... and " << nSame << " more";
					nSame = 0;
				}
				
				if (i % pagesize == 0) log << "\nDumpWire: [" << i << "]";
				for (auto value: buffer) log << '\t' << value;
				buffer.swap(prev_buffer);
			}
		} // while
		if (nSame > 0) {
			log << "\nDumpWire: [" << i << "]  ... and " << nSame << " more to the end";
			nSame = 0;
		}
	} // DumpWire()
}
*/

///creation of calibrated signals on wires
namespace caldata {

  class CalWireROI : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireROI(fhicl::ParameterSet const& pset); 
    virtual ~CalWireROI();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
    void reconfFFT(int temp_fftsize);
  private:
    
    std::string  fDigitModuleLabel;  ///< module that made digits

    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data
    unsigned short fNoiseSource;     ///< 0 = none, 1 = digitVec->GetSigma(), 2 = sss->GetRawNoise(channel)
    unsigned short fNumBinsHalf;    ///< Make running average using 2 * fNumBinsHalf + 1 ticks
    std::vector<float> fThreshold;   ///< abs(threshold) ADC counts for ROI
    int fFFTSize;        ///< FFT size for ROI deconvolution
    std::vector<unsigned short> fPreROIPad; ///< ROI padding
    std::vector<unsigned short> fPostROIPad; ///< ROI padding
    bool fDoBaselineSub;  ///< Do baseline subtraction after deconvolution?
    int  fSaveWireWF;     ///< Save recob::wire object waveforms
    size_t fEventCount;  ///< count of event processed
    int  fMaxAllowedChanStatus;
    
    void doDecon(std::vector<float>& holder, 
      raw::ChannelID_t channel, unsigned int thePlane,
      std::vector<std::pair<unsigned int, unsigned int>> rois,
      std::vector<std::pair<unsigned int, unsigned int>> holderInfo,
      recob::Wire::RegionsOfInterest_t& ROIVec,
      art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss);
    float SubtractBaseline(std::vector<float>& holder, float basePre,
			   float basePost,unsigned int roiStart,unsigned int roiLen,
			   unsigned int dataSize);

    bool fDoBaselineSub_WaveformPropertiesAlg;
    util::WaveformPropertiesAlg<float> fROIPropertiesAlg;
    float SubtractBaseline(const std::vector<float>& holder);
    
  protected: 
    
  }; // class CalWireROI

  DEFINE_ART_MODULE(CalWireROI)
  
  //-------------------------------------------------
  CalWireROI::CalWireROI(fhicl::ParameterSet const& pset):
    fROIPropertiesAlg(pset.get<fhicl::ParameterSet>("ROIPropertiesAlg"))
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireROI::~CalWireROI()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireROI::reconfigure(fhicl::ParameterSet const& p)
  {
    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    bool doInducedChargeDeconv = false;
    std::vector<std::vector<size_t> > respNums = sss->GetNResponses();
    for (size_t i = 0; i < respNums.at(1).size(); i++) {
      if (respNums.at(1).at(i) > 1) {
        doInducedChargeDeconv = true;
      }
    }

    // Throw exception if deconvolution should include dynamic induced charge effects (not yet implemented in CalROI) - M. Mooney
    if (doInducedChargeDeconv == true) {
      throw art::Exception(art::errors::Configuration)
        << "CalWireROI can not yet handle deconvolution with dynamic induced charge effects turned on.  Please use CalWireMicroBooNE instead.";
    }

    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;

    fDigitModuleLabel     = p.get< std::string >                   ("DigitModuleLabel", "daq");
    fNoiseSource          = p.get< unsigned short >                ("NoiseSource", 3);
    fNumBinsHalf          = p.get< unsigned short >                ("NumBinsHalf", 3);
    fThreshold            = p.get< std::vector<float > >           ("Threshold");
    uin                   = p.get< std::vector<unsigned short> >   ("uPlaneROIPad");
    vin                   = p.get< std::vector<unsigned short> >   ("vPlaneROIPad");
    zin                   = p.get< std::vector<unsigned short> >   ("zPlaneROIPad");
    fDoBaselineSub        = p.get< bool >                          ("DoBaselineSub");
    fFFTSize              = p.get< int  >                          ("FFTSize");
    fSaveWireWF           = p.get< int >                           ("SaveWireWF");
    fMaxAllowedChanStatus = p.get< int >                           ("MaxAllowedChannelStatus");

    fDoBaselineSub_WaveformPropertiesAlg = p.get< bool >("DoBaselineSub_WaveformPropertiesAlg", "");
        
    if(uin.size() != 2 || vin.size() != 2 || zin.size() != 2) {
      throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }
    
    if(fNoiseSource > 2) {
      throw art::Exception(art::errors::Configuration)
      << "Invalid NoiseSource "<<fNoiseSource<<" Should be 0, 1 or 2";
    }

    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = uin[0];
    fPostROIPad[0] = uin[1];
    fPreROIPad[1]  = vin[0];
    fPostROIPad[1] = vin[1];
    fPreROIPad[2]  = zin[0];
    fPostROIPad[2] = zin[1];
    
    for(unsigned short ipl = 0; ipl < 3; ++ipl) {
      if(fPreROIPad[ipl] < 20 || fPostROIPad[ipl] < 20) throw art::Exception(art::errors::Configuration)
        <<"ROIPads must be at least 20";
    }
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  } // reconfigure

  //////////////////////////////////////////////////////
  void CalWireROI::reconfFFT(int temp_fftsize){
    // re-initialize the FFT service for the requested size
    art::ServiceHandle<util::LArFFT> fFFT;

    if(fFFT->FFTSize() >= temp_fftsize) return;

    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT(temp_fftsize, options, fitbins);
  } // reconfFFT

  //////////////////////////////////////////////////////
  void CalWireROI::beginJob()
  {
    fEventCount = 0;
  } // beginJob

  //////////////////////////////////////////////////////
  void CalWireROI::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireROI::produce(art::Event& evt)
  {      
    
    //get pedestal conditions
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();
  
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    reconfFFT(fFFTSize);
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    unsigned int bin(0);     // time bin loop variable
    
    std::unique_ptr<filter::ChannelFilter> chanFilt(new filter::ChannelFilter());

    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    double DeconNorm = sss->GetDeconNorm();

    // number of bins in the running average
    float nBins = 2 * fNumBinsHalf + 1;
    float sqrt_nBins = sqrt(nBins);

    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      
      // vector that will be moved into the Wire object
      recob::Wire::RegionsOfInterest_t ROIVec;
      
      // the starting position and length of each ROI in the packed holder vector
      std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
      // vector of ROI begin and end bins
      std::vector<std::pair<unsigned int, unsigned int>> rois;
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      
      // Testing an idea about rejecting channels
      if (digitVec->GetPedestal() < 0.) continue;

      unsigned int dataSize = digitVec->Samples();
      // vector holding uncompressed adc values
      std::vector<short> rawadc(dataSize);
      std::vector<float> adcPedSub(dataSize);
      
      std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
      unsigned int thePlane = wids[0].Plane;
      //unsigned int theWire = wids[0].Wire;
      
      // skip bad channels identified by ChannelFilter or use already filtered channel?
      if(fMaxAllowedChanStatus >= 0) {
        // The following test is meant to be temporary until the "correct" solution is implemented
        if(chanFilt->GetChannelStatus(channel) > fMaxAllowedChanStatus) continue;
      }

      // uncompress the data
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      adcPedSub.resize(rawadc.size());

      float pdstl = pedestalRetrievalAlg.PedMean(channel);
      
      float noiseRMS = 0;
      if(fNoiseSource == 1) {
        noiseRMS = digitVec->GetSigma();
     } else if(fNoiseSource == 2) {
        noiseRMS = sss->GetRawNoise(channel);
      } else if(fNoiseSource == 3) {
        noiseRMS = std::max((float)digitVec->GetSigma(), (float)sss->GetRawNoise(channel));
      }
      
      // loop over all adc values and subtract the pedestal
      for(bin = 0; bin < adcPedSub.size(); ++bin) adcPedSub[bin] = (float)rawadc[bin] - pdstl;

      // calculate average RMS using all bins in the sum
      noiseRMS /= sqrt_nBins;
      // Convert running average into a sum to minimize calculation
      // Set the noise floor
      float threshold = nBins * (3 * noiseRMS + fThreshold[thePlane]);
      
      double fourOverSqrt10 = 4 / sqrt(10.);
      
      // On the first entry, when bin = fNumBinsHalf, sum the ADC's between bin 0 and 2 * fNumBinsHalf.
      // After bin is incremented by 1, subtract the ADC value in bin 0:
      unsigned short subBin = 0;
      // and add the ADC value in bin 2 * fNumBinsHalf + 1:
      unsigned short addBin = 2 * fNumBinsHalf + 1;
      // These variables will be incremented on each iteration rather than recalculating
      
      // sum of pedestal subtracted ADC values in the range (bin-fNumBinsHalf) to (bin+fNumBinsHalf)
      float adcSum = 0;
     
      // The starting bin of a ROI. It is set to 0 after leaving a ROI
      unsigned int roiStart = 0;
      
      // Find the starting ADC sum
      for(bin = 0; bin < addBin; ++bin) adcSum += std::abs(adcPedSub[bin]);
     
      for(bin = fNumBinsHalf + 1; bin < dataSize - fNumBinsHalf - 1; ++ bin) {
        // subtract the ADC value that is now outside the range
        ++subBin;
        adcSum -= std::abs(adcPedSub[subBin]);
        // add the ADC value that is now in the range
        ++addBin;
        adcSum += std::abs(adcPedSub[addBin]);
        if(roiStart > 0) {
          // In a ROI and leaving it?
          if(adcSum < threshold) {
            rois.push_back(std::make_pair(roiStart, bin));
            roiStart = 0;
          } // adcSum < threshold
        } else {
          // Not in a ROI but entering?
          if(adcSum > threshold) roiStart = bin;
        }
      } // bin
      // add the last ROI if necessary
      if (roiStart!=0) rois.push_back(std::make_pair(roiStart, dataSize-1));

        // skip deconvolution if there are no ROIs
      if(rois.size() == 0) continue;
      
      holderInfo.clear();
      
      // pad the ROIs
      for(unsigned int ii = 0; ii < rois.size(); ++ii) {
        // low ROI end
        int low = rois[ii].first - fPreROIPad[thePlane];
        if(low < 0) low = 0;
        rois[ii].first = low;
        // high ROI end
        unsigned int high = rois[ii].second + fPostROIPad[thePlane];
        if(high >= dataSize) high = dataSize-1;
        rois[ii].second = high;
      }
      
      // merge the ROIs?
      if(rois.size() > 1) {
        // temporary vector for merged ROIs
        std::vector<std::pair<unsigned int, unsigned int>> trois;
        for (unsigned int ii = 0; ii<rois.size();ii++){
          unsigned int roiStart = rois[ii].first;
          unsigned int roiEnd = rois[ii].second;
          int flag1 = 1;
          unsigned int jj=ii+1;
          while(flag1){
            if (jj<rois.size()) {
              if(rois[jj].first <= roiEnd ) {
                roiEnd = rois[jj].second;
                ii = jj;
                jj = ii+1;
              } else {
                flag1 = 0;
              }
            } else {
              flag1 = 0;
            }
          } // flag1
          trois.push_back(std::make_pair(roiStart,roiEnd));
        } // ii
        rois = trois;
      } // rois.size() > 1
      
      // pack the ROI's and deconvolve
      for (unsigned int ir = 0; ir < rois.size(); ++ir) {
        unsigned int roiLen = rois[ir].second - rois[ir].first + 1;
        unsigned int roiStart = rois[ir].first;
        int flag = 1;
        float tempPre=0,tempPost=0;
        std::vector<float> holder;
        while(flag){
          unsigned int transformSize = fFFTSize; //current transformsize
          //if ROI length is longer, take ROI length
          if (roiLen > transformSize) transformSize = roiLen;
          // Get signal shaping service.
          sss->SetDecon(transformSize);
          transformSize = fFFT->FFTSize();
          // temporary vector of signals
          holder.resize(transformSize,0);
          unsigned int hBin = 0;
          for(unsigned int bin = roiStart; bin < roiStart + holder.size(); ++bin) {
            if (bin < dataSize){
              holder[hBin] = adcPedSub[bin];
            } else {
              holder[hBin] = rawadc[bin-dataSize]-pdstl;
            }
            if (bin>=dataSize-1) flag = 0;
            ++hBin;
          } // bin
          sss->Deconvolute(channel,holder);
          for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;
          // 1. Check Baseline match
          // If not, include next ROI(if none, go to the end of signal)
          // If yes, proceed
          tempPre=0,tempPost=0;
          for(unsigned int bin = 0; bin < 20; ++bin) {
            tempPre  += holder[bin];
            tempPost += holder[roiLen -1 - bin];
          } // bin
          tempPre = tempPre/20.;
          tempPost = tempPost/20.;
          
          // Xin: Average over 10 bins and make a 4 sigma cut
          double deconNoise = sss->GetDeconNoise(channel) * fourOverSqrt10;
          
          if (fabs(tempPost-tempPre) < deconNoise){
            flag = 0;
          } else {
            if (tempPre > tempPost && roiStart <= 2){
              //if (tempPre > tempPost){
              flag = 0;
            } else {
              //		ir++;
              if (ir+1<rois.size()) {
                roiLen += 100;
                if (roiLen >= rois[ir+1].first - roiStart + 1)
                  roiLen = rois[++ir].second - roiStart + 1;
              }else {
                roiLen += 100;
                if (roiLen>dataSize-roiStart)
                  roiLen = dataSize - roiStart;
              }
            } // !tempPre > tempPost && roiStart <= 2
          } // !fabs(tempPost-tempPre) < deconNoise
        } // while(flag)
        // transfer the ROI and start bins into the vector that will be
        // put into the event
        std::vector<float> sigTemp;
        unsigned int bBegin = 0;
        //unsigned int theROI =ir;
        unsigned int bEnd = bBegin + roiLen;
        float basePre = 0., basePost = 0.;
        
        float base=0;
        if(fDoBaselineSub && fPreROIPad[thePlane] > 0 ) {
          basePre =tempPre;
          basePost=tempPost;
          base = SubtractBaseline(holder, basePre,basePost,roiStart,roiLen,dataSize);
        } // fDoBaselineSub ...
        else if(fDoBaselineSub_WaveformPropertiesAlg)
        {
          holder.resize(roiLen);
          base = fROIPropertiesAlg.GetWaveformPedestal(holder);
        }
        for(unsigned int jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]-base);
        // add the range into ROIVec
        ROIVec.add_range(roiStart, std::move(sigTemp));
      } // ir
      
      // create the new wire directly in wirecol
      wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
      // add an association between the last object in wirecol
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
        throw art::Exception(art::errors::InsertFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
      } // if failed to add association
    //  DumpWire(wirecol->back()); // for debugging
    }

//    if(wirecol->size() == 0) mf::LogWarning("CalWireROI") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
      art::ServiceHandle<art::TFileService> tfs;
      for (size_t wireN = 0; wireN < wirecol->size(); ++wireN){
        std::vector<float> sigTMP = wirecol->at(wireN).Signal();
        TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);",
                                      sigTMP.size(),-0.5,sigTMP.size()-0.5);
        for (size_t tick = 0; tick < sigTMP.size(); ++tick){
          fWire->SetBinContent(tick+1, sigTMP.at(tick) );
        }
      }
    } // fSaveWireWF
    

    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
  } // produce


  float CalWireROI::SubtractBaseline(std::vector<float>& holder, float basePre,
				     float basePost,unsigned int roiStart,
				     unsigned int roiLen,unsigned int dataSize)
  {
    float base=0;

    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20) {
      base = basePost;
      // can not trust the later part
    } else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20) {
      base = basePre;
      // can trust both
    } else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20) {
      if (fabs(basePre - basePost)<3){
        base = (basePre + basePost)/2.;
      } else {
        if (basePre < basePost) {
          base = basePre;
        } else {
          base = basePost;
        }
      }
      // can not use both
    } else {
      // Find the most probable value
      float min = 0,max = 0;
      for (unsigned int bin = 0; bin < roiLen; bin++) {
        if (holder[bin] > max) max = holder[bin];
        if (holder[bin] < min) min = holder[bin];
      }
      int nbin = max - min;
      if (nbin!=0){
        TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
        for (unsigned int bin = 0; bin < roiLen; bin++) {
          h1->Fill(holder[bin]);
        }
        float ped = h1->GetMaximum();
        float ave=0,ncount = 0;
        for (unsigned short bin = 0; bin < roiLen; bin++) {
          if (fabs(holder[bin]-ped)<2){
            ave +=holder[bin];
            ncount ++;
          }
        }
        if (ncount==0) ncount=1;
        ave = ave/ncount;
        h1->Delete();
        base = ave;
      }
    }
   
    return base;
  }


  void CalWireROI::doDecon(std::vector<float>& holder, 
    raw::ChannelID_t channel, unsigned int thePlane,
    std::vector<std::pair<unsigned int, unsigned int> > rois,
    std::vector<std::pair<unsigned int, unsigned int> > holderInfo,
    recob::Wire::RegionsOfInterest_t& ROIVec,
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss)
  {
    
    sss->Deconvolute(channel,holder);
   
    // transfer the ROIs and start bins into the vector that will be
    // put into the event
    for(unsigned int jr = 0; jr < holderInfo.size(); ++jr) {
      std::vector<float> sigTemp;
      unsigned int bBegin = holderInfo[jr].first;
      unsigned int theROI = holderInfo[jr].second;
      unsigned int roiLen = rois[theROI].second - rois[theROI].first;
      unsigned int bEnd = bBegin + roiLen;
      float basePre = 0., basePost = 0.;
      // Baseline subtraction if requested and the ROIs are padded.
      // Can't baseline subtract signals when the ROI start is close to 0 either
      if(fDoBaselineSub && fPreROIPad[thePlane] > 0 && 
          rois[theROI].first > fPreROIPad[thePlane]) {
        // find the baseline from the first few bins in the leading Pad region
        unsigned short bbins = fPreROIPad[thePlane];
        unsigned int bin;
        if(bbins > 5) bbins = 5;
        for(bin = 0; bin < bbins; ++bin) {
          basePre  += holder[bBegin + bin];
          basePost += holder[bEnd - bin];
        }
        basePre /= (float)bbins;
        basePost /= (float)bbins;
        float slp = (basePost - basePre) / (float)(roiLen - bbins);
        float base;
        for(unsigned int jj = bBegin; jj < bEnd; ++jj) {
          base = basePre + slp * (jj - bBegin);
          sigTemp.push_back(holder[jj] - base);
        } // jj
      } // fDoBaselineSub ...
      else {
        for(unsigned int jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]);
      } // !fDoBaselineSub ...
      
      // add the range into ROIVec 
      ROIVec.add_range(rois[theROI].first, std::move(sigTemp));
    } // jr
  } // doDecon


} // end namespace caldata
