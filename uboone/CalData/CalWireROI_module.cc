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
#include <fstream>
#include <random>

// ROOT libraries
#include "TH1D.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "art/Utilities/make_tool.h"
#include "uboone/CalData/DeconTools/IROIFinder.h"
#include "uboone/CalData/DeconTools/IDeconvolution.h"
#include "uboone/CalData/DeconTools/IBaseline.h"

///creation of calibrated signals on wires
namespace caldata {

class CalWireROI : public art::EDProducer
{
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
    
    std::string                                              fDigitModuleLabel;           ///< module that made digits
    std::string                                              fSpillName;                  ///< nominal spill is an empty string
                                                                                          ///< it is set by the DigitModuleLabel
                                                                                          ///< ex.:  "daq:preSpill" for prespill data
    unsigned short                                           fNoiseSource;                ///< Used to determine ROI threshold
    size_t                                                   fFFTSize;                    ///< FFT size for ROI deconvolution
    int                                                      fSaveWireWF;                 ///< Save recob::wire object waveforms
    size_t                                                   fEventCount;                 ///< count of event processed
    int                                                      fMinAllowedChanStatus;       ///< Don't consider channels with lower status
    
    std::unique_ptr<uboone_tool::IROIFinder>                 fROIFinder;
    std::unique_ptr<uboone_tool::IDeconvolution>             fDeconvolution;
    
    const geo::GeometryCore*                                 fGeometry = lar::providerFrom<geo::Geometry>();
    art::ServiceHandle<util::LArFFT>                         fFFT;
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> fSignalShaping;
    
  protected: 
    
}; // class CalWireROI

DEFINE_ART_MODULE(CalWireROI)
  
//-------------------------------------------------
CalWireROI::CalWireROI(fhicl::ParameterSet const& pset)
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
    // Get the handle for the FFT
    fFFT = art::ServiceHandle<util::LArFFT>();
    
    // Get signal shaping service.
    fSignalShaping = art::ServiceHandle<util::SignalShapingServiceMicroBooNE>();
    
    bool doInducedChargeDeconv = false;
    std::vector<std::vector<size_t> > respNums = fSignalShaping->GetNResponses();
    for (size_t i = 0; i < respNums.at(1).size(); i++)
    {
        if (respNums.at(1).at(i) > 1) {
            doInducedChargeDeconv = true;
        }
    }

    // Throw exception if deconvolution should include dynamic induced charge effects (not yet implemented in CalROI) - M. Mooney
    if (doInducedChargeDeconv == true)
    {
        throw art::Exception(art::errors::Configuration)
            << "CalWireROI can not yet handle deconvolution with dynamic induced charge effects turned on.  Please use CalWireMicroBooNE instead.";
    }
    
    fROIFinder = art::make_tool<uboone_tool::IROIFinder>        (p.get<fhicl::ParameterSet>("ROIFinder"));
    fDeconvolution = art::make_tool<uboone_tool::IDeconvolution>(p.get<fhicl::ParameterSet>("Deconvolution"));
    
    fDigitModuleLabel           = p.get< std::string >   ("DigitModuleLabel", "daq");
    fNoiseSource                = p.get< unsigned short >("NoiseSource",          3);
    fFFTSize                    = p.get< size_t >        ("FFTSize"                );
    fSaveWireWF                 = p.get< int >           ("SaveWireWF"             );
    fMinAllowedChanStatus       = p.get< int >           ("MinAllowedChannelStatus");
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos )
    {
        fSpillName = fDigitModuleLabel.substr( pos+1 );
        fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
	
    return;
}


void CalWireROI::reconfFFT(int temp_fftsize)
{
    if(fFFT->FFTSize() >= temp_fftsize) return;

    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT(temp_fftsize, options, fitbins);
}

//-------------------------------------------------
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
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();

    // get the FFT service to have access to the FFT size
    reconfFFT(fFFTSize);
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;
      
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();

        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!chanFilt.IsPresent(channel)) continue;

        // Testing an idea about rejecting channels
        if (digitVec->GetPedestal() < 0.) continue;
        
        // skip bad channels
        if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
            size_t dataSize = digitVec->Samples();
            
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            // uncompress the data
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            // loop over all adc values and subtract the pedestal
            // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
            float pedestal = pedestalRetrievalAlg.PedMean(channel);
            
            // Get the pedestal subtracted data, centered in the deconvolution vector
            std::vector<float> rawAdcLessPedVec(dataSize);
            
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),[pedestal](const short& adc){return std::round(float(adc) - pedestal);});
            
            // Recover a measure of the noise on the channel for use in the ROI finder
            float raw_noise = digitVec->GetSigma();
            
            // vector of candidate ROI begin and end bins
            uboone_tool::IROIFinder::CandidateROIVec candRoiVec;

            // Now find the candidate ROI's
            fROIFinder->FindROIs(rawAdcLessPedVec, channel, raw_noise, candRoiVec);
            
            fDeconvolution->Deconvolve(rawAdcLessPedVec, channel, candRoiVec, ROIVec);
       } // end if not a bad channel

        // create the new wire directly in wirecol
        wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
        
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName))
        {
            throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Can't associate wire #" << (wirecol->size() - 1)
                << " with raw digit #" << digitVec.key();
        } // if failed to add association
        //  DumpWire(wirecol->back()); // for debugging
    }

    if(wirecol->size() == 0)
      mf::LogWarning("CalWireROI") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
        art::ServiceHandle<art::TFileService> tfs;
        for (size_t wireN = 0; wireN < wirecol->size(); wireN++){
            std::vector<float> sigTMP = wirecol->at(wireN).Signal();
            TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);",
				      sigTMP.size(),-0.5,sigTMP.size()-0.5);
            for (size_t tick = 0; tick < sigTMP.size(); tick++){
                fWire->SetBinContent(tick+1, sigTMP.at(tick) );
            }
        }
    }
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce


} // end namespace caldata
