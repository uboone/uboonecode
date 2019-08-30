////////////////////////////////////////////////////////////////////////
//
// CalWireZS class - modification of CalWireROI that deconvolves zero suppressed 
// waveforms such as those produced by the SN stream
//
// jcrespo@nevis.columbia.edu
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <memory> // std::unique_ptr<>

// ROOT libraries
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
//#include "TCanvas.h"

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
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "art/Utilities/make_tool.h"
#include "uboone/CalData/DeconTools/IROIFinder.h"
#include "uboone/CalData/DeconTools/IDeconvolution.h"

///creation of calibrated signals on wires from zero-suppressed data
namespace caldata {

class CalWireZS : public art::EDProducer
{
  public:
    // create calibrated signals on wires from zero-suppressed data. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireZS(fhicl::ParameterSet const& pset); 
    virtual ~CalWireZS();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
    void reconfFFT(int temp_fftsize);
    
  private:
    // Subtract baseline from zero-suppressed ROIs and copy them to a vector of floats spanning the full wire readout. 
    // Also return a vector of pairs with the beginning/end of each ROI needed by the Deconvolve function
    void extractZSROIs(const recob::Wire::RegionsOfInterest_t&, const raw::ChannelID_t, uboone_tool::IROIFinder::CandidateROIVec&, std::vector<float>&) const;

    // Get parameters for baseline modeled as linear
    void getBaselineParams(const recob::Wire::RegionsOfInterest_t::range_const_iterator, const raw::ChannelID_t, float&, float&) const;
    
    std::string                                              fDigitModuleLabel;           ///< module that made zero-suppressed digits (stored as recob::Wires)
    std::string                                              fSpillName;                  ///< nominal spill is an empty string
                                                                                          ///< it is set by the DigitModuleLabel
                                                                                          ///< ex.:  "daq:preSpill" for prespill data
    size_t                                                   fFFTSize;                    ///< FFT size for ROI deconvolution
    int                                                      fSaveWireWF;                 ///< Save recob::wire object waveforms
    size_t                                                   fEventCount;                 ///< count of event processed
    size_t                                                   fEventID;                    ///< event number
    int                                                      fMinAllowedChanStatus;       ///< Don't consider channels with lower status // jcrespo: defined for SN?
    
    bool                                                     fOutputHistograms;           ///< Output histograms?
    
    std::unique_ptr<uboone_tool::IDeconvolution>             fDeconvolution;
    
    const geo::GeometryCore*                                 fGeometry = lar::providerFrom<geo::Geometry>();
    art::ServiceHandle<util::LArFFT>                         fFFT;

    bool fSubtractBaseline; ///< Subtract baseline?
    size_t fZSPresamples; ///< Number of samples saved before crossing the zero-suppression threshold
    size_t fZSPostsamples; ///< Number of samples saved after crossing the zero-suppression threshold
    float fMaxADCMedianCut; ///< Maximum abs ADC difference between a presample/postsample and the median to be considered as baseline
    bool fFixFlippingBits; ///< Try to fix flipping bits?
    float fMaxADCInterpolDiff; ///< Maximum abs ADC difference to the interpolation using nearest neighbors to be considered non-flipped bits
    float fMaxADCNNDiff; ///< Maximum abs ADC difference between nearest neighbors to be considered for interpolation
    
    // Define here a temporary set of histograms...
    std::vector<TH1D*> fNumROIsHistVec; // Number of ROIs after zero-suppression
    std::vector<TH1D*> fROILenHistVec; // Lenght of ROI after zero-suppression
    std::vector<TH1D*> fPedPresampleHistVec; // Pedestal distributions per plane using first presample
    std::vector<TH1D*> fPedPostsampleHistVec; // Pedestal distributions per plane using last postsample
    TH2D* fPedPresampleHist; // Pedestal distribution as a function of channel using first presample
    TH2D* fPedPostsampleHist; // Pedestal distribution as a function of channel using last postsample
    TH2D* fGoodPresamplesHist; // Number of presamples passing quality cuts to compute median
    TH2D* fGoodPostsamplesHist; // Number of postsamples passing quality cuts to compute median
    TH2D* fDiffPresamplesMedianHist; // Difference between presamples passing quality cuts and median
    TH2D* fDiffPostsamplesMedianHist; // Difference between postsamples passing quality cuts and median
}; // class CalWireZS

DEFINE_ART_MODULE(CalWireZS)
  
//-------------------------------------------------
CalWireZS::CalWireZS(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);

  produces< std::vector<recob::Wire> >(fSpillName);
  produces< art::Assns<raw::RawDigit, recob::Wire> >(fSpillName); // dummy association needed by hit finder
}

//-------------------------------------------------
CalWireZS::~CalWireZS()
{
}

//////////////////////////////////////////////////////
void CalWireZS::reconfigure(fhicl::ParameterSet const& p)
{

    fDeconvolution = art::make_tool<uboone_tool::IDeconvolution>(p.get<fhicl::ParameterSet>("Deconvolution"));
    
    fDigitModuleLabel           = p.get< std::string >   ("DigitModuleLabel", "sndaq"); // Module that produced the zero-suppressed digits (stored as recob::Wires)
    fFFTSize                    = p.get< size_t >        ("FFTSize"                );
    fSaveWireWF                 = p.get< int >           ("SaveWireWF"             );
    fMinAllowedChanStatus       = p.get< int >           ("MinAllowedChannelStatus");
    
    fOutputHistograms           = p.get< bool  >         ("OutputHistograms",   true);
    
    fSpillName.clear();

    // Add to FHiCL file
    fSubtractBaseline = p.get< bool >("SubtractBaseline", true);
    fZSPresamples = p.get< size_t > ("ZSPresamples", 7);
    fZSPostsamples = p.get< size_t > ("ZSPostsamples", 8);
    fMaxADCMedianCut = p.get< float > ("MaxADCMedianCut", 15);
    fFixFlippingBits = p.get< bool >("FixFlippingBits", true);
    fMaxADCInterpolDiff = p.get< float > ("MaxADCInterpolDiff", 32); 
    fMaxADCNNDiff = p.get< float > ("MaxADCNNDiff", 64); 
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos )
    {
        fSpillName = fDigitModuleLabel.substr( pos+1 );
        fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
    
        fNumROIsHistVec.resize(3);
        fROILenHistVec.resize(3);
	fPedPresampleHistVec.resize(3);   
	fPedPostsampleHistVec.resize(3);   
    
        for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
        {
            fNumROIsHistVec[planeIdx] = tfs->make<TH1D>( Form("NROISplane_%02zu",planeIdx), ";# ROIs;", 100, 0, 100 );
            fROILenHistVec[planeIdx] = tfs->make<TH1D>( Form("ROISIZEplane_%02zu",planeIdx), ";ROI size;", 500, 0, 500 );
	    fPedPresampleHistVec[planeIdx] = tfs->make<TH1D>( Form("PedPresamplePlane_%02zu", planeIdx), 
							      Form("Presample pedestal plane %02zu; ADC; Entries", planeIdx), 4096, 0, 4096);   
	    fPedPostsampleHistVec[planeIdx] = tfs->make<TH1D>( Form("PedPostsamplePlane_%02zu", planeIdx), 
							       Form("Postsample pedestal plane %02zu; ADC; Entries", planeIdx), 4096, 0, 4096);
       }
       fPedPresampleHist = tfs->make<TH2D>( "PedPresampleChannel", "Presample pedestal; Channel; ADC", 8256, 0, 8256, 4096, 0, 4096 );
       fPedPostsampleHist = tfs->make<TH2D>( "PedPostsampleChannel", "Postsample pedestal; Channel; ADC", 8256, 0, 8256, 4096, 0, 4096 );
       fGoodPresamplesHist = tfs->make<TH2D>( "GoodPresamplesChannel", "Good presamples; Channel; ADC", 
					     8256, 0, 8256, fZSPresamples + 1, 0, fZSPresamples + 1 );
       fGoodPostsamplesHist = tfs->make<TH2D>( "GoodPostsamplesChannel", "Good postsamples; Channel; ADC", 
					       8256, 0, 8256, fZSPostsamples + 1, 0, fZSPostsamples + 1) ;
       fDiffPresamplesMedianHist = tfs->make<TH2D>("DiffPresamplesMedianChannel", "Good presample - median; Channel; ADC", 8256, 0, 8256, 8192, -4096, 4096);
       fDiffPostsamplesMedianHist = tfs->make<TH2D>("DiffPostsamplesMedianChannel", "Good postsample - median; Channel; ADC", 8256, 0, 8256, 8192, -4096, 4096);
     }
    
    return;
}


void CalWireZS::reconfFFT(int temp_fftsize)
{
    if(fFFT->FFTSize() >= temp_fftsize) return;

    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT(temp_fftsize, options, fitbins);
}

//-------------------------------------------------
void CalWireZS::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void CalWireZS::endJob()
{
}


//////////////////////////////////////////////////////
void CalWireZS::produce(art::Event& evt)
{
    // get event number
    fEventID = evt.id().event();
    mf::LogPrint("CalWireZS") << "Producing event " << fEventID;

    // get the FFT service to have access to the FFT size
    reconfFFT(fFFTSize);

    // make a collection of (deconvolved) wires
    std::unique_ptr< std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set with dummy raw digits (needed by hit finder)
    std::unique_ptr< art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the non-deconvolved wire List object(s). 
    art::Handle< std::vector<recob::Wire> > digitVecHandle;
    
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if(!digitVecHandle->size()) return;
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec; // 
      
        // get the reference to the current non-deconvolved recob::Wire
        art::Ptr<recob::Wire> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();

        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!chanFilt.IsPresent(channel)) continue;

	// skip channels without ticks
	size_t dataSize = digitVec->NSignal();
	if( !dataSize ) continue;

	// skip bad channels 
        if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
	  // Vector to hold the expanded recob::Wire after baseline subtraction
	  std::vector<float> rawAdcLessPedVec(dataSize, 0); 

	  // vector of zero-suppressed ROI begin and end bins
	  uboone_tool::IROIFinder::CandidateROIVec candRoiVec;
	  
	  //mf::LogPrint("CalWireZS") << "Extracting zero-suppressed ROIs from channel " << channel;

          // Subtract baseline from zero-suppressed ROIs, try to correct flippling bits, 
	  // and copy them to a vector of floats, also filling a CandidateRoiVec with the ROI limits
	  extractZSROIs(digitVec->SignalROI(), channel, candRoiVec, rawAdcLessPedVec);
	  
	  if( !candRoiVec.size() ) continue;

	  // Do the deconvolution
	  fDeconvolution->Deconvolve(rawAdcLessPedVec, channel, candRoiVec, ROIVec);

	  // Make some histograms?
	  if (fOutputHistograms)
            {
	      // First up, determine what kind of wire we have
	      std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
	      const geo::PlaneID&      planeID = wids[0].planeID();
              
	      fNumROIsHistVec.at(planeID.Plane)->Fill(candRoiVec.size(), 1.);
              
	      for(const auto& pair : candRoiVec)
		fROILenHistVec.at(planeID.Plane)->Fill(pair.second-pair.first, 1.);
            }

       } // end if not a bad channel

	// Ensure wires keep the correct number of nominal samples. This is required because we do not have a RawDigit object to pass to the recob::WireCreator
	ROIVec.resize(dataSize); 

        // create the new wire directly in wirecol
        wirecol->push_back(recob::WireCreator(std::move(ROIVec), channel, fGeometry->View(channel)).move());

	// Dummy raw digits
	// Needed by producers downstream that require raw::RawDigit - recob::Wire associations
	art::Ptr< raw::RawDigit > dummyRawDigit;

        // add an association between the last object in wirecol
        // (that we just inserted) and dummyRawDigit
        if (!util::CreateAssn(*this, evt, *wirecol, dummyRawDigit, *WireDigitAssn, fSpillName))
        {
            throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Can't associate wire #" << (wirecol->size() - 1)
                << " with dummy raw digit #" << dummyRawDigit.key();
        } // if failed to add association

    } // end of loop over all wires

    if(wirecol->size() == 0)
      mf::LogWarning("CalWireZS") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
      art::ServiceHandle<art::TFileService> tfs;
      // std::cout << "Saving waveforms for event " << fEventID << std::endl;
      for (size_t wireN = 0; wireN < wirecol->size(); wireN++){
	recob::Wire& iwire = (*wirecol)[wireN];
	//std::cout << "Channel " << iwire.Channel() << " has " << iwire.NSignal() << " samples and " 
	//	  << iwire.SignalROI().n_ranges() << " ROIs" << std::endl;   
	if(iwire.NSignal() != 0){
	  std::vector<float> sigTMP = iwire.Signal();
	  TH1D* hWire = tfs->make<TH1D>(Form("WF_Evt%06zu_Ch%04u_decon", fEventID, iwire.Channel()), ";Tick;ADC",
					sigTMP.size(), 0, sigTMP.size());
	  for (size_t tick = 0; tick < sigTMP.size(); tick++){
	    hWire->Fill(tick, sigTMP[tick]);
	  }

	  // Plot deconvolved ROIs
	  // Loop over deconvolved ROIs within a wire
	  // size_t ctrROI = 0;
	  // for (auto dROI = iwire.SignalROI().begin_range(); dROI != iwire.SignalROI().end_range(); ++dROI){
	  //   const size_t firstTick = (*dROI).begin_index(); const size_t pastEndTick = (*dROI).end_index();
	  //   TH1D* hROI = tfs->make<TH1D>(Form("ROI_Evt%06zu_Ch%04u_deconroi%04zu", fEventID, iwire.Channel(), ctrROI),
	  // 				 ";Tick;ADC", pastEndTick - firstTick, firstTick, pastEndTick);
	  //   for (size_t iTick = firstTick; iTick < pastEndTick; iTick++){
	  //     hROI->Fill((int)iTick, (*dROI)[iTick]);
	  //   } 
	  //   // Uncomment this block to print deconvolved ROIs to png files
	  //   // TCanvas canvas(Form("c_Evt%06zu_Ch%04u_deconroi%04zu", fEventID, iwire.Channel(), ctrROI),
	  //   // 	      Form("c_Evt%06zu_Ch%04u_deconroi%04zu", fEventID, iwire.Channel(), ctrROI));
	  //   // hROI->Draw("hist ][");
	  //   // canvas.Modified();
	  //   // canvas.Update();
	  //   // canvas.Print(".png");
	  //   // ctrROI++;
	  // } // End of plot deconvolved ROIs
	}
      }
    }
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    mf::LogPrint("CalWireZS") << "Produced event " << fEventID;

    return;
} // produce

// Subtract baseline from zero-suppressed ROIs, try to correct flippling bits, and copy them to a vector of floats, also filling a CandRoiVec with the ROI limits
void CalWireZS::extractZSROIs( const recob::Wire::RegionsOfInterest_t & zsROIs, const raw::ChannelID_t channel, 
			       uboone_tool::IROIFinder::CandidateROIVec& roiLimits, std::vector<float>& waveform) const
{
    if(zsROIs.n_ranges() == 0) return; // do not waste time with empty waveforms

    art::ServiceHandle<art::TFileService> tfs;

    roiLimits.reserve(zsROIs.n_ranges());
    size_t ctrROI = 0;    

    // Loop over zero-suppressed ROIs within a wire
    for(auto iROI = zsROIs.begin_range(); iROI != zsROIs.end_range(); ++iROI){
      auto ROI = *iROI;
      const size_t firstTick = ROI.begin_index();
      const size_t pastEndTick = ROI.end_index();

      // Check for minimum ROI length
      if( pastEndTick - firstTick < fZSPresamples + 1 + fZSPostsamples ) continue;

      //roiLimits.push_back(uboone_tool::IROIFinder::CandidateROI(firstTick, pastEndTick - 1));
      roiLimits.push_back(uboone_tool::IROIFinder::CandidateROI(firstTick, pastEndTick));

      // Linear interpolation for baseline
      float slope = 0.;
      float intercept = 0.;
      if( fSubtractBaseline ) getBaselineParams(iROI, channel, slope, intercept);
      
      // Loop over ADCs within a zero-suppressed ROI and subtract baseline
      // Also correct flipping bits
      for (size_t iTick = firstTick; iTick < pastEndTick; iTick++ ){
	// --- Add cuts on bad ADCs here ---
	// Get rid of weird ADC values (swizzler errors?)
	if(ROI[iTick] > 1 && ROI[iTick] < 4096) waveform[iTick] = ROI[iTick]- slope*iTick - intercept;

	if( fFixFlippingBits ){
	  // Check for flipping bits and correct if needed
	  // Get neighbor ADCs for comparison to interpolation
	  float preADC = iTick > 0? waveform[iTick - 1] : 0.; // Assume 0 ADC for the sample preceding the 0th sample
	  float postADC = 0.;
	  if( iTick < pastEndTick - 1 ){ // If not at the last sample
	    // Check the neighbor ADC has not flipped bits
	    if( fabs(ROI[iTick + 1] - slope*(iTick + 1) - intercept - preADC) < 2.*fMaxADCNNDiff ){
	      postADC = ROI[iTick + 1]- slope*(iTick + 1) - intercept;
	      // Clara's bug: first sample over threshold is copied as last postsample
	      if( iTick + 1 == (pastEndTick - 1) && ROI[pastEndTick - 1] == ROI[firstTick + fZSPresamples] ) postADC = preADC; 
	    } else postADC = preADC;
	  }
	  // Interpolate to correct for flipped bits
	  waveform[iTick] = fabs(waveform[iTick] - (preADC + postADC)/2.) < fMaxADCInterpolDiff? waveform[iTick] : (preADC + postADC)/2.;
	  // Clara's bug: first sample over threshold is copied as last postsample
	  if(iTick == (pastEndTick - 1) && ROI[pastEndTick - 1] == ROI[firstTick + fZSPresamples]) waveform[iTick] = waveform[iTick - 1]; 
	}
      }

      // Low-level plot of each ROI before and after background subtraction
      // if(fSaveWireWF){ 
      // 	TH1D* hOrigROI = tfs->make<TH1D>(Form("ROI_Evt%06zu_Ch%04u_origroi%04zu", fEventID, channel, ctrROI),
      // 				       ";Tick;ADC", pastEndTick - firstTick, firstTick, pastEndTick);
      // 	TH1D* hSubROI = tfs->make<TH1D>(Form("ROI_Evt%06zu_Ch%04u_bgsubroi%04zu", fEventID, channel, ctrROI),
      // 				       ";Tick;ADC", pastEndTick - firstTick, firstTick, pastEndTick);
      // 	hOrigROI->SetLineColor(kBlack);
      // 	hSubROI->SetLineColor(kRed);	
      // 	for (size_t iTick = firstTick; iTick < pastEndTick; iTick++ ){
      // 	  //hOrigROI->Fill((int)iTick, ROI[iTick] - slope*firstTick - intercept); // Shift original waveform to compare to baseline subtracted
      // 	  // Shift original waveform by the average baseline to compare to baseline subtracted
      // 	  hOrigROI->Fill((int)iTick, ROI[iTick] - slope*(firstTick + pastEndTick - 1)/2. - intercept); 
      // 	  hSubROI->Fill((int)iTick, waveform[iTick]);
      // 	}
      // 	// Uncomment this block to print ROIs to png files	
      // 	// TCanvas canvas(Form("c_Evt%06zu_Ch%04u_zsroi%04zu", fEventID, channel, ctrROI), Form("c_Evt%06zu_Ch%04u_zsroi%04zu", fEventID, channel, ctrROI));
      // 	// hOrigROI->Draw("hist ][");
      // 	// hSubROI->Draw("hist ][ same");
      // 	// double ymax = (fabs(hSubROI->GetMaximum()) > fabs(hSubROI->GetMinimum()))? 1.1*hSubROI->GetMaximum() : 1.1*fabs(hSubROI->GetMinimum());
      // 	// hOrigROI->SetMaximum( ymax );
      // 	// hOrigROI->SetMinimum( -ymax );
      // 	// canvas.Modified();
      // 	// canvas.Update();
      // 	// canvas.Print(".png");
      // 	// canvas.Print(".root");
      // }
      
      ctrROI++;
    } // End of loop over zero-suppressed ROIs within a wire

    if(fSaveWireWF){
      TH1D* hWF = tfs->make<TH1D>(Form("WF_Evt%06zu_Ch%04u_bgsub", fEventID, channel),
				  ";Tick;ADC", waveform.size(), 0, waveform.size());
      for(size_t i = 0; i < waveform.size(); i++){
	hWF->Fill((int)i, waveform[i]);
      }
    }
} // end of CalWireZS::extractZSROIs

// Get parameters for baseline modeled as linear
void CalWireZS::getBaselineParams( const recob::Wire::RegionsOfInterest_t::range_const_iterator itROI, const raw::ChannelID_t channel,
				   float& slope, float& intercept ) const
{
    auto ROI = *itROI;
    const size_t firstTick = ROI.begin_index();
    const size_t pastEndTick = ROI.end_index();

    // Find the median of presamples/postsamples
    std::vector<float> presamples; 
    std::vector<float> postsamples; 
    presamples.reserve(fZSPresamples);
    postsamples.reserve(fZSPostsamples);
    for(size_t isample = 0; isample < fZSPresamples; isample++){
      // Get rid of weird ADC values (swizzler errors?)
      if(ROI[firstTick + isample] > 1 && ROI[firstTick + isample] < 4096){ 
	presamples.push_back(ROI[firstTick + isample]); 
      }
    }
    for(size_t isample = 0; isample < fZSPostsamples; isample++){
      // Get rid of weird ADC values (swizzler errors?)
      if(ROI[pastEndTick - 1 - isample] > 1 && ROI[pastEndTick - 1 - isample] < 4096){
	postsamples.push_back(ROI[pastEndTick - 1 - isample]);
      }
    }
    std::sort(presamples.begin(), presamples.end()); 
    std::sort(postsamples.begin(), postsamples.end()); 
    const size_t goodPresamples = presamples.size();
    const size_t goodPostsamples = postsamples.size();
    const float medianPresample = ((goodPresamples % 2) == 0)?
      (presamples[goodPresamples/2 - 1] + presamples[goodPresamples/2])/2. : presamples[goodPresamples/2];
    const float medianPostsample = ((goodPostsamples % 2) == 0)?
      (postsamples[goodPostsamples/2 - 1] + postsamples[goodPostsamples/2])/2. : postsamples[goodPostsamples/2];
    
    // Get the earliest presample and latest postsample which are close to the medians
    float prebaseline = -4095; // baseline estimated from presamples
    float postbaseline = -4095; // baseline estimated from postsamples
    size_t pretick = -1; // tick with baseline estimated from presamples
    size_t postick = -1; // tick with baseline estimated from postsamples
    for(size_t isample = 0; isample < fZSPresamples; isample++){
      if(fOutputHistograms && ROI[firstTick + isample] > 1 && ROI[firstTick + isample] < 4096) 
	fDiffPresamplesMedianHist->Fill(channel, ROI[firstTick + isample]- medianPresample);
      
      if(fabs(ROI[firstTick + isample]- medianPresample) < fMaxADCMedianCut){
	pretick = firstTick + isample;
	prebaseline = ROI[pretick];
	break;
      }
    }
    for(size_t isample = 0; isample < fZSPostsamples; isample++){
      // Clara's bug: first sample over threshold is copied as last postsample
      if( isample == 0 && ROI[pastEndTick - 1] == ROI[firstTick + fZSPresamples] ) continue;
      
      if(fOutputHistograms && ROI[pastEndTick - 1 - isample] > 1 && ROI[pastEndTick - 1 - isample] < 4096)
	fDiffPostsamplesMedianHist->Fill(channel, ROI[pastEndTick - isample]- medianPostsample);
      
      if(fabs(ROI[pastEndTick - 1 -  isample] - medianPostsample) < fMaxADCMedianCut){
	postick = pastEndTick - 1 - isample;
	postbaseline = ROI[postick];
	break;
      }
    }
      
    // Linear interpolation for baseline
    slope = (postbaseline - prebaseline)/(postick - pretick);
    // intercept = (postick*prebaseline - postbaseline*pretick)/(postick - pretick);
    intercept = prebaseline - slope*pretick;;
    
    // Make some histograms?
    if (fOutputHistograms){
      // First up, determine what kind of wire we have
      std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
      const geo::PlaneID&      planeID = wids[0].planeID();
      fPedPresampleHistVec.at(planeID.Plane)->Fill(prebaseline);
      fPedPostsampleHistVec.at(planeID.Plane)->Fill(postbaseline);
      fPedPresampleHist->Fill(channel, prebaseline);
      fPedPostsampleHist->Fill(channel, postbaseline);
      fGoodPresamplesHist->Fill(channel, goodPresamples);
      fGoodPostsamplesHist->Fill(channel, goodPostsamples);
    }
} //end of CalWireZS::getBaselineParams

} // end namespace caldata
