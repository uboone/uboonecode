#ifndef CHANNELSTATUS_H
#define CHANNELSTATUS_H

#include <string>
#include <vector>
#include <iostream>
#include <TProfile2D.h>
#include <TTree.h>
#include <TMath.h>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "RawData/RawDigit.h"
#include "Utilities/DetectorProperties.h"

namespace calibration {

  class ChannelStatus : public art::EDAnalyzer {

  public:
    explicit ChannelStatus(fhicl::ParameterSet const& pset);
    virtual ~ChannelStatus();

    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void beginJob();
    void endJob();
    void GetChannelStatus(unsigned int chan, std::vector<short> const& adcVector, bool& isCut);

    //likely we will need begin/end run and subrun functions
    void beginRun(art::Run const& run);
    void endRun(art::Run const& run);
    void beginSubRun(art::SubRun const& subrun);
    void endSubRun(art::SubRun const& subrun);
    
  private:

    //******************************
    //Variables Taken from FHICL File
    std::string       fRawDigitModuleLabel;   //label for rawdigit module

    //Other variables
    //TTree* tStatus;
    //short fRun, fSubrun, fEvent, fChan, fStatus;

    TH1F *hData, *hFftData;
    TH2F *hRmsVsChan,*hFreqSumVsChan;
    TProfile2D *pFFTVsChan;
    TProfile *pMean, *pRms, *pLimit, *pFreqSum, *pStatus;
  }; //end class Noise


  //-------------------------------------------------------------------
  ChannelStatus::ChannelStatus(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
  }

  //-------------------------------------------------------------------
  ChannelStatus::~ChannelStatus(){}

  //-------------------------------------------------------------------
  void ChannelStatus::reconfigure(fhicl::ParameterSet const& pset){
    fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
  }

  //-------------------------------------------------------------------
  void ChannelStatus::beginJob(){
	art::ServiceHandle<art::TFileService> tfs;
	int fMaxTicks = 9594;
	hData = tfs->make<TH1F>("hData","",fMaxTicks,-0.5,fMaxTicks-0.5);
	hFftData = tfs->make<TH1F>("hFftData","",fMaxTicks,-0.5,fMaxTicks-0.5);
	pFFTVsChan = tfs->make<TProfile2D>("pFFTVsChan","",8256,-0.5,8256-0.5,100,0,0.5);

	hRmsVsChan = tfs->make<TH2F>("hRmsVsChan","",8256,-0.5,8256-0.5,200,0,20);
	hFreqSumVsChan = tfs->make<TH2F>("hFreqSumVsChan","",8256,-0.5,8256-0.5,200,0,20);

	pMean = tfs->make<TProfile>("pMean","",8256,-0.5,8256-0.5);
	pRms = tfs->make<TProfile>("pRms","",8256,-0.5,8256-0.5);
	pFreqSum = tfs->make<TProfile>("pFreqSum","",8256,-0.5,8256-0.5);
	pLimit = tfs->make<TProfile>("pLimit","",8256,-0.5,8256-0.5);
	pStatus = tfs->make<TProfile>("pStatus","",8256,-0.5,8256-0.5);
  }

  //-------------------------------------------------------------------
  void ChannelStatus::endJob(){
  }

  //-------------------------------------------------------------------
  void ChannelStatus::beginRun(art::Run const& run){

    art::ServiceHandle<art::TFileService> tfs;
    /*
    tStatus = tfs->make<TTree>("StatusTree","StatusTree");
    tStatus->Branch("run", &fRun, "run/s");
    tStatus->Branch("subrun", &fSubrun, "subrun/s");
    tStatus->Branch("event", &fEvent, "event/s");
    tStatus->Branch("chan", &fChan, "chan/s");
    tStatus->Branch("status", &fStatus, "status/s");
    */
  }

  //-------------------------------------------------------------------
  void ChannelStatus::endRun(art::Run const& run){
    //art::ServiceHandle<art::TFileService> tfs;
    //return;
  }


  //-------------------------------------------------------------------
  void ChannelStatus::beginSubRun(art::SubRun const& subrun){
  }

  //-------------------------------------------------------------------
  void ChannelStatus::endSubRun(art::SubRun const& subrun){
  }
  
  //-------------------------------------------------------------------
  void ChannelStatus::analyze(art::Event const& evt){

    art::ServiceHandle<art::TFileService> tfs;    
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    //loop over channels, get overall average waveform for event
    const unsigned int n_channels = rawDigitVector.size();
    if( n_channels == 0 )
	return;
    for(unsigned int ich=0; ich<n_channels; ich++){
	const size_t n_samp = rawDigitVector.at(ich).NADC();
	if( n_samp == 0 ) 
		continue;

	bool isCut = 0;
	GetChannelStatus( rawDigitVector.at(ich).Channel(), rawDigitVector.at(ich).ADCs(), isCut );
    }
  }
  
  void ChannelStatus::GetChannelStatus(unsigned int chan, std::vector<short> const& adcVector, bool& isCut){
	//calculate mean
    	double mean = 0;
    	int count = 0;
    	for(unsigned int it = 0 ; it < adcVector.size() ; it++ ){
		mean += adcVector.at(it);
		count++;
    	}
    	if( count > 0)
		mean = mean / (double) count;

	//calculate RMS
	count = 0;
    	double rms = 0;
    	for(unsigned int it = 0 ; it < adcVector.size() ; it++ ){
		rms += (adcVector.at(it) - mean)*(adcVector.at(it) - mean);
		count++;
    	}
    	if(count - 1 > 0)
		rms = TMath::Sqrt( rms /(double)(  count - 1  ) );
	
	//load vector into histogram for FFT
	hData->Reset();
	for(unsigned int it = 0 ; it < adcVector.size() ; it++ )
		hData->SetBinContent(it+1, adcVector.at(it) - mean );

	//perform FFT on histogram, get 0-100 kHz frequency sum
	hFftData->Reset();
	hData->FFT(hFftData,"MAG");
	double freqSum = 0;
    	for(int it = 1 ; it < hFftData->GetNbinsX() ; it++ ){
		double freq = 2.* it / (double) hFftData->GetNbinsX() ;
		pFFTVsChan->Fill( chan, freq, hFftData->GetBinContent(it+1) );
		if( freq < 0.1 )
			freqSum += hFftData->GetBinContent(it+1);
    	}

	//apply selection
	bool isCut1 = 0;
	bool isCut2 = 0;
	bool isCut3 = 0;

	//selection cut 1 : pedestal mean
	if( chan < 4800 && (mean < 2046 - 50 || mean > 2046 + 50) )
		isCut1 = 1;
	if( chan >= 4800 && (mean < 474 - 50 || mean > 474 + 50) )
		isCut1 = 1;

	//selection cut 2: noisy channel
	if( rms > 30 )
		isCut2 = 1;

	//selection cut 3: low RMS
	double limit = -1;
	if( chan < 2400 ){
		if( chan < 680 )
			limit = (chan - 0.)*(8.6 - 1.8 )/(680. - 0.) + 1.8;
		if( chan >= 680 && chan < 1728 )
			limit = 8.6;
		if( chan >= 1728 ) 
			limit = (chan - 1728)*( 1.7 - 8.6 )/(2400. - 1728.) + 8.6;
	}
	if( chan >= 2400 && chan < 4800){
		if( chan < 2976 )
			limit = (chan - 2400.)*(4.6 - 1.3 )/(2976.-2400.) + 1.3;
		if( chan >= 2976 && chan < 4136 )
			limit = 7.0;
		if( chan >= 4136 ) 
			limit = (chan - 4136)*( 1.5 - 6.8 )/(4800. - 4136.) + 6.8;
	}
	if( chan >= 4800  )
		limit = 4.6;
	limit = limit - 2;
	if( limit < 0 )	
		limit = 0;
	if( rms < limit )
		isCut3 = 1;

	hRmsVsChan->Fill(chan, rms);
	hFreqSumVsChan->Fill(chan, freqSum);

	pMean->Fill(chan, mean);
	pRms->Fill(chan, rms);
	pLimit->Fill(chan, limit);
	pFreqSum->Fill(chan, freqSum);
	pStatus->Fill(chan, isCut1+isCut2+isCut3);

	isCut = isCut1+isCut2+isCut3;
  }

  DEFINE_ART_MODULE(ChannelStatus)

} //end namespace ChannelStatus

#endif //CHANNELSTATUS_H
