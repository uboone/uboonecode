#ifndef PEDESTALCALCULATORMODULE_H
#define PEDESTALCALCULATORMODULE_H

#include <string>
#include <vector>
#include <map>
#include <utility>

#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "SimpleTypesAndConstants/RawTypes.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

/*!
 * Title:   PedestalCalculator class
 * Authors:  eberly@fnal.gov
 * Inputs:  raw::RawDigit
 * Outputs: Text file of pedestals, and optionally new RawDigits with pedestal fields filled
 *
 * Description:
 * This producer calculates electronics pedestals.  A single value is produced  for each channel by using all RawDigits in 
 * the subrun.  Optionally, a new set of RawDigits can be created that stores the pedestal, calculated event-by-event, for 
 * each channel.
 */

namespace calibration {

  class PedestalCalculator : public art::EDProducer {
  
    public:
    
      explicit PedestalCalculator(fhicl::ParameterSet const& pset);
      ~PedestalCalculator();
      
      void produce (art::Event& evt) override;
      
      void beginJob() override;
      //void endJob() override;
      
      void beginRun(art::Run& run) override;
      //void endRun(art::Run const& run) override;
      
      void beginSubRun(art::SubRun& subrun) override;
      void endSubRun(art::SubRun& subrun) override;
      
    private:
    
      void FindPedestal(std::map<short, int>& adc_counts, double& ped, double& rms) const;
    
      std::map< std::uint32_t, std::map<short, int> > fRawPedestalInfo;
      std::map< std::uint32_t, std::vector< std::pair<double,double> > > fEventPedestalInfo;
      
      bool fMakeRawDigits;
      std::string fRawDigitModuleLabel;
  };
  
  //-----------------------------------------------------
  PedestalCalculator::PedestalCalculator(fhicl::ParameterSet const& pset)
    : EDProducer() {
    
    fEventPedestalInfo.clear();
    fRawPedestalInfo.clear();
    
    fMakeRawDigits       = pset.get<bool>("MakeRawDigits");
    fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
    
    produces< std::vector<raw::RawDigit> >();
  }
  
  //-----------------------------------------------------
  PedestalCalculator::~PedestalCalculator(){}
   
  //-----------------------------------------------------
  void PedestalCalculator::beginJob() {
    fEventPedestalInfo.clear();
    fRawPedestalInfo.clear();
  }
  
  
  //-----------------------------------------------------
  void PedestalCalculator::beginRun(art::Run& run) {
    fEventPedestalInfo.clear();
    fRawPedestalInfo.clear();
  }
  
  
  //-----------------------------------------------------
  void PedestalCalculator::beginSubRun(art::SubRun& subrun) {
    fEventPedestalInfo.clear();
    fRawPedestalInfo.clear();
  }
  
  void PedestalCalculator::endSubRun(art::SubRun& subrun) {
  
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg
      = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();
      
    const lariov::IChannelStatusProvider& chanFilt 
      = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();

    
    art::ServiceHandle<art::TFileService> tfs;
    std::string tag = "r"+std::to_string(subrun.run())+"_s"+std::to_string(subrun.subRun());
    TH1D* means_histo = tfs->make<TH1D>( ("means_"+tag).c_str(),"Event-wise Pedestal Means for All Channels", 4000, 0.0, 4000.0);
    TH1D* rmss_histo = tfs->make<TH1D>( ("RMSs_"+tag).c_str(),"Event-wise Pedestal RMSs for All Channels", 500, 0.0, 10.0);
    
    TH1D* single_mean_histo = tfs->make<TH1D>( ("single_mean_"+tag).c_str(),"Subrun Pedestal Mean for All Channels", 4000, 0.0, 4000.0);
    TH1D* single_rms_histo = tfs->make<TH1D>( ("single_RMS_"+tag).c_str(),"Subrun Pedestal RMS for All Channels", 500, 0.0, 10.0);
    
    TH1D* mean_diffs_histo = tfs->make<TH1D>( ("mean_diffs_"+tag).c_str(), "Subrun and Event-wise Pedestal Mean Differences",500, -100.0, 100.0);
    TH1D* rms_diffs_histo = tfs->make<TH1D>( ("rms_diffs_"+tag).c_str(), "Subrun and Event-wise Pedestal RMS Differences",500, -10.0, 10.0);
    
    TH1D* mean_diffs_goodhisto = tfs->make<TH1D>( ("mean_diffs_good_"+tag).c_str(), "Good Channel Subrun and Event-wise Pedestal Mean Differences",500, -100.0, 100.0);
    TH1D* rms_diffs_goodhisto = tfs->make<TH1D>( ("rms_diffs_good_"+tag).c_str(), "Good Channel Subrun and Event-wise Pedestal RMS Differences",500, -10.0, 10.0);
    
    TH1D* mean_dbdiffs_histo = tfs->make<TH1D>( ("mean_dbdiffs_"+tag).c_str(), "Subrun and DB Pedestal Mean Differences",500, -100.0, 100.0);
    TH1D* rms_dbdiffs_histo = tfs->make<TH1D>( ("rms_dbdiffs_"+tag).c_str(), "Subrun and DB Pedestal RMS Differences",500, -10.0, 10.0);
    
    for (auto itC = fEventPedestalInfo.begin(); itC != fEventPedestalInfo.end(); ++itC) {
      
      //subrun pedestals
      raw::ChannelID_t channel = itC->first;
      std::map<short, int>& tmpMap = fRawPedestalInfo[channel];
      double pedestal, rms;
      FindPedestal(tmpMap, pedestal, rms);
      single_mean_histo->Fill(pedestal);
      single_rms_histo->Fill(rms);
      
      //compare subrun pedestals to current valid database pedestal
      double db_pedestal = pedestalRetrievalAlg.PedMean(channel);
      double db_rms = pedestalRetrievalAlg.PedRms(channel);
      mean_dbdiffs_histo->Fill(pedestal-db_pedestal);
	  
           
      //the rms of the mean of event-wise pedestals
      //this should be equivalent to the pedestal rms stored in the db (should be equiv. to mean error in future)
      /*double mean_rms = 0.0;
      double mean = 0.0;
      double sq_mean = 0.0;
      for (auto itP = itC->second.begin(); itP != itC->second.end(); ++itP) {
	mean += (long)(itP->first);
	sq_mean += (long)(itP->first)*(long)(itP->first);
      }
      mean_rms = sqrt(sq_mean/itC->second.size() - pow(mean/itC->second.size(),2.0));
      rms_dbdiffs_histo->Fill(mean_rms-db_rms);*/
      rms_dbdiffs_histo->Fill(rms-db_rms);
      
      //event-wise pedestals
      for (auto itP = itC->second.begin(); itP != itC->second.end(); ++itP) {
        means_histo->Fill(itP->first);
	rmss_histo->Fill(itP->second);
	
	mean_diffs_histo->Fill(pedestal - itP->first);
	//rms_diffs_histo->Fill(rms - sqrt(mean_rms*mean_rms + itP->second*itP->second));
	rms_diffs_histo->Fill(rms - itP->second);
	
	if (chanFilt.IsGood(channel)) {
	  mean_diffs_goodhisto->Fill(pedestal - itP->first);
	  //rms_diffs_goodhisto->Fill(rms - sqrt(mean_rms*mean_rms + itP->second*itP->second));
	  rms_diffs_goodhisto->Fill(rms - itP->second);
	  
	  if (fabs(pedestal - itP->first) > 3.5) {
	    std::cout<<"Channel "<<channel<<" pedestal diff is "<<pedestal-itP->first<<std::endl;
	  }
	}
      }
           
      //std::cout<<"Channel "<<channel<<"  Rms: "<<mean_rms<<" DB Rms: "<<db_rms<<"  Mean: "<<pedestal<<"  DB Ped: "<<db_pedestal<<std::endl;
      //std::cout<<"Channel "<<channel<<"  Rms: "<<rms<<" DB Rms: "<<db_rms<<"  Mean: "<<pedestal<<"  DB Ped: "<<db_pedestal<<std::endl;
    }
    
    //todo - make text file output

    
  }
    
  
  //-----------------------------------------------------
  void PedestalCalculator::produce(art::Event& evt) {
  
    //get raw digits
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel, rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);
    
    //new rawdigit collection
    std::unique_ptr< std::vector<raw::RawDigit> > new_rawdigits(new std::vector<raw::RawDigit>);
    new_rawdigits->reserve(rawDigitVector.size());
    
    //loop over channels (RawDigits)
    for (auto itR = rawDigitVector.begin(); itR != rawDigitVector.end(); ++itR) {
      raw::RawDigit rd = *itR;
      raw::ChannelID_t channel = rd.Channel();
      
      std::vector<short> raw_adc(rd.Samples());
      raw::Uncompress(rd.ADCs(), raw_adc, rd.Compression());
      
      std::map<short, int> num_counts;
      for (auto itADC = raw_adc.begin(); itADC != raw_adc.end(); ++itADC) {
	
	if ( !num_counts.emplace(*itADC, 1).second ) {
	  num_counts[*itADC]++;
	  fRawPedestalInfo[channel][*itADC]++;
	}
	else {
	  std::map<short, int> data;
	  data.emplace(*itADC, 1);
	  if (!fRawPedestalInfo.emplace(channel,data).second) {
	    fRawPedestalInfo[channel].emplace(*itADC,1);
	  }
	}
      }
      
      double pedestal, rms;
      FindPedestal(num_counts, pedestal, rms);
      fEventPedestalInfo[channel].push_back( std::pair<double,double>(pedestal,rms));
      
      if (fMakeRawDigits) {
        raw::RawDigit new_rd(channel, rd.Samples(), rd.ADCs(), rd.Compression());
	new_rd.SetPedestal(pedestal);
	new_rawdigits->push_back(std::move(new_rd));
      } 
          
    }//end loop over RawDigits
      
    //if (fMakeRawDigits) {
      evt.put(std::move(new_rawdigits));
    //}
    
    return;
  }
        

  //determine the pedestal to be the most frequent ADC value
  //If there is a tie, average if the difference is <= 2, 
  //choose the lowest if difference is > 2
  void PedestalCalculator::FindPedestal(std::map<short, int>& num_counts, double& pedestal, double& rms) const {

    pedestal = 0.0;
    rms = 0.0;
    
    //insert missing zeroes into numcounts.
    int old_adc = num_counts.begin()->first;
    for (int i=0; i!=3 && old_adc > 0; ++i) {
      num_counts[old_adc - 1] = 0;
      old_adc--;
    }
    for (auto itM = num_counts.begin(); itM != num_counts.end(); ++itM) {
      if (itM == num_counts.begin()) continue;

      int new_adc = itM->first;
      while (new_adc - old_adc > 1) {
	++old_adc;
	num_counts[old_adc] = 0;
      }
      old_adc = new_adc;
    } 
    for (int i=0; i!=4; ++i) {
      old_adc++;
      num_counts[old_adc] = 0;
    }

    //find candidate bin for pedestal mean
    int count = 0;
    short ped_candidate = -1;
    std::vector<short> pedestal_candidates;
    for (auto itNC = num_counts.begin(); itNC != num_counts.end(); ++itNC) {
      if (itNC->second > count) {
	pedestal_candidates.clear();
	count = itNC->second;
	pedestal_candidates.push_back(itNC->first);
      }
      else if (itNC->second == count) {
	pedestal_candidates.push_back(itNC->first);
      }
    }
    
    if (pedestal_candidates.size() == 1) {
      ped_candidate = pedestal_candidates.front();
    }
    else {
      std::sort(pedestal_candidates.begin(), pedestal_candidates.end()); 
      short diff = pedestal_candidates.back() - pedestal_candidates.front();
      if ( diff <= 2 ) {
	ped_candidate = pedestal_candidates.front() + diff/2;
      }
      else ped_candidate = pedestal_candidates.front();
    }

    //Find first bins below half max of peak on either side to 
    //define fit region.  Fit to determine pedestal and rms.
    std::map<short, int>::const_iterator itLow, itHigh;
    unsigned int low_counter = 0;
    for (auto itM = num_counts.find(ped_candidate); ;--itM) {
      itLow = itM;
      if ( itM->second < count/2 && low_counter >=3) {     
	break;
      }
      low_counter++;
    }
    unsigned int high_counter = 0;
    for (auto itM = num_counts.find(ped_candidate); itM != num_counts.end();++itM) {
      itHigh = itM;
      if (itM->second < count/2 && high_counter >=4) {      
	break;
      }
      high_counter++;
    }
    const unsigned int NBINS = itHigh->first - itLow->first;

    double xpts[NBINS];
    double ex[NBINS];
    double ypts[NBINS];
    double ey[NBINS];
    unsigned int itM_counter = 0;
    for (auto itM = itLow; itM != itHigh; ++itM) {
      xpts[itM_counter] = (double)(itM->first);
      ypts[itM_counter] = (double)(itM->second);
      ex[itM_counter] = 0.5;
      ey[itM_counter] = sqrt(ypts[itM_counter]);
      itM_counter++;
    }
    TGraphErrors g(NBINS, xpts, ypts, ex, ey);
    TFitResultPtr result = g.Fit("gaus","Q N EX0 S");
    pedestal = result->Value(1);
    rms = result->Value(2);
    
  }
 
  DEFINE_ART_MODULE(PedestalCalculator)
  
}//end namespace calibration



#endif
