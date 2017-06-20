#ifndef UBOONEPMTGAINPROVIDER_CXX
#define UBOONEPMTGAINPROVIDER_CXX

#include "UboonePmtGainProvider.h"
#include "larevt/CalibrationDBI/Providers/WebError.h"

// art/LArSoft libraries
#include "cetlib/exception.h"
#include "larcore/Geometry/Geometry.h"


#include <fstream>

namespace lariov {

  //constructor      
  UboonePmtGainProvider::UboonePmtGainProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
      
  void UboonePmtGainProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    this->DatabaseRetrievalAlg::Reconfigure(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg"));
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());

    bool UseDB      = p.get<bool>("UseDB", false);
    bool UseFile    = p.get<bool>("UseFile", false);
    std::string fileName = p.get<std::string>("FileName", "");

    //priority:  (1) use db, (2) use table, (3) use defaults
    //If none are specified, use defaults
    if ( UseDB )      fDataSource = DataSource::Database;
    else if (UseFile) fDataSource = DataSource::File;
    else              fDataSource = DataSource::Default;

    if (fDataSource == DataSource::Default) {
      float default_gain         = p.get<float>("DefaultGain");
      float default_gain_err     = p.get<float>("DefaultGainErr");
      std::string default_gfit_status = "CONVERGED";
      
      float default_amp_gain     = p.get<float>("DefaultAmplitudeGain");
      float default_amp_gain_err = p.get<float>("DefaultAmplitudeGainErr");
      std::string default_afit_status = "CONVERGED";
      
      float default_ped_mean     = p.get<float>("DefaultBaselineMean");
      float default_ped_mean_err = 0.0;
      
      float default_ped_rms      = p.get<float>("DefaultBaselineRms");
      float default_ped_rms_err  = 0.0;

      PmtGain defaultGain(0);
      CalibrationExtraInfo extra_info("PmtGain");

      defaultGain.SetGain(default_gain);
      defaultGain.SetGainErr(default_gain_err);
      extra_info.AddOrReplaceStringData("gain_fit_status",default_gfit_status);
      extra_info.AddOrReplaceFloatData("amplitude_gain",default_amp_gain);
      extra_info.AddOrReplaceFloatData("amplitude_gain_err",default_amp_gain_err);
      extra_info.AddOrReplaceStringData("amplitude_gain_fit_status",default_afit_status);
      extra_info.AddOrReplaceFloatData("pedestal_mean",default_ped_mean);
      extra_info.AddOrReplaceFloatData("pedestal_mean_err",default_ped_mean_err);
      extra_info.AddOrReplaceFloatData("pedestal_rms",default_ped_rms);
      extra_info.AddOrReplaceFloatData("pedestal_rms_err",default_ped_rms_err);
      
      defaultGain.SetExtraInfo(extra_info);
      
      art::ServiceHandle<geo::Geometry> geo;
      for (unsigned int od=0; od!=geo->NOpDets(); ++od) {
        if (geo->IsValidOpChannel(od)) {
	  defaultGain.SetChannel(od);
	  fData.AddOrReplaceRow(defaultGain);
	}
      }
      
    }
    else if (fDataSource == DataSource::File) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string abs_fp = sp.find_file(fileName);
      std::cout << "Using pmt gains from local file: "<<abs_fp<<"\n";
      std::ifstream file(abs_fp);
      if (!file) {
        throw cet::exception("UboonePmtGainProvider")
          << "File "<<abs_fp<<" is not found.";
      }

      std::string line;
      PmtGain dp(0);
      while (std::getline(file, line)) {
        if (line[0] == '#') continue;
        size_t current_comma = line.find(',');
        DBChannelID_t ch = (DBChannelID_t)std::stoi(line.substr(0, current_comma));   
        float gain           = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
        
        current_comma = line.find(',',current_comma+1);
        float gain_err       = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
	
	current_comma = line.find(',',current_comma+1);
        std::string gain_fit_status     = line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1));
	
	current_comma = line.find(',',current_comma+1);
        float amp_gain       = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
        
        current_comma = line.find(',',current_comma+1);
        float amp_gain_err   = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
	
	current_comma = line.find(',',current_comma+1);
        std::string amp_gain_fit_status = line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1));
	
	current_comma = line.find(',',current_comma+1);
        float ped_mean       = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
        
        current_comma = line.find(',',current_comma+1);
        float ped_mean_err   = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
	
	current_comma = line.find(',',current_comma+1);
        float ped_rms        = std::stof( line.substr(current_comma+1, line.find(',',current_comma+1)-(current_comma+1)) );
        
        current_comma = line.find(',',current_comma+1);
        float ped_rms_err    = std::stof( line.substr(current_comma+1) );	

        CalibrationExtraInfo extra_info("PmtGain");
	extra_info.AddOrReplaceStringData("gain_fit_status",gain_fit_status);
	extra_info.AddOrReplaceFloatData("amplitude_gain",amp_gain);
	extra_info.AddOrReplaceFloatData("amplitude_gain_err",amp_gain_err);
	extra_info.AddOrReplaceStringData("amplitude_gain_fit_status",amp_gain_fit_status);
	extra_info.AddOrReplaceFloatData("pedestal_mean",ped_mean);
	extra_info.AddOrReplaceFloatData("pedestal_mean_err",ped_mean_err);
	extra_info.AddOrReplaceFloatData("pedestal_rms",ped_rms);
	extra_info.AddOrReplaceFloatData("pedestal_rms_err",ped_rms_err);

        dp.SetChannel(ch);
        dp.SetGain(gain);
        dp.SetGainErr(gain_err);
	dp.SetExtraInfo(extra_info);
	       
        fData.AddOrReplaceRow(dp);
      }
    }
    else {
      std::cout << "Using pmt gains from conditions database"<<std::endl;
    }
  }

  bool UboonePmtGainProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    std::vector<DBChannelID_t> channels;
    fFolder->GetChannelList(channels);
    for (auto it = channels.begin(); it != channels.end(); ++it) {

      double gain, gain_err, amp_gain, amp_gain_err, ped_mean, ped_mean_err, ped_rms, ped_rms_err;
      std::string gain_fit_status, amp_gain_fit_status;
      fFolder->GetNamedChannelData(*it, "gain",              gain);
      fFolder->GetNamedChannelData(*it, "gainerror",         gain_err); 
      fFolder->GetNamedChannelData(*it, "gainfitstatus",     gain_fit_status);
      fFolder->GetNamedChannelData(*it, "ampgain",           amp_gain);
      fFolder->GetNamedChannelData(*it, "ampgainerror",      amp_gain_err); 
      fFolder->GetNamedChannelData(*it, "ampgainfitstatus",  amp_gain_fit_status);
      fFolder->GetNamedChannelData(*it, "baselinemean",      ped_mean);
      fFolder->GetNamedChannelData(*it, "baselinemeanerror", ped_mean_err); 
      fFolder->GetNamedChannelData(*it, "baselinerms",       ped_rms);
      fFolder->GetNamedChannelData(*it, "baselinermserror",  ped_rms_err); 
      
      PmtGain pg(*it);
      CalibrationExtraInfo extra_info("PmtGain");
      extra_info.AddOrReplaceStringData("gain_fit_status",gain_fit_status);
      extra_info.AddOrReplaceFloatData("amplitude_gain",amp_gain);
      extra_info.AddOrReplaceFloatData("amplitude_gain_err",amp_gain_err);
      extra_info.AddOrReplaceStringData("amplitude_gain_fit_status",amp_gain_fit_status);
      extra_info.AddOrReplaceFloatData("pedestal_mean",ped_mean);
      extra_info.AddOrReplaceFloatData("pedestal_mean_err",ped_mean_err);
      extra_info.AddOrReplaceFloatData("pedestal_rms",ped_rms);
      extra_info.AddOrReplaceFloatData("pedestal_rms_err",ped_rms_err);
           
      pg.SetGain( (float)gain );
      pg.SetGainErr( (float)gain_err );
      pg.SetExtraInfo(extra_info);

      fData.AddOrReplaceRow(pg);
    }

    return true;
  }
  
  const PmtGain& UboonePmtGainProvider::PmtGainObject(DBChannelID_t ch) const { 
    return fData.GetRow(ch);
  }
      
  float UboonePmtGainProvider::Gain(DBChannelID_t ch) const {
    return this->PmtGainObject(ch).Gain();
  }
  
  float UboonePmtGainProvider::GainErr(DBChannelID_t ch) const {
    return this->PmtGainObject(ch).GainErr();
  }
  
  CalibrationExtraInfo const& UboonePmtGainProvider::ExtraInfo(DBChannelID_t ch) const {
    return this->PmtGainObject(ch).ExtraInfo();
  }


}//end namespace lariov
	
#endif
        
