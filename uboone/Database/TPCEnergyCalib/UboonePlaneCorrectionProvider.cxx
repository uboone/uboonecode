#ifndef UBOONEPLANECORRECTIONPROVIDER_CXX
#define UBOONEPLANECORRECTIONPROVIDER_CXX

#include "UboonePlaneCorrectionProvider.h"
#include "larevt/CalibrationDBI/Providers/WebError.h"

// art/LArSoft libraries
#include "cetlib/exception.h"
#include "larcore/Geometry/Geometry.h"


#include <fstream>

namespace lariov {

  //constructor      
  UboonePlaneCorrectionProvider::UboonePlaneCorrectionProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
      
  void UboonePlaneCorrectionProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
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
      float default_correction     = p.get<float>("DefaultCorrection");
      float default_correction_err = p.get<float>("DefaultCorrectionErr");

      DqdxCorrection defaultCalib(0);
      defaultCalib.SetCorrection(default_correction);
      defaultCalib.SetCorrectionErr(default_correction_err);
      
      art::ServiceHandle<geo::Geometry> geo;
      geo::plane_id_iterator itP = geo->begin_plane_id();
      for (; itP != geo->end_plane_id(); ++itP) {
	DBChannelID_t plane = itP.get()->ID().Plane;
	defaultCalib.SetChannel(plane);
	fData.AddOrReplaceRow(defaultCalib);
      }
      
    }
    else if (fDataSource == DataSource::File) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string abs_fp = sp.find_file(fileName);
      std::cout << "Using plane correction calibrations from local file: "<<abs_fp<<"\n";
      std::ifstream file(abs_fp);
      if (!file) {
        throw cet::exception("UboonePlaneCorrectionProvider")
          << "File "<<abs_fp<<" is not found.";
      }
      
      std::string line;
      DqdxCorrection dp(0);
      while (std::getline(file, line)) {
        size_t current_comma = line.find(',');
        DBChannelID_t plane = (DBChannelID_t)std::stoi(line.substr(0, current_comma));     
        float correction = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
        
        current_comma = line.find(',',current_comma+1);
        float correction_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	

        dp.SetChannel(plane);
        dp.SetCorrection(correction);
        dp.SetCorrectionErr(correction_err);
        
        fData.AddOrReplaceRow(dp);
      }
    }
    else {
      std::cout << "Using dqdx plane corrections from conditions database"<<std::endl;
    }
  }

  bool UboonePlaneCorrectionProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    std::vector<DBChannelID_t> planes;
    fFolder->GetChannelList(planes);
    for (auto it = planes.begin(); it != planes.end(); ++it) {

      double correction, correction_err;
      fFolder->GetNamedChannelData(*it, "correction",     correction);
      fFolder->GetNamedChannelData(*it, "correction_err", correction_err); 
      
      DqdxCorrection pg(*it);
      
      pg.SetCorrection( (float)correction );
      pg.SetCorrectionErr( (float)correction_err );

      fData.AddOrReplaceRow(pg);
    }

    return true;
  }
  
  const DqdxCorrection& UboonePlaneCorrectionProvider::DqdxCorrectionObject(int plane) const { 
    return fData.GetRow(plane);
  }
      
  float UboonePlaneCorrectionProvider::Correction(int plane) const {
    return this->DqdxCorrectionObject(plane).Correction();
  }
  
  float UboonePlaneCorrectionProvider::CorrectionErr(int plane) const {
    return this->DqdxCorrectionObject(plane).CorrectionErr();
  }


}//end namespace lariov
	
#endif
        
