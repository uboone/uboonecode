#ifndef UBOONEPMTCALIBRATIONPROVIDER_CXX
#define UBOONEPMTCALIBRATIONPROVIDER_CXX

#include "UboonePmtCalibrationProvider.h"
#include "CalibrationDBI/Providers/WebError.h"
#include "CalibrationDBI/IOVData/IOVDataConstants.h"

// art/LArSoft libraries
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "Geometry/Geometry.h"
#include "cetlib/exception.h"

namespace lariov {

  //constructor      
  UboonePmtCalibrationProvider::UboonePmtCalibrationProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
      
  void UboonePmtCalibrationProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    this->DatabaseRetrievalAlg::Reconfigure(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg"));
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());

    bool UseDB      = p.get<bool>("UseDB", false);
    bool UseFile    = p.get<bool>("UseFile", false);

    //priority:  (1) use db, (2) use table, (3) use defaults
    //If none are specified, use defaults
    if ( UseDB )      fDataSource = DataSource::Database;
    else if (UseFile) fDataSource = DataSource::File;
    else              fDataSource = DataSource::Default;

    if (fDataSource == DataSource::Default) {
      float default_amplitude     = p.get<float>("DefaultAmplitude");
      //float default_amplitude_err = p.get<float>("DefaultAmplitudeErr");
      float default_width         = p.get<float>("DefaultWidth");
      //float default_width_err     = p.get<float>("DefaultWidthErr");
      float default_area          = p.get<float>("DefaultArea");
      //float default_area_err      = p.get<float>("DefaultAreaErr");

      PmtCalibrationContainer DefaultPmt(0);

      DefaultPmt.SetAmplitude(default_amplitude);
      //DefaultPmt.SetAmplitudeErr(default_amplitude_err);
      DefaultPmt.SetWidth(default_width);
      //DefaultPmt.SetWidthErr(default_width_err);
      DefaultPmt.SetArea(default_area);
      //DefaultPmt.SetAreaErr(default_area_err); 
      
      art::ServiceHandle<geo::Geometry> geo;
      std::cout<<"Number of OpDets: "<<geo->NOpDets()<<std::endl;
      std::cout<<"Number of OpChannels: "<<geo->NOpChannels()<<std::endl;
      std::cout<<"Max OpChannel: "<<geo->MaxOpChannel()<<std::endl;
      for (unsigned int i=0; i<= geo->MaxOpChannel(); ++i) {
        if (geo->IsValidOpChannel(i)) {
	  std::cout<<"  OpChannel "<<i<<" is valid;"<<std::endl;
	}
      }
      
      for (unsigned int i=0; i<36; ++i) {
        std::cout<<"    OpDet"<<i<<":";
	for (unsigned int j=0; j<4; ++j) {
	  std::cout<<" "<<geo->OpChannel(i,j);
	}
	std::cout<<std::endl;
      }
      
      /*need loop over pmts in geomtry to add row to fData*/
           
    }
    else if (fDataSource == DataSource::File) {
      throw cet::exception("UboonePmtCalibrationProvider")
        << "UboonePmtCalibrationProvider: input from file not implemented yet\n";
      //need to implement
    }
  }


  bool UboonePmtCalibrationProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    std::vector<DBChannelID_t> channels;
    fFolder->GetChannelList(channels);
    for (auto it = channels.begin(); it != channels.end(); ++it) {
      
      double amplitude, width, area;
      std::vector<double> av_waveform;
      
      //double amplitude, amplitude_err, width, width_err, area, area_err;
      //std::vector<double> av_waveform, av_waveform_err;
      
      fFolder->GetNamedChannelData(*it, "amplitude",     amplitude);
      //fFolder->GetNamedChannelData(*it, "amplitude_err", amplitude_err);
      fFolder->GetNamedChannelData(*it, "width",      width);
      //fFolder->GetNamedChannelData(*it, "width_err",  width_err);
      fFolder->GetNamedChannelData(*it, "area",       area);
      //fFolder->GetNamedChannelData(*it, "area_err",   area_err);   
      fFolder->GetNamedChannelData(*it, "avwaveform",       av_waveform);
      //fFolder->GetNamedChannelData(*it, "avwaveform_err",   av_waveform_err);    

      PmtCalibrationContainer pg(*it);
      pg.SetAmplitude( (float)amplitude );
      //pg.SetAmplitudeErr( (float)amplitude_err );
      pg.SetWidth( (float)width );
      //pg.SetWidthErr( (float)width_err );
      pg.SetArea( (float)area );
      //pg.SetAreaErr( (float)area_err );
      pg.SetAvWaveform( av_waveform );
      //pg.SetAvWaveformErr( av_waveform_err );

      fData.AddOrReplaceRow(pg);
    }

    return true;

  }
  
  const PmtCalibrationContainer& UboonePmtCalibrationProvider::ChannelInfo(DBChannelID_t ch) const { 
    return fData.GetRow(ch);
  }
      
  float UboonePmtCalibrationProvider::Amplitude(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).Amplitude();
  }
  
  /*float UboonePmtCalibrationProvider::AmplitudeErr(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).AmplitudeErr();
  }*/
  
  float UboonePmtCalibrationProvider::Width(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).Width();
  }
  
  /*float UboonePmtCalibrationProvider::WidthErr(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).WidthErr();
  }*/

  float UboonePmtCalibrationProvider::Area(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).Area();
  }
  
  /*float UboonePmtCalibrationProvider::AreaErr(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).AreaErr();
  }*/
  
  const std::vector<double>& UboonePmtCalibrationProvider::AvWaveForm(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).AvWaveForm();
  }
  
  /*const std::vector<double>& UboonePmtCalibrationProvider::AvWaveFormErr(DBChannelID_t ch) const {
    return this->ChannelInfo(ch).AvWaveFormErr();
  }*/

}//end namespace lariov
	
#endif
        
