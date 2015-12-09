/**
 * \file UboonePmtCalibrationProvider.h
 *
 * \brief Class def header for a class UboonePmtCalibrationProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONEPMTCALIBRATIONPROVIDER_H
#define UBOONEPMTCALIBRATIONPROVIDER_H

#include "PmtCalibrationContainer.h"
#include "CalibrationDBI/IOVData/Snapshot.h"
#include "CalibrationDBI/IOVData/IOVDataConstants.h"
#include "CalibrationDBI/Providers/DatabaseRetrievalAlg.h"

namespace lariov {

  /**
   * @brief Retrieves information: pmt calibrations
   * 
   * Configuration parameters
   * =========================
   * 
   * - *DatabaseRetrievalAlg* (parameter set, mandatory): configuration for the
   *   database; see lariov::DatabaseRetrievalAlg
   * - *UseDB* (boolean, default: false): retrieve information from the database
   * - *UseFile* (boolean, default: false): retrieve information from a file;
   *   not implemented yet
   * - *DefaultAmplitude* (real, default: ): SPE amplitude returned 
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultAmplitudeErr* (real, default: ): SPE amplitude uncertainty returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultWidth* (real, default: ): SPE width returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultWidthErr* (real, default: 0.3): SPE width uncertainty returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultArea* (real, default: 0.0): SPE area returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultAreaErr* (real, default: 0.0): SPE area uncertainty returned
   *   when /UseDB/ and /UseFile/ parameters are false
   */
  class UboonePmtCalibrationProvider : public DatabaseRetrievalAlg {
  
    public:
    
      /// Constructors
      UboonePmtCalibrationProvider(fhicl::ParameterSet const& p);
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Default destructor
      ~UboonePmtCalibrationProvider() {}
      
      /// Update Snapshot and inherited DBFolder if using database.  Return true if updated
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve calibration information
      const PmtCalibrationContainer& ChannelInfo(DBChannelID_t ch) const;      
      float Amplitude(DBChannelID_t ch) const;
      float AmplitudeErr(DBChannelID_t ch) const;
      float Width(DBChannelID_t ch) const;
      float WidthErr(DBChannelID_t ch) const;
      float Area(DBChannelID_t ch) const;
      float AreaErr(DBChannelID_t ch) const;
      const std::vector<float>& AvWaveForm(DBChannelID_t ch) const;
      const std::vector<float>& AvWaveFormErr(DBChannelID_t ch) const;
           
      //hardcoded information about database folder - useful for debugging cross checks
      const unsigned int NCOLUMNS = 9;    
      const std::vector<std::string> FIELD_NAMES = {"channel", "amplitude", "amplitude_err", "area", "area_err" "width", "width_err", "avwaveform", "avwaveform_err"};
      const std::vector<std::string> FIELD_TYPES = {"unsigned int", "float", "float", "float", "float", "float", "float", "float[]", "float[]"};
      
    private:
    
      DataSource::ds fDataSource;
          
      Snapshot<PmtCalibrationContainer> fData;
  };
}//end namespace lariov

#endif

