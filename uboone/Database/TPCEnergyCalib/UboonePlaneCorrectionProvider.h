/**
 * \file UboonePlaneCorrectionProvider.h
 * 
 * \brief Class def header for a class UboonePlaneCorrectionProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONEPLANECORRECTIONPROVIDER_H
#define UBOONEPLANECORRECTIONPROVIDER_H

#include "DqdxCorrection.h"
#include "larevt/CalibrationDBI/IOVData/Snapshot.h"
#include "larevt/CalibrationDBI/IOVData/IOVDataConstants.h"
#include "larevt/CalibrationDBI/Providers/DatabaseRetrievalAlg.h"

namespace lariov {

  /**
   * @brief Retrieves information: electronics calibrations, specifically gain and shaping time
   * 
   * Configuration parameters
   * =========================
   * 
   * - *DatabaseRetrievalAlg* (parameter set, mandatory): configuration for the
   *   database; see lariov::DatabaseRetrievalAlg
   * - *UseDB* (boolean, default: false): retrieve information from the database
   * - *UseFile* (boolean, default: false): retrieve information from a file;
   *   not implemented yet
   * - *DefaultGain* (real, default: ): Gain returned 
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultGainErr* (real, default: ): Gain uncertainty returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultShapingTime* (real, default: ): Shaping Time returned 
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultShapingTimeErr* (real, default: ): Shaping Time uncertainty returned
   *   when /UseDB/ and /UseFile/ parameters are false
   */
  class UboonePlaneCorrectionProvider : public DatabaseRetrievalAlg {
  
    public:
    
      /// Constructor
      UboonePlaneCorrectionProvider(fhicl::ParameterSet const& p);      
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Update Snapshot and inherited DBFolder if using database.  Return true if updated
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve calibration object
      const DqdxCorrection& DqdxCorrectionObject(int plane) const;      
      
      /// Retrieve calibration info
      float Correction(int plane) const;
      float CorrectionErr(int plane) const;
      
    private:
    
      DataSource::ds fDataSource;
          
      Snapshot<DqdxCorrection> fData;
  };
}//end namespace lariov

#endif
