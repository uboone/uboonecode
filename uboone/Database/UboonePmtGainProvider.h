/**
 * \file UboonePmtGainProvider.h
 *
 * \ingroup WebDBI
 * 
 * \brief Class def header for a class UboonePmtGainProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONEPMTGAINPROVIDER_H
#define UBOONEPMTGAINPROVIDER_H

#include "larevt/CalibrationDBI/IOVData/PmtGain.h"
#include "larevt/CalibrationDBI/IOVData/Snapshot.h"
#include "larevt/CalibrationDBI/IOVData/IOVDataConstants.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "larevt/CalibrationDBI/Providers/DatabaseRetrievalAlg.h"

namespace lariov {

  /**
   * @brief Retrieves information: pmt gain
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
   */
  class UboonePmtGainProvider : public DatabaseRetrievalAlg, public PmtGainProvider {
  
    public:
    
      /// Constructors
      UboonePmtGainProvider(fhicl::ParameterSet const& p);
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Update Snapshot and inherited DBFolder if using database.  Return true if updated
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve gain information
      const PmtGain& PmtGainObject(DBChannelID_t ch) const;      
      float Gain(DBChannelID_t ch) const override;
      float GainErr(DBChannelID_t ch) const override;
      CalibrationExtraInfo const& ExtraInfo(DBChannelID_t ch) const override;
      
    private:
    
      DataSource::ds fDataSource;
          
      Snapshot<PmtGain> fData;
  };
}//end namespace lariov

#endif

