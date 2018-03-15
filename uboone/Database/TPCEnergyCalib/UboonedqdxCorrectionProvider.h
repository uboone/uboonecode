/**
 * \file UboonedqdxCorrectionProvider.h
 * 
 * \brief Class def header for a class UboonedqdxCorrectionProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONEDQDXCORRECTIONPROVIDER_H
#define UBOONEDQDXCORRECTIONPROVIDER_H

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
  class UboonedqdxCorrectionProvider : public DatabaseRetrievalAlg {
  
    public:
    
      /// Constructor
      UboonedqdxCorrectionProvider(fhicl::ParameterSet const& p);      
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Update Snapshot and inherited DBFolder if using database.  Return true if updated
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve calibration object
      const DqdxCorrection& DqdxCorrectionObject(const std::vector<float>&) const;      
      
      /// Retrieve calibration info
      float Correction(const std::vector<float>&) const;
      float CorrectionErr(const std::vector<float>&) const;
      
    private:
      
      /// Compares two std::vector<floats> and returns true if each element of v1 is 
      /// less than or equal to its corresponding element in v2
      static bool CompareFunction(const std::vector<float>& v1, const std::vector<float>& v2);
      
      /// Determine which bin number (element position in fCoordToBinMap) corresponds to the 
      /// input coord
      int GetBinNumber(const std::vector<float>& coord) const;
    
      DataSource::ds fDataSource;
          
      Snapshot<DqdxCorrection> fData;
      
      //convert from input coordinate to bin number used by Snapshot
      std::vector<float> fCoord_min;
      std::vector<float> fCoord_max;
      std::vector<unsigned int>   fNbins;
      std::vector<std::string> fBinEdgeNames;
      bool fIsFixedBinSize;
      std::vector<std::vector<float> > fCoordToBinMap;
  };
}//end namespace lariov

#endif
