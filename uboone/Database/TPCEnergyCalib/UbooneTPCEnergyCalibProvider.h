/**
 * \file UbooneTPCEnergyCalibProvider.h
 * 
 * \brief Class def header for a class UbooneTPCEnergyCalibProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONETPCENERGYCALIBPROVIDER_H
#define UBOONETPCENERGYCALIBPROVIDER_H

#include "larevt/CalibrationDBI/IOVData/IOVDataConstants.h"
#include "UboonedqdxCorrectionProvider.h"
#include "UboonePlaneCorrectionProvider.h"
#include "TPCEnergyCalibProvider.h"

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
  class UbooneTPCEnergyCalibProvider : public TPCEnergyCalibProvider {
  
    public:
    
      /// Constructors
      UbooneTPCEnergyCalibProvider(fhicl::ParameterSet const& p);
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Update each source provider
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve TPC dq/dx and dE/dx calibration information     
      float YZdqdxCorrection(int plane, float y, float z) const override;
      float YZdqdxCorrectionErr(int plane, float y, float z) const override;
      float DriftdqdxCorrection(int plane, float drift_coord) const override;
      float DriftdqdxCorrectionErr(int plane, float drift_coord) const override;
      
      float XdqdxCorrection(int plane, float x) const override;
      float XdqdxCorrectionErr(int plane, float x) const override;
      float YdqdxCorrection(int plane, float y) const override;
      float YdqdxCorrectionErr(int plane, float y) const override;
      float ZdqdxCorrection(int plane, float z) const override;
      float ZdqdxCorrectionErr(int plane, float z) const override;
      
      float XNormdqdxCorrection(int plane) const override;
      float XNormdqdxCorrectionErr(int plane) const override;
      float XShapedqdxCorrection(int plane, float x) const override;
      float XShapedqdxCorrectionErr(int plane, float x) const override;
      
      float ThetadqdxCorrection(int plane, float theta) const override;
      float ThetadqdxCorrectionErr(int plane, float theta) const override;
      float PhidqdxCorrection(int plane, float phi) const override;
      float PhidqdxCorrectionErr(int plane, float phi) const override;
      
      /// total dqdx correction
      float TotaldqdxCorrection(int plane, float x, float y, float z, float theta=0.0, float phi=0.0) const override;
      float TotaldqdxCorrectionErr(int plane, float x, float y, float z, float theta=0.0, float phi=0.0) const override;
      
      /// dEdx correction as a function of plane number
      float dEdxCorrection(int plane) const override;
      float dEdxCorrectionErr(int plane) const override;
      
      
    private:
      
      std::vector<std::shared_ptr<UboonedqdxCorrectionProvider> >     fXShapeProvider; 
      std::vector<std::shared_ptr<UboonedqdxCorrectionProvider> >     fYZProvider;
      
      UboonePlaneCorrectionProvider fXNormProvider; 
      
      UboonePlaneCorrectionProvider fdEdxProvider;
  };
}//end namespace lariov

#endif
