#ifndef UBOONEPMTCALIBRATIONSERVICE_H
#define UBOONEPMTCALIBRATIONSERVICE_H

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "UboonePmtCalibrationProvider.h"

namespace lariov{

  /**
     \class UboonePmtCalibrationService
     art service implementation of IPmtGainService.  Implements 
     a pmt gain retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UboonePmtCalibrationService {
  
    public:
    
      UboonePmtCalibrationService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UboonePmtCalibrationService(){}
      
      void PreProcessEvent(const art::Event& evt) {
        fProvider.Update(evt.time().value());
      }
      
      UboonePmtCalibrationProvider const& GetProvider() const {
        return fProvider;
      }
     
    private:
    
      UboonePmtCalibrationProvider fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE(lariov::UboonePmtCalibrationService, LEGACY)

#endif
