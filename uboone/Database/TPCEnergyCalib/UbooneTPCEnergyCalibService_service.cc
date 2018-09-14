#ifndef UBOONETPCENERGYCALIBSERVICE_CC
#define UBOONETPCENERGYCALIBSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "TPCEnergyCalibService.h"
#include "UbooneTPCEnergyCalibProvider.h"
#include "uboone/Database/UbooneCalibrationServiceHelper.h"

namespace lariov{

  /**
     \class UbooneTPCEnergyCalibService
     art service implementation of TPCEnergyCalibService.  Implements 
     a TPC Energy (dq/dx and dE/dx) calibration retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UbooneTPCEnergyCalibService : public TPCEnergyCalibService {
  
    public:
    
      UbooneTPCEnergyCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UbooneTPCEnergyCalibService(){}
      
      void PreProcessEvent(const art::Event& evt);
     
    private:
    
      TPCEnergyCalibProvider const& DoGetProvider() const override {
        return fProvider;
      }   
      
      TPCEnergyCalibProvider const* DoGetProviderPtr() const override {
        return &fProvider; 
      }
    
      UbooneTPCEnergyCalibProvider fProvider;
      UbooneCalibrationServiceHelper fHelper;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneTPCEnergyCalibService, lariov::TPCEnergyCalibService, LEGACY)
      

namespace lariov{

  UbooneTPCEnergyCalibService::UbooneTPCEnergyCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("TPCEnergyCalibProvider")),
    fHelper(pset.get<fhicl::ParameterSet>("CalibrationHelper"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneTPCEnergyCalibService::PreProcessEvent);
  }
  
  void UbooneTPCEnergyCalibService::PreProcessEvent(const art::Event& evt) {
    
    fProvider.Update( fHelper.GetTimeStamp(evt, "TPC Energy Calibrations") );
  } 

}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneTPCEnergyCalibService, lariov::TPCEnergyCalibService)

#endif
