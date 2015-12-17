#ifndef UBOONEPMTCALIBRATIONSERVICE_CC
#define UBOONEPMTCALIBRATIONSERVICE_CC

#include "UboonePmtCalibrationService.h"

namespace lariov{

  UboonePmtCalibrationService::UboonePmtCalibrationService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("PmtCalibrationProvider"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UboonePmtCalibrationService::PreProcessEvent);
  }

}//end namespace lariov

DEFINE_ART_SERVICE(lariov::UboonePmtCalibrationService)

#endif
