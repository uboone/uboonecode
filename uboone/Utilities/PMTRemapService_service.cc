/**
 * @file   PMTRemapService_service.cc
 * @brief  MicroBooNE service that provides a map to correct PMT OpChannel numbers
 * @author Brandon Eberly (eberly@slac.stanford.edu)
 */

#include "uboone/Utilities/PMTRemapService.h"

namespace util {

  PMTRemapService::PMTRemapService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg) :
    fProvider(pset.get<fhicl::ParameterSet>("PMTRemapProvider"))
  {
  
    //no callbacks or fhicl configuration required at this time
  }
}

DEFINE_ART_SERVICE(util::PMTRemapService)
  
