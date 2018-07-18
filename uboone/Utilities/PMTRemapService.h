/**
 * @file   PMTRemapService.h
 * @brief  MicroBooNE service that provides a map to correct PMT OpChannel numbers
 * @author Brandon Eberly (eberly@slac.stanford.edu)
 */

#ifndef PMTREMAPSERVICE_H
#define PMTREMAPSERVICE_H

// Framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "uboone/Utilities/PMTRemapProvider.h"

namespace util {

  /**
   \class PMTRemapService
   This service provides only a simple interface to a provider class
   */
  class PMTRemapService {
      
    public:
      using provider_type = PMTRemapProvider;
   
      /// Constructor
      PMTRemapService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
   
      /// Destructor
      ~PMTRemapService() = default;

      //retrieve provider
      PMTRemapProvider const& GetProvider() const
      { return fProvider; }
	
      PMTRemapProvider const* GetProviderPtr() const
      { return &fProvider; }
    
    private:

      PMTRemapProvider fProvider;
       
  }; // class PMTRemapService
} // namespace util


DECLARE_ART_SERVICE(util::PMTRemapService, LEGACY)

#endif
