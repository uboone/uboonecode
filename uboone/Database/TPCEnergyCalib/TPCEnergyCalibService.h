/**
 * @file   TPCEnergyCalibService.h
 * @brief  Interface for experiment-specific service for TPC dqdx and dEdx c alibrations
 * @author Brandon Eberly (eberly@slac.stanford.edu)
 */

#ifndef TPCENERGYCALIBSERVICE_H
#define TPCENERGYCALIBSERVICE_H

// Framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"

//forward declarations
namespace lariov {
  class TPCEnergyCalibProvider;
}

namespace lariov {

  /**
   \class TPCEnergyCalibService
   This service provides only a simple interface to a provider class
   */
  class TPCEnergyCalibService {
      
    public:
      using provider_type = TPCEnergyCalibProvider;
   
      /// Destructor
      virtual ~TPCEnergyCalibService() = default;

      //retrieve provider
      TPCEnergyCalibProvider const& GetProvider() const
      { return DoGetProvider(); }
	
      TPCEnergyCalibProvider const* GetProviderPtr() const
      { return DoGetProviderPtr(); }
    
    private:

      /// Returns a reference to the service provider
      virtual TPCEnergyCalibProvider const& DoGetProvider() const = 0;
      
      virtual TPCEnergyCalibProvider const* DoGetProviderPtr() const = 0;
    
    
    
  }; // class TPCEnergyCalibService
} // namespace lariov


DECLARE_ART_SERVICE_INTERFACE(lariov::TPCEnergyCalibService, LEGACY)

#endif
