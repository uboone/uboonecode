/**
 * @file   SpaceChargeMicroBooNETestHelpers.h
 * @brief  Helper functions for support of SpaceChargeMicroBooNEService in
 *         LArSoft tests
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 24, 2017
 * 
 * This library is a pure header.
 * It requires linkage with:
 * 
 * * `uboone_SpaceCharge`
 * * `mf_MessageLogger`
 * * `mf_Utilities`
 * * `fhiclcpp`
 * 
 */

#ifndef UBOONE_SPACECHARGE_SPACECHARGEMICROBOONETESTHELPERS_H
#define UBOONE_SPACECHARGE_SPACECHARGEMICROBOONETESTHELPERS_H

// LArSoft libraries
#include "uboone/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "larcore/TestUtils/ProviderTestHelpers.h"

// framework and utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


namespace testing {
  
  /**
   * @brief Set up a `spacecharge::SpaceChargeMicroBooNE` from a parameter set.
   * @return pointer to a newly created and set up `SpaceChargeMicroBooNE`
   * 
   * This function specialization enables the support of `SetupProvider()`
   * methods of `testing::TesterEnvironment`.
   */
  template <>
  struct ProviderSetupClass<spacecharge::SpaceChargeMicroBooNE> {
    
    static std::unique_ptr<spacecharge::SpaceChargeMicroBooNE> setup
      (fhicl::ParameterSet const& pset)
      {
        // some feedback about whether we are using the right configuration;
        // remember that art calls the implementation of a service interface
        // a "service provider"; LArSoft's service provider is a different thing
        std::string ServiceImplPath;
        if (pset.get_if_present("service_provider", ServiceImplPath)) {
          std::string ServiceImplName = ServiceImplPath;
          size_t iSlash = ServiceImplPath.rfind('/');
          if (iSlash != std::string::npos)
            ServiceImplName.erase(0, iSlash + 1);
          
          if (ServiceImplName == "SpaceChargeServiceMicroBooNE") {
            LOG_TRACE("setupProvider")
              << "Verified service implementation for SpaceChargeService: '"
              << ServiceImplPath << "'";
          }
          else {
            mf::LogWarning("setupProvider")
              << "This set up is for a SpaceChargeMicroBooNE provider.\n"
              "Your configuration specifies a '" << ServiceImplPath
              << "' service implementation"
              " that is not known to use that provider.";
          }
        }
        
        //
        // create the new DetectorClocksStandard service provider
        //
        return std::make_unique<spacecharge::SpaceChargeMicroBooNE>(pset);
      } // setup()
    
  }; // ProviderSetupClass<SpaceChargeMicroBooNE>
  
  
  /**
   * @brief Environment setup helper for SpaceChargeMicroBooNE.
   * @tparam TestEnv type of environment to set up
   * @see simpleEnvironmentSetup
   * 
   * A service provider is set up in the environment, associated with the types
   * `spacecharge::SpaceChargeMicroBooNE` and `spacecharge::SpaceCharge`.
   * Its configuration is read from "services.SpaceChargeService".
   * 
   * The environment is expected to expose an interface equivalent to the one
   * of `testing::TesterEnvironment`.
   * 
   * This class specialisation enables the support of `SimpleProviderSetup()`
   * methods of `testing::TesterEnvironment`.
   * It should be possible to set up a testing environment by calling:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * env.SimpleProviderSetup<spacecharge::SpaceChargeMicroBooNE>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *     
   * The provider will be available from any of these two calls:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * env.Provider<spacecharge::SpaceChargeMicroBooNE>();
   * env.Provider<spacecharge::SpaceCharge>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  template <typename TestEnv>
  struct SimpleEnvironmentSetupClass<spacecharge::SpaceChargeMicroBooNE, TestEnv> {
    static spacecharge::SpaceChargeMicroBooNE* setup(TestEnv& env)
      {
        return SimpleEnvironmentStandardSetupByName
          <spacecharge::SpaceChargeMicroBooNE, spacecharge::SpaceCharge, TestEnv>
          (env, "SpaceChargeService");
      }
  };
  

} // namespace testing

  
#endif // UBOONE_SPACECHARGE_SPACECHARGEMICROBOONETESTHELPERS_H
