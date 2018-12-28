/**
 * @file   SpaceChargeMicroBooNE_test.cc
 * @brief  "Unit test" for MicroBooNE implementation of space charge provider.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 24, 2017
 * 
 * This is not really a unit test, in that it does not test the results of the
 * service.
 * It still tests the correct loading of the test, and prints some information.
 * 
 * No configuration is passed from the configuration FHiCL file to the test
 * itself (the service still needs to be configured in there though).
 * If this test needs to be extended to get parameters from the FHiCL
 * configuration, please contact the author for directions.
 * 
 * @note This is *not* a good example on how to set up MicroBooNE geometry;
 *       I am getting away with using the standard mapping
 *       (`geo::ChannelMapStandardAlg`) instead of MicroBooNE's specific one
 *       and it happens to work _correctly_ as long as no specific feature of
 *       MicroBooNE's mapping is used. If instead those features are used, the
 *       results will be _silently_ wrong.
 * 
 */

// LArSoft libraries
#include "ubevt/SpaceCharge/SpaceChargeMicroBooNETestHelpers.h"
#include "ubevt/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandardTestHelpers.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcorealg/Geometry/geo_vectors_fhicl.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "test/Geometry/geometry_unit_test_base.h"

// framework libraries
#include "fhiclcpp/types/Table.h" 
#include "fhiclcpp/types/Atom.h" 
#include "fhiclcpp/types/Sequence.h" 
#include "fhiclcpp/types/Name.h" 
#include "fhiclcpp/types/Comment.h" 

// C/C++ standard libraries
#include <string>
#include <vector>


//------------------------------------------------------------------------------
//--- The test class
//---
class SCEtestAlgo {
  
  std::vector<geo::Point_t> fTestPoints;
  std::string fOutputCat;
  
  spacecharge::SpaceCharge const* pSC = nullptr;
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<lar::fhicl::geo::Point_t> testPoints{
      Name("testPoints"),
      Comment("points to be tested; each is a triplet [ x, y, z ] (cm)")
      };
    
    fhicl::Atom<std::string> outputCategory{
      Name("outputCategory"),
      Comment("name of the message facility stream for the test output"),
      "SCETest" // default
      };
    
  }; // Config
  
  using Parameters = fhicl::Table<Config>;
  
  /// Constructor from fhicl validation
  SCEtestAlgo(Parameters const& config)
    : fTestPoints(config().testPoints())
    , fOutputCat(config().outputCategory())
    {}
  
  /// Set up the algorithm.
  void Setup(spacecharge::SpaceCharge const& sc) { pSC = &sc; }
  
  /// Print an introduction.
  void Intro() const;
  
  /// Runs the test and returns the number of errors.
  unsigned int Run();
  
  /// Print a summary.
  void Summary() {}
  
  
  template <typename Stream>
  void PrintPoint(Stream&& out, geo::Point_t const& point);
  
}; // class SCEtestAlgo


void SCEtestAlgo::Intro() const {
  mf::LogVerbatim(fOutputCat)
    <<   "Space charge displacement correction: "
      << (pSC->EnableSimSpatialSCE()? "enabled": "disabled")
    << "\nSpace charge field correction:        "
      << (pSC->EnableSimEfieldSCE()? "enabled": "disabled")
    ;
} // SCEtestAlgo::Intro()


unsigned int SCEtestAlgo::Run() {
  
  mf::LogVerbatim log(fOutputCat);
  for (geo::Point_t const& point: fTestPoints) PrintPoint(log, point);
  return 0;
} // SCEtestAlgo::Run()


template <typename Stream>
void SCEtestAlgo::PrintPoint(Stream&& out, geo::Point_t const& point) {
  
  auto const displacement = pSC->GetPosOffsets(point);
  auto const EfieldDistortion = pSC->GetEfieldOffsets(point);
  
  out
    << "\nPoint: " << point
    << "  displacement: " << displacement
    << "  E field distortion: " << EfieldDistortion
    ;
  
} // SCEtestAlgo::PrintPoint()


//------------------------------------------------------------------------------
//---  The test environment
//---

/*
 * GeometryTesterEnvironment, configured with a geometry-aware configuration
 * object, is used in a non-Boost-unit-test context.
 * After proper setup it will provide:
 * - `geo::GeometryCore`
 * - all the other services configured as dependencies
 * 
 * BUG: again, MicroBooNE geometry requires a custom mapping: using
 * `geo::ChannelMapStandardAlg` may (and should be expected to) yield wrong
 * results.
 */
using TesterConfiguration
  = testing::BasicGeometryEnvironmentConfiguration<geo::ChannelMapStandardAlg>;
using TestEnvironment = testing::GeometryTesterEnvironment<TesterConfiguration>;


//------------------------------------------------------------------------------
//---  The tests
//---


/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 * 
 * The arguments in argv are:
 * 0. name of the executable ("SpaceChargeMicroBooNE_test")
 * 1. (mandatory) path to the FHiCL configuration file
 * 2. (default: `scetest`) path to the FHiCL configuration file 
 * 3. (default: `services.SpaceChargeService`) FHiCL path to the configuration
 *    of SpaceCharge service
 * 
 * The configuration of the test consist of an element "testPoints" with
 * a list of 3-element lists; each one describes a test point. For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * physics: {
 *   analyzers: {
 *     scetest: {
 *       testPoints: [
 *           [ -100.0, -100.0, -100.0 ],
 *           [    0.0,    0.0,    0.0 ],
 *           [  100.0,  100.0,  100.0 ]
 *         ] # testPoints
 *     }
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Coordinate values are in the usual centimeters.
 * 
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv) {
  
  // ***
  // *** test environment setup
  // ***
  TesterConfiguration config("spacecharge_uboone_test");
  config.SetMainTesterParameterSetName("scetest");
  
  //
  // parameter parsing
  //
  int iParam = 0;
  
  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);
  else {
    std::cerr << "FHiCL configuration file path required as first argument!"
      << std::endl;
    return 1;
  }
  
  // second argument: path of the parameter set for geometry test configuration
  // (optional; default does not have any tester)
  if (++iParam < argc) config.SetMainTesterParameterSetPath(argv[iParam]);
  
  // third argument: path of the parameter set for DetectorClocks configuration
  // (optional; default: "services.DetectorClocks" from the inherited object)
  if (++iParam < argc)
    config.SetServiceParameterSetPath("SpaceChargeService", argv[iParam]);
  
  
  //
  // testing environment setup
  //
  TestEnvironment testEnv(config);
  
  // SpaceChargeMicroBooNE requires `DetectorProperties`.
  // Note that we must choose a specific implementation of each of those service
  // providers, e.g. `detinfo::DetectorPropertiesStandard` rather than just
  // `detinfo::DetectorProperties`. We go with the current MicroBooNE defaults.
  // DetectorPropertiesStandard and all its dependencies support the simple set
  // up (see testing::TesterEnvironment::SimpleProviderSetup()), except for
  // Geometry, that has been configured already in the geometry-aware
  // environment. So we invoke a simple set up for each of the dependencies
  testEnv.SimpleProviderSetup<detinfo::LArPropertiesStandard>();
  testEnv.SimpleProviderSetup<detinfo::DetectorClocksStandard>();
  testEnv.SimpleProviderSetup<detinfo::DetectorPropertiesStandard>();
  
  // SpaceChargeMicroBooNE supports the simple set up; so we invoke it
  testEnv.SimpleProviderSetup<spacecharge::SpaceChargeMicroBooNE>();
  
  SCEtestAlgo tester(testEnv.TesterParameters());
  tester.Setup(*(testEnv.Provider<spacecharge::SpaceCharge>()));
  
  // ***
  // *** "test"
  // ***
  
  tester.Intro();
  unsigned int nErrors = tester.Run();
  tester.Summary();
  
  return nErrors;
} // main()
