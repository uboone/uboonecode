# ======================================================================
#  uboonecode main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardataobj/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(uboonecode VERSION 10.09.00 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package( ubobj REQUIRED EXPORT )
find_package( lardata REQUIRED EXPORT )
find_package( nurandom REQUIRED EXPORT )
find_package( nugen REQUIRED EXPORT )
find_package( GENIE REQUIRED EXPORT )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(fcl)
add_subdirectory(test)
add_subdirectory(tools)
add_subdirectory(releaseDB)

# ups - table and config files

# packaging utility
cet_cmake_config() ### MIGRATE-ADDED (migrate-3.21.01) - required cmd
