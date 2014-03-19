# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic)
#set( UBOONECODE_LIB -L$ENV{UBOONECODE_LIB} )

#link_libraries( ${LIB_NAME} ${ROOTLIB} ${LARLIGHT_LIB}) # All execs need this.

#art_make_library(LIBRARY_NAME ${LIB_NAME} SOURCE ${LIB_SOURCES} LIBRARIES DataScanner_module)
art_make( BASENAME_ONLY
	  LIBRARY_NAME     MCShowerReco
          LIB_LIBRARIES    Geometry
			   Geometry_service
			   Simulation
			   SimulationBase
#			   BackTracker_service
			   Utilities
#			   Filters
			   RawData
#			   RecoAlg
			   RecoBase
			   RecoObjects
			   LArProperties_service
			   DetectorProperties_service
			   SignalShapingServiceMicroBooNE_service
			   ${UBOONECODE_LIB}
			   ${SIMULATIONBASE}
			   ${LARRECO_LIB}
			   ${LARDATA_LIB}
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
			   ${ART_PERSISTENCY_PROVENANCE}
			   ${ART_UTILITIES}
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
			   ${FHICLCPP}
			   ${CETLIB}
			   ${ROOT_GEOM}
			   ${ROOT_XMLIO}
			   ${ROOT_GDML}
			   ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES MCShowerReco
			   Geometry
			   Geometry_service
			   Simulation
			   SimulationBase
#			   BackTracker_service
			   Utilities
#			   Filters
			   RawData
#			   RecoAlg
			   RecoBase
			   RecoObjects
			   LArProperties_service
			   DetectorProperties_service
			   SignalShapingServiceMicroBooNE_service
			   ${UBOONECODE_LIB}
			   ${SIMULATIONBASE}
			   ${LARRECO_LIB}
			   ${LARDATA_LIB}
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
			   ${ART_PERSISTENCY_PROVENANCE}
			   ${ART_UTILITIES}
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
			   ${FHICLCPP}
			   ${CETLIB}
			   ${ROOT_GEOM}
			   ${ROOT_XMLIO}
			   ${ROOT_GDML}
			   ${ROOT_BASIC_LIB_LIST}
			)

# install_headers()
install_fhicl()
install_source()
