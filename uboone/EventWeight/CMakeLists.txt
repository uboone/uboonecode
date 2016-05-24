file(GLOB WGHCALCFILES *WeightCalc* Factory*)
include_directories ( $ENV{GENIE_INC}/GENIE )

art_make_library(LIBRARY_NAME WeightCalculators 
                 SOURCE	      ${WGHCALCFILES}
	         LIBRARIES
		       lardata_Utilities
		       NuReweightArt	
                       NuReweight
		       larsim_Simulation
		       ${SIMULATIONBASE}
              	      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                      ${ART_PERSISTENCY_COMMON}
                      ${ART_PERSISTENCY_PROVENANCE}
#
                      ${ART_UTILITIES}
		      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                      ${CETLIB}
		      ${FHICLCPP}
		      ${CLHEP}
                      ${ROOT_BASIC_LIB_LIST}
		      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make(LIBRARY_NAME   EventWeight
         EXCLUDE        ${WGHCALCFILES}
	 LIB_LIBRARIES
		      WeightCalculators
#              	      ${ART_FRAMEWORK_CORE}
#		      ${ART_FRAMEWORK_IO_SOURCES}
#                      ${ART_FRAMEWORK_PRINCIPAL}
#                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                      ${ART_PERSISTENCY_COMMON}
#                      ${ART_PERSISTENCY_PROVENANCE}
#                      ${ART_UTILITIES}
		      larsim_RandomUtils_LArSeedService_service
                      ${MF_MESSAGELOGGER}
                      ${MF_UTILITIES}
#                      ${FHICLCPP}
#                      ${CETLIB}
#                      ${ROOT_GEOM}
#                      ${ROOT_XMLIO}
#                      ${ROOT_GDML}
#                      ${ROOT_BASIC_LIB_LIST}
        MODULE_LIBRARIES EventWeight
			 WeightCalculators
	)


install_headers()
install_fhicl()
install_source()
