# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

set( SWTRIGGER_LIB -L$ENV{SWTRIGGER_LIBDIR} -lSWTriggerBase -lFEMBeamTrigger )

include_directories ( $ENV{SWTRIGGER_INCDIR} )

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} ${SWTRIGGER_LIB}) # All execs need this.

add_subdirectory(fcl)

art_make( LIBRARY_NAME SPEcalibrationCode
	  EXCLUDE SPEcalibration_module.cc CalibrationFit.C
	  RawData
	  UBOpReadoutMap_service
	  LIB_LIBRARIES ${ART_FRAMEWORK_CORE} 
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
	       ${ROOT_BASIC_LIB_LIST}	     
	       BASENAME_ONLY  
          MODULE_LIBRARIES 	  		   
	  	   SPEcalibrationCode
		   uboone_RawData
		   Geometry
		   uboone_Geometry
		   Geometry_service
		   UBOpReadoutMap_service
		   TimeService_service
		   Simulation
		   SimulationBase
		   Utilities
		   MCBase
		   SummaryData
		   RawData
		   RecoBase
		   RecoObjects
		   LArProperties_service
		   DetectorProperties_service
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
		

simple_plugin( SPEcalibration "module"
	       SPEcalibrationCode
	       OpticalSubEvents
	       UBOpReadoutMap_service
	       RawData
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

#simple_plugin( FEMemulator "module"
#	       SPEcalibrationCode
#	       OpticalSubEvents
#	       UBOpReadoutMap_service
#	       RawData
#	       ${ART_FRAMEWORK_CORE}
#	       ${ART_FRAMEWORK_IO_SOURCES}
#              ${ART_FRAMEWORK_PRINCIPAL}
#              ${ART_FRAMEWORK_SERVICES_REGISTRY}
#              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#              ${ART_PERSISTENCY_COMMON}
#              ${ART_PERSISTENCY_PROVENANCE}
#              ${ART_UTILITIES}
#	       ${ROOT_BASIC_LIB_LIST}	       
#              BASENAME_ONLY
#             )


install_fhicl()
install_headers()
install_source()

