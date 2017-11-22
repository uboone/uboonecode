include_directories ( $ENV{GENIE_INC}/GENIE )
file(GLOB INTERFACEFILES *Interface* )
simple_plugin( FluxReader
               "source"
	       FluxReader
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
               BASENAME_ONLY
             )

art_make_library(LIBRARY_NAME FluxInterface 
                 SOURCE	      ${INTERFACEFILES}
                 LIBRARIES
		       nusimdata_SimulationBase
#GENIE STUFF########################################
			xml2
			${LOG4CPP}
			${GALGORITHM}
			${GBASE}
			${GEVGCORE}
			${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
			${GINTERACTION}
			${GMESSENGER}
			${GHEP}
####################################################
                      cetlib cetlib_except
		      ${FHICLCPP}
                      #${MF_MESSAGELOGGER}
		      ${CLHEP}
                      ${ROOT_BASIC_LIB_LIST}
                      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make( EXCLUDE FluxReader_source.cc ${INTERFACEFILES}
	  BASENAME_ONLY
          LIBRARY_NAME FluxReader
          LIB_LIBRARIES	nusimdata_SimulationBase
	  		FluxInterface
			larcoreobj_SummaryData
	  		${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
			${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}

          MODULE_LIBRARIES FluxReader
	  		FluxInterface
                        nusimdata_SimulationBase
			larcoreobj_SummaryData
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
#GENIE STUFF########################################
			xml2
			${LOG4CPP}
			${GALGORITHM}
			${GBASE}
			${GEVGCORE}
			${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
			${GINTERACTION}
			${GMESSENGER}
			${GHEP}
####################################################
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )                                                           

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
