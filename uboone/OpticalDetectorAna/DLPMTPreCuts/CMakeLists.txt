# Either set these here or set environment variables

include_directories ( $ENV{DLPMTPRECUTALGO_INCDIR} )
find_library( DLPMTPRECUTALGO_LIBBASE NAMES LEEPreCutsAlgo PATHS $ENV{DLPMTPRECUTALGO_LIBDIR} )

#add_subdirectory(fcl)

simple_plugin( DLPMTPreCuts "module"
	       lardataobj_RecoBase
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
	       ${DLPMTPRECUTALGO_LIBBASE}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

install_fhicl()
install_headers()
