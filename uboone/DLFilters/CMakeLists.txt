
art_make( BASENAME_ONLY
	  LIBRARY_NAME  uboone_DLFilters
	  LIB_LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES 
			larcorealg_Geometry
	  	        #larcore_Geometry_Geometry_service
		        lardata_Utilities
	                uboone_Geometry
			lardataobj_RecoBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
			nusimdata_SimulationBase
        )							   

install_headers()
install_source()
install_fhicl()
#add_subdirectory(job)
