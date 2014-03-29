
art_make( BASENAME_ONLY
          LIBRARY_NAME  Looter 
          
          LIB_LIBRARIES     Geometry
                            Geometry_service
                      			Simulation
                      			BackTracker_service
                      	    Utilities
                      			Filters
                      			RawData
                      			RecoBase
                      			RecoObjects   
                      			SignalShapingServiceMicroBooNE_service
                            art_Framework_IO 
                            art_Framework_IO_Root
                            art_Framework_IO_Root_RootInput_source
                            art_Framework_Services_System_FileCatalogMetadata_service
                            art_Utilities ${CETLIB} art_Framework_Art
                            ${SIMULATIONBASE}
                            ${LARRECO_LIB}
                            ${LARDATA_LIB} 
                            ${ART_FRAMEWORK_CORE}
                      			${ART_FRAMEWORK_PRINCIPAL}
                      			${ART_FRAMEWORK_SERVICES_REGISTRY}
                      			${ART_FRAMEWORK_SERVICES_OPTIONAL}
                      			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                      			${ART_IO_ROOT}
                      			${ART_PERSISTENCY_COMMON}
                      			${ART_PERSISTENCY_PROVENANCE}
                      			${ART_UTILITIES}
                      			${MF_MESSAGELOGGER}
                      			${MF_UTILITIES}
	                          ${FHICLCPP}
			                       ${CETLIB}
                             ${ART_EXEC_ROOT_LIBS}
                             
			                       ${ROOT_BASIC_LIB_LIST}
                             ${ROOT_TREE}
                             ${ROOT_GEOM}
                             ${ROOT_XMLIO}
                             ${ROOT_GDML}
            MODULE_LIBRARIES Looter
            SERVICE_LIBRARIES Looter
 			                       ${ROOT_BASIC_LIB_LIST}
                              

        )

# install_headers()
install_fhicl()
install_source()

