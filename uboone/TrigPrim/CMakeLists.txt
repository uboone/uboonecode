
art_make( 
          LIB_LIBRARIES    lardataobj_RecoBase
	  		   larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larevt_Filters
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
			   larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
			   nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
   			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			   art_Persistency_Common canvas
		       	   art_Persistency_Provenance canvas
			   art_Utilities canvas
			   ${MF_MESSAGELOGGER}
			   
	                   ${FHICLCPP}
			   cetlib cetlib_except
                           ${CLHEP}
			   ${Boost_SYSTEM_LIBRARY}
			   ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_GEOM}


			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			lardataobj_RawData
			lardataobj_RecoBase 
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}



         MODULE_LIBRARIES 
                           lardataobj_RecoBase
                           lardata_ArtDataHelper
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
			   larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                           lardata_Utilities
			    ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas
			   art_Utilities canvas
			   ${MF_MESSAGELOGGER}
			   
	                   ${FHICLCPP}
			   cetlib cetlib_except
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_GEOM}

			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}	
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}	
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			${CETLIB}
			larreco_HitFinder
                        larsim_MCCheater_BackTrackerService_service
	                larsim_MCCheater_ParticleInventoryService_service
			lardataobj_MCBase
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Framework_Services_System_TriggerNamesService_service


         )

install_headers()
install_fhicl()
install_source()
