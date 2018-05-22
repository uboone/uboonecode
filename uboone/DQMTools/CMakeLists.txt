
art_make( 
          LIB_LIBRARIES
			lardataobj_RawData
			lardataobj_RecoBase 
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES
     uboone_DQMTools
	    larreco_HitFinder
                        larsim_MCCheater_ParticleInventoryService_service
                        larsim_MCCheater_BackTrackerService_service
			larreco_RecoAlg
                        lardataobj_RecoBase
			lardataobj_MCBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
