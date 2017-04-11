file(GLOB WGHCALCFILES *WeightCalc* Factory* )
include_directories ( $ENV{GENIE_INC}/GENIE )

art_make_library(LIBRARY_NAME WeightCalculators 
                 SOURCE	      ${WGHCALCFILES}
                 LIBRARIES
		       lardata_Utilities
                       nutools_NuReweight_art    # Uncomment this to use GENIE weights
                       nutools_NuReweight        # Uncomment this to use GENIE weights
		       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
		       nusimdata_SimulationBase
              	      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
		      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		      art_Persistency_Common canvas canvas
                      art_Persistency_Provenance canvas canvas
                      art_Utilities canvas canvas
                      cetlib cetlib_except
		      ${FHICLCPP}
                      ${MF_MESSAGELOGGER}
		      ${CLHEP}
                      ${ROOT_BASIC_LIB_LIST}
                      ${Boost_SYSTEM_LIBRARY}
	        ) 

art_make(LIBRARY_NAME   EventWeight
         EXCLUDE        
                      ${WGHCALCFILES}
                      #GenieWeightCalc.cxx        # Comment this to use GENIE weights
	 LIB_LIBRARIES
                      BNBPrimaryHadronWeightCalculators
		      WeightCalculators
                      ${MF_MESSAGELOGGER}
                      ${MF_UTILITIES}
                      ${Boost_SYSTEM_LIBRARY}
        MODULE_LIBRARIES EventWeight
			 WeightCalculators
                         BNBPrimaryHadronWeightCalculators

	)


install_headers()
install_fhicl()
install_source()
add_subdirectory(BNBPrimaryHadron)
add_subdirectory(jobs)
