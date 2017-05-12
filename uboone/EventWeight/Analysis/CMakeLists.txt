art_make(LIBRARY EventWeightAnalysis
         LIB_LIBRARIES
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${Boost_SYSTEM_LIBRARY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              lardata_Utilities
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${CLHEP}
              ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES
              ${ROOT_BASIC_LIB_LIST}
              lardata_Utilities
              larsim_Simulation
              nutools_RandomUtils_NuRandomService_service
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              art_Persistency_Common canvas canvas
              art_Persistency_Provenance canvas canvas
              art_Utilities canvas canvas
              cetlib cetlib_except
              ${CETLIB}
              ${FHICLCPP}
              ${CLHEP}
              ${Boost_SYSTEM_LIBRARY}
        ) 

install_headers()
install_fhicl()
install_source()

