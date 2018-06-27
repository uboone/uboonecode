#include_directories ( $ENV{GENIE_INC}/GENIE )

art_make(LIBRARY_NAME  uboone_EventWeight_Calculators_BNBPrimaryHadron
         LIB_LIBRARIES
                       larsim_EventWeight_Base
                       lardata_Utilities
                       cetlib cetlib_except
                       ${ART_FRAMEWORK_CORE}
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                       ${CLHEP}
                       ${ROOT_BASIC_LIB_LIST}
                       ${Boost_SYSTEM_LIBRARY}
                       nutools_NuReweight_art   
                       nutools_NuReweight       
                       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
                       nusimdata_SimulationBase
                       art_Persistency_Common canvas canvas
                       art_Persistency_Provenance canvas canvas
                       art_Utilities canvas canvas
                       ${FHICLCPP}
        )

install_headers()
install_fhicl()
install_source()


