cet_find_library( BOOST_SYSTEM              NAMES boost_system     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_PROGRAM_OPTIONS     NAMES boost_program_options     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_REGEX               NAMES boost_regex     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( FHCIL                     NAMES fhiclcpp     PATHS ENV FHICLCPP_LIB NO_DEFAULT_PATH )

# cet_parse_args(AE "LIBRARIES" "" ${ARGN})
art_make_exec(NAME argophoenix
              SOURCE argo_phoenix_main.cc
              ${AE_DEFAULT_ARGS}
              LIBRARIES
              ${AE_LIBRARIES}
              # art_Version
              art_Framework_Core              
              art_Framework_IO
              art_Framework_IO_Root
              art_Framework_IO_Root_RootInput_source
              art_Framework_Services_System_FileCatalogMetadata_service
              art_Framework_Services_System_CurrentModule_service
              art_Framework_Services_Registry
              art_Framework_Services_FileServiceInterfaces
              art_Framework_EventProcessor
              art_Utilities 
              art_Framework_Art
              ${ART_EXEC_ROOT_LIBS}
              ${MF_MESSAGELOGGER}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${BOOST_SYSTEM}
              ${BOOST_PROGRAM_OPTIONS}
              ${BOOST_REGEX}
              ${FHCIL}
              ${CETLIB} 
              )

install_headers()
install_source()