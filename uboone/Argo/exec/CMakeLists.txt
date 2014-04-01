  add_subdirectory(Looter)
  add_subdirectory(ArgoEventProcessor)
  
# Build an art exec.
# macro(art_exec TARGET_STEM IN_STEM MAIN_FUNC)
#   cet_parse_args(AE "LIBRARIES" "" ${ARGN})
#   set(ART_MAIN_FUNC ${MAIN_FUNC})
#   # configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${IN_STEM}.cc.in
#   #   ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_STEM}.cc @ONLY
#   #   )
#   art_make_exec(NAME ${TARGET_STEM}
#     SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_STEM}.cc
#     ${AE_DEFAULT_ARGS}
#     LIBRARIES
#     ${AE_LIBRARIES}
#     art_Framework_Art
#     # ROOT libraries added for use convenience ONLY -- they will go away
#     # eventually.
#     ${ART_EXEC_ROOT_LIBS}
#     ${MF_MESSAGELOGGER}
#     )
# endmacro()

cet_find_library( BOOST_SYSTEM    NAMES boost_system     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

cet_parse_args(AE "LIBRARIES" "" ${ARGN})
art_make_exec(NAME argo 
              SOURCE argo.cc
              ${AE_DEFAULT_ARGS}
              LIBRARIES
              ArgoEventProcessor
              ${AE_LIBRARIES}
              # art_Version
              art_Framework_IO
              art_Framework_IO_Root
              art_Framework_IO_Root_RootInput_source
              art_Framework_Services_System_FileCatalogMetadata_service
              art_Framework_Services_System_CurrentModule_service
              art_Framework_Services_Registry
              art_Framework_Services_FileServiceInterfaces
              art_Utilities ${CETLIB} art_Framework_Art
              
              ${ART_EXEC_ROOT_LIBS}
              ${MF_MESSAGELOGGER}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${BOOST_SYSTEM}
              RecoBase
              RecoBase_dict
              RecoBase_map
              Geometry
              Geometry_service
              LArProperties_service
              DetectorProperties_service
              Utilities
              Simulation
              )

# Standard execs
# art_exec(argo art artapp)

# Execs with Boost unit testing enabled for modules.
#art_exec(argo_ut argo_ut artapp LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# art_make_exec(NAME check_libs
#   SOURCE ${CMAKE_CURRENT_BINARY_DIR}/check_libs.cc
#   LIBRARIES
#   art_Utilities
#   ${CETLIB}
#   )

install_headers()
install_source()
