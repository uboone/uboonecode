

cet_find_library( BOOST_SYSTEM    NAMES boost_system     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_parse_args(AE "LIBRARIES" "" ${ARGN})
art_make_exec(NAME argoLinked
              SOURCE argo-backend.cc
ColorMap.cpp          MakePng.cpp           SocketServer.cpp      cencode.cpp
FormString.cpp        RecordComposer.cpp    Timer.cpp
FormulaLooter.cpp     ResultComposer.cpp    TreeElementLooter.cpp
JsonElement.cpp       RootToJson.cpp        TreeReader.cpp
              ${AE_DEFAULT_ARGS}
              LIBRARIES
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
              ${ROOT_TREE}
              ${ROOT_TREEPLAYER}
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
              png
              )


install_headers()
install_source()
