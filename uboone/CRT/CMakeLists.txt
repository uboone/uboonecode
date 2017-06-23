# canvas_utlities
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( $ENV{NUTOOLS_INC} )

add_subdirectory(CRTProducts)

art_make( BASENAME_ONLY
    EXCLUDE 
        CRTChannelMapAlg.cc 
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTDaqSim_module.cc
	CRTCalibration_module.cc
	CRTCoincidence_module.cc
	CRTHitProducer_module.cc
	CRTTrackProducer_module.cc
	TimeExtraction_module.cc
	CRTMerger_module.cc
    LIBRARY_NAME  
        CRTData
    LIB_LIBRARIES
        larcore_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common 
        canvas
        art_Persistency_Provenance 
        canvas
        art_Utilities 
        canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)

art_make_library( LIBRARY_NAME CRT
    SOURCE
    	CRTGeoObjectSorter.cc
    	CRTChannelMapAlg.cc
    	CRTAuxFunctions.cc
    LIBRARIES
	CRTProducts
	larcore_Geometry
        CRTData
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)

#simple_plugin( CRTGeometryHelper service
#        larcore_Geometry
#        CRT
#        ${ART_FRAMEWORK_CORE}
#        ${ART_FRAMEWORK_IO_SOURCES}
#        ${ART_FRAMEWORK_PRINCIPAL}
#        art_Persistency_Common 
#        canvas
#        art_Persistency_Provenance 
#        canvas
#        art_Utilities 
#        canvas
#        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#        ${MF_MESSAGELOGGER}
#        ${MF_UTILITIES}
#        ${ROOT_BASIC_LIB_LIST}
#        cetlib cetlib_except
#    BASENAME_ONLY
#)

simple_plugin( CRTDetSim module
    CRT
    CRTProducts
    art_Utilities 
    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
    nutools_RandomUtils_NuRandomService_service
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${MF_MESSAGELOGGER}
    ${ROOT_BASIC_LIB_LIST}
    ${CLHEP}
    ${CETLIB}
  BASENAME_ONLY
)

#simple_plugin( CRTDaqSim module
#    larcore_Geometry
#    CRTProducts
#    CRT
#    art_Utilities 
#    canvas_Utilities
#    ${ART_FRAMEWORK_SERVICES_REGISTRY}
#    ${MF_MESSAGELOGGER}
#    ${ROOT_BASIC_LIB_LIST}
#    ${CLHEP}
#    ${CETLIB}
#    bernfebdaq-core_Overlays
#    ${ROOT_BASIC_LIB_LIST}
#    ${ART_UTILITIES}
#    ${ARTDAQ_CORE}
#    ${FHICLCPP}
#    ${MF_MESSAGELOGGER}
#    ${CETLIB}
#    ${GALLERY}  
#  BASENAME_ONLY
#)

simple_plugin( CRTGeometryHelper service
    larcore_Geometry
    CRT
    art_Utilities 
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${ROOT_BASIC_LIB_LIST}
    ${CETLIB}
  BASENAME_ONLY
)

simple_plugin(CRTCalibration "module"
        CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTCoincidence "module"
            CRT
              bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )


simple_plugin(CRTHitProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTTrackProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(TimeExtraction "module"
            CRT
              bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

#simple_plugin( CRTMerger module
#  bernfebdaq-core_Overlays
#  ${ROOT_BASIC_LIB_LIST}
#  ${ART_UTILITIES}
#  ${ART_FRAMEWORK_SERVICES_REGISTRY}
#  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
#  ${ARTDAQ_CORE}
#  ${FHICLCPP}
#  ${MF_MESSAGELOGGER}
#  ${MF_UTILITIES}
#  ${CETLIB}
#  ${GALLERY}
#  CRTProducts
#=======
#		CRTData
#        ${ART_FRAMEWORK_CORE}
#        ${ART_FRAMEWORK_IO_SOURCES}
#        ${ART_FRAMEWORK_PRINCIPAL}
#        art_Persistency_Common 
#        canvas
#        art_Persistency_Provenance 
#        canvas
#        art_Utilities 
#        canvas
#        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
#        nutools_RandomUtils_NuRandomService_service
#        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#        ${MF_MESSAGELOGGER}
#       ${MF_UTILITIES}
#        ${ROOT_BASIC_LIB_LIST}
#        ${CLHEP}
#        cetlib cetlib_except
#    BASENAME_ONLY
#>>>>>>> v06_36_00
#)


install_headers()
install_fhicl()
install_source()
