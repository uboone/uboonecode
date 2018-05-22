if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage -Wno-undefined-var-template )
endif()

# canvas_utlities
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

add_subdirectory(CRTProducts)
add_subdirectory(CRTBernFEBDAQCore)
#add_subdirectory(CRTBinarySource)

art_make( BASENAME_ONLY
    EXCLUDE 
        CRTChannelMapAlg.cc 
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTDaqSim_module.cc
        CRTHitProducer_module.cc
        CRTCalibration_module.cc
        CRTTrackProducer_module.cc
	FlashExt_module.cc	
	CRTAuxFunctions.cc
  	CRTMerger_module.cc
	TimeExtraction_module.cc
	MergeTest_module.cc
	CRTRawtoCRTHit_source.cc
	CRTRawtoCRTHitGPScorr_source.cc
	CRTAnalyser_module.cc	
    LIBRARY_NAME  
        CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
	lardataobj_RecoBase
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)

art_make_library( 
    LIBRARY_NAME 
    	CRT
    SOURCE
        CRTChannelMapAlg.cc 
	CRTAuxFunctions.cc	
    LIBRARIES 
    	#larcorealg_Geometry
        CRTData
        ${MF_MESSAGELOGGER}
        
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        cetlib cetlib_except
)

simple_plugin( CRTGeometryHelper service
        larcorealg_Geometry
        CRT
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        cetlib cetlib_except
)

simple_plugin( CRTDetSim module
        CRT
	CRTData
	CRTProducts
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
    BASENAME_ONLY
)

simple_plugin(CRTCalibration "module"
        CRT
	CRTData
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common canvas
        art_Persistency_Provenance canvas
        art_Utilities canvas
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
      nutools_RandomUtils_NuRandomService_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
        ${MF_MESSAGELOGGER}
        
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        cetlib cetlib_except
    BASENAME_ONLY
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
          ${ROOT_BASIC_LIB_LIST} 
	  cetlib cetlib_except
)

simple_plugin(CRTTrackProducer "module"
            CRT
              bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} 
	  cetlib cetlib_except
)

simple_plugin(FlashExt "module"
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

simple_plugin(MergeTest "module"
            CRT
	  bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTRawtoCRTHit "source"
            CRT
	  bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTRawtoCRTHitGPScorr "source"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin(CRTAnalyser "module"
            CRT
          bernfebdaq-core_Overlays
          art_Framework_Services_Optional
          art_Framework_Services_Registry
          art_Framework_Services_Optional_TFileService_service
          ${ROOT_BASIC_LIB_LIST} )

simple_plugin( CRTMerger "module"
  CRT
  bernfebdaq-core_Overlays
  ${CETLIB}	
  ${ROOT_BASIC_LIB_LIST}
  ${ART_UTILITIES}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
  ${ARTDAQ_CORE}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  
  ${CETLIB}
  ${GALLERY}
  CRTProducts
  ifdh
  gallery
  boost_date_time
  cetlib cetlib_except	
)


install_headers()
install_fhicl()
install_source()


