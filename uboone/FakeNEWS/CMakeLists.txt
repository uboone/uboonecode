

art_make( EXCLUDE SingleSimChannel_module.cc SingleSimChannelAna_module.cc
          BASENAME_ONLY
          LIBRARY_NAME     uboonecode_uboone_FakeNEWS
          LIB_LIBRARIES
	)

simple_plugin(SingleSimChannel "module"
              BASENAME_ONLY
	      uboonecode_uboone_FakeNEWS
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
	      lardataobj_Simulation
	      larcoreobj_SummaryData
              lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              art_Persistency_Common canvas
              art_Persistency_Provenance canvas
              art_Utilities canvas
              ${MF_MESSAGELOGGER}
              
              ${FHICLCPP}
              ${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(SingleSimChannelAna "module"
              BASENAME_ONLY
	      uboonecode_uboone_FakeNEWS
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
	      lardataobj_RecoBase
	      lardataobj_RawData
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              art_Persistency_Common canvas
              art_Persistency_Provenance canvas
              art_Utilities canvas
              ${MF_MESSAGELOGGER}
              
              ${FHICLCPP}
              ${ROOT_BASIC_LIB_LIST}
        )

install_fhicl()

