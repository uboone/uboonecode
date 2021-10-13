simple_plugin(KalmanFilterTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
              art_Persistency_Common
              art_Persistency_Provenance
	      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
	      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
              art_Utilities
              canvas
	      ${MF_MESSAGELOGGER}
	      ${FHICLCPP}
              cetlib
              cetlib_except
	      PQ
              BASENAME_ONLY
             )

simple_plugin(PropTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      nurandom_RandomUtils_NuRandomService_service
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
	      ${MF_MESSAGELOGGER}
	      CLHEP::CLHEP
              art_Persistency_Provenance
	      ${FHICLCPP}
              cetlib
              cetlib_except
	      PQ
              BASENAME_ONLY
             )

simple_plugin(InteractTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
              art_Persistency_Common
              art_Persistency_Provenance
	      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
	      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art_Utilities canvas
	      ${MF_MESSAGELOGGER}
	      ${FHICLCPP}
              cetlib
              cetlib_except
	      PQ
              BASENAME_ONLY
             )

#cet_test( PropTest HANDBUILT
#  	  TEST_EXEC lar
#	  TEST_ARGS -c PropTest.fcl
#	  DATAFILES PropTest.fcl )

cet_test( InteractTest HANDBUILT
  	  TEST_EXEC lar
	  TEST_ARGS -c InteractTest.fcl
	  DATAFILES PropTest.fcl )

cet_test( KalmanFilterTest HANDBUILT
  	  TEST_EXEC lar
	  TEST_ARGS -c KalmanFilterTest.fcl
	  DATAFILES KalmanFilterTest.fcl )

install_headers()
install_fhicl()
install_source()
