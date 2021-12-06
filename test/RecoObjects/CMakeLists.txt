simple_plugin(KalmanFilterTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      art::Framework_Core
	      art::Framework_Principal
	      art::Framework_Services_Registry
              art::Persistency_Common
              art::Persistency_Provenance
	      art_root_io::tfile_support ROOT::Core
	      art_root_io::TFileService_service
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
              art::Utilities
              canvas::canvas
	      messagefacility::MF_MessageLogger
	      fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
	      PostgreSQL::PostgreSQL
              BASENAME_ONLY
             )

simple_plugin(PropTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      nurandom::RandomUtils_NuRandomService_service
	      art::Framework_Core
	      art::Framework_Services_Registry
	      art::Framework_Services_Optional_RandomNumberGenerator_service
	      messagefacility::MF_MessageLogger
	      CLHEP::CLHEP
              art::Persistency_Provenance
	      fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
	      PostgreSQL::PostgreSQL
              BASENAME_ONLY
             )

simple_plugin(InteractTest "module"
              lardata_RecoObjects
	      larcorealg_Geometry
	      art::Framework_Core
	      art::Framework_Principal
	      art::Framework_Services_Registry
              art::Persistency_Common
              art::Persistency_Provenance
	      art_root_io::tfile_support ROOT::Core
	      art_root_io::TFileService_service
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art::Utilities canvas::canvas
	      messagefacility::MF_MessageLogger
	      fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
	      PostgreSQL::PostgreSQL
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
