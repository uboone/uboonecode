cet_build_plugin(
  InteractTest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  lardata::RecoObjects
)

cet_build_plugin(
  KalmanFilterTest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  lardata::RecoObjects
  larcore::Geometry_Geometry_service
)

cet_build_plugin(
  PropTest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  lardata::RecoObjects
  nurandom::RandomUtils_NuRandomService_service
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
