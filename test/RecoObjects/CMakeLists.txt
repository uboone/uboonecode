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

cet_test( InteractTest.sh PREBUILT )

cet_test( KalmanFilterTest.sh PREBUILT )

install_headers()
install_fhicl()
install_source()
