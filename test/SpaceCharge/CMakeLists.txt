cet_test( SpaceChargeMicroBooNE_test
  LIBRARIES
    ubevt_SpaceCharge
    lardataalg_DetectorInfo
    larcorealg_Geometry
    mf_MessageLogger
    fhiclcpp
    cetlib cetlib_except
    ${ROOT_GENVECTOR}
    ${ROOT_BASIC_LIB_LIST}
    ${GSL} ${GSLCBLAS}
  DATAFILES spacechargetest_uboone.fcl
  TEST_ARGS ./spacechargetest_uboone.fcl
)
