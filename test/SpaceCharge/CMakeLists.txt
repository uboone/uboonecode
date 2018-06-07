cet_test( SpaceChargeMicroBooNE_test
  LIBRARIES
    ubevt_SpaceCharge
    mf_MessageLogger
    fhiclcpp
    ${ROOT_GENVECTOR}
    ${GSL} ${GSLCBLAS}
  DATAFILES spacechargetest_uboone.fcl
  TEST_ARGS ./spacechargetest_uboone.fcl
)
