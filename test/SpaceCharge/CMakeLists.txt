cet_test( SpaceChargeMicroBooNE_test
  LIBRARIES
    uboone_SpaceCharge
    mf_MessageLogger
    mf_Utilities
    fhiclcpp
    ${ROOT_GENVECTOR}
    ${GSL} ${GSLCBLAS}
  DATAFILES spacechargetest_uboone.fcl
  TEST_ARGS ./spacechargetest_uboone.fcl
)
