simple_plugin(GENIETweaker "module"
  BASENAME_ONLY
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${GENIE_LIB_LIST}
      ${MF_MESSAGELOGGER}
      ${FHICLCPP}
      ${LOG4CPP}
      cetlib cetlib_except
      ${CLHEP}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}
      ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(GENIEExtractor "module"
  BASENAME_ONLY
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${GENIE_LIB_LIST}
      ${MF_MESSAGELOGGER}
      ${FHICLCPP}
      ${LOG4CPP}
      cetlib cetlib_except
      ${CLHEP}
      ${ROOT_EGPYTHIA6}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}
      ${ROOT_BASIC_LIB_LIST}
      nusimdata_SimulationBase
      nutools_EventGeneratorBase_GENIE
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
)

# install_headers()
install_fhicl()
install_source()
