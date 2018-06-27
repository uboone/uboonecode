simple_plugin(SparseRawDigitTest "module"
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_BASIC}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
              art_Persistency_Common canvas
              art_Persistency_Provenance canvas
              art_Utilities canvas
              ${MF_MESSAGELOGGER}
              
              ${FHICLCPP}
              ${CETLIB}
              ${PQ}
              BASENAME_ONLY
             )

install_headers()
install_fhicl()
install_source()
