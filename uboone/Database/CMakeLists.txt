art_make( BASENAME_ONLY
          LIBRARY_NAME      UbooneDBI_Providers
          LIB_LIBRARIES     CalibrationDBI_Providers
	  		    CalibrationDBI_IOVData
			    Geometry
			    Geometry_service
			    ${FHICLCPP}
			    ${CETLIB}
			    ${Boost_SYSTEM_LIBRARY}
			    ${ROOT_BASIC_LIB_LIST}
			    ${ART_UTILITIES}
	  SERVICE_LIBRARIES CalibrationDBI_Providers
	                    CalibrationDBI_IOVData
			    UbooneDBI_Providers
			    RawData
			    Geometry
			    Geometry_service
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                            ${ART_PERSISTENCY_COMMON}
                            ${ART_PERSISTENCY_PROVENANCE}
			    ${ART_UTILITIES}
			    ${FHICLCPP}
			    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

