
art_make(
	  LIB_LIBRARIES
			#${MF_MESSAGELOGGER}
			#${MF_UTILITIES}
			#cetlib cetlib_except
                        #${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
		        ${ART_FRAMEWORK_CORE}	
		        ${FHICLCPP}		
        )


install_headers()
install_fhicl()
install_source()

