
art_make(
	  LIB_LIBRARIES
			#${MF_MESSAGELOGGER}
			#${MF_UTILITIES}
			#${CETLIB}
                        #${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
        )


install_headers()
install_fhicl()
install_source()

