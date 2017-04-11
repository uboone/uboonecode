
art_make(
	  LIB_LIBRARIES
	                uboone_SNStreamSim_Fmwk
			#${MF_MESSAGELOGGER}
			#${MF_UTILITIES}
			#cetlib cetlib_except
                        #${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
		        cetlib cetlib_except
        )


install_headers()
install_fhicl()
install_source()

