art_make(LIBRARY_NAME     uboone_EventWeight
	 LIB_LIBRARIES
                          larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          uboone_EventWeight_Calculators
                          uboone_EventWeight_Calculators_BNBPrimaryHadron
                          #${MF_MESSAGELOGGER}
                          #
                          ${Boost_SYSTEM_LIBRARY}
         MODULE_LIBRARIES larsim_EventWeight_Base
                          larsim_EventWeight_Calculators
                          uboone_EventWeight_Calculators
                          uboone_EventWeight_Calculators_BNBPrimaryHadron
			  ${MF_MESSAGELOGGER}
                     	  
	)


install_headers()
install_fhicl()
install_source()
