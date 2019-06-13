# project name
project(ConvertNtuple)

include_directories( $ENV{ROOT_INC}
		     $ENV{GENIE_INC}/GENIE )

                       
cet_make_exec(ConvertNtuple SOURCE ConvertNtuple.cxx BooNEtoGSimple.cxx 
			    LIBRARIES ${ROOT_BASIC_LIB_LIST} 
			    	      xml2
				      ${LOG4CPP} 
				      ${GENIE_LIB_LIST}
                                      ${GSL}
                                      ${GSLCBLAS}
				      #more root stuff
				      ${ROOT_EG}
				      ${ROOT_EGPYTHIA6}
				      )

