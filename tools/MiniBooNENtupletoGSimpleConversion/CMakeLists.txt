# project name
project(ConvertNtuple)

include_directories( $ENV{ROOT_INC}
		     $ENV{GENIE_INC}/GENIE )

                       
cet_make_exec(ConvertNtuple SOURCE ConvertNtuple.cxx BooNEtoGSimple.cxx 
			    LIBRARIES ${ROOT_BASIC_LIB_LIST} 
			    	      xml2
				      ${GFWPARDAT}
				      ${LOG4CPP} 
				      ${GFWALG}
				      ${GFWPARDAT}
				      ${GPHDISXS}
				      ${GFWEG}
				      ${GFWEG}
	                              ${GTLFLX}
				      ${GFWGHEP}
				      ${GFWINT}
				      ${GFWMSG}
				      ${GPHNUCLST}
				      ${GFWNUM}
				      ${GFWPARDAT}
				      ${GFWREG}
				      ${GFWUTL}
				      #more root stuff
				      ${ROOT_EG}
				      ${ROOT_EGPYTHIA6}
				      )

