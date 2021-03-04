# project name
project(ConvertNtuple)

include_directories( $ENV{ROOT_INC} )

find_library(XML2LIB NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib)

cet_make_exec(ConvertNtuple SOURCE ConvertNtuple.cxx BooNEtoGSimple.cxx 
			    LIBRARIES ${ROOT_BASIC_LIB_LIST} 
			    	      ${XML2LIB}
				      ${GFWPARDAT}
				      ${LOG4CPP} 
				      ${GFWALG}
				      ${GFWPARDAT}
				      ${GPHDISXS}
                                     ${GPHPDF}
                                     ${GPHXSIG}
                                     ${GPHDISEG}
                                     ${GPHHADNZ}
                                     ${GPHDCY}
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
                                     ${GPHCMN}
                                      ${GFWUTL}
                                     ${GSL}
                                     ${GSLCBLAS}
				      #more root stuff
				      ${ROOT_EG}
				      ${ROOT_EGPYTHIA6}
				      ${ROOT_GENVECTOR}
				      )
