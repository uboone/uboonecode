# project name
project(ConvertNtuple)

include_directories( $ENV{ROOT_INC}
		     $ENV{GENIE_INC}/GENIE )

find_library(XML2LIB NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib)

cet_make_exec(ConvertNtuple SOURCE ConvertNtuple.cxx BooNEtoGSimple.cxx 
			    LIBRARIES ${ROOT_BASIC_LIB_LIST} 
			    	      ${XML2LIB}
				      ${GPDG}
				      ${LOG4CPP} 
				      ${GALGORITHM}
				      ${GBARYONRESONANCE}
				      ${GBASE}
				      ${GEVGCORE}
				      ${GEVGDRIVERS}
	                              ${GFLUXDRIVERS}
				      ${GHEP}
				      ${GINTERACTION}
				      ${GMESSENGER}
				      ${GNUCLEAR}
				      ${GNUMERICAL}
				      ${GPDG}
				      ${GREGISTRY}
				      ${GUTILS}
				      #more root stuff
				      ${ROOT_EG}
				      ${ROOT_EGPYTHIA6}
				      )

