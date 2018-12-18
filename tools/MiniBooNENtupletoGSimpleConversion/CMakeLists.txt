# project name
project(ConvertNtuple)

include_directories( $ENV{ROOT_INC} )
include_directories( $ENV{GENIE_INC}/GENIE )

cet_make_exec(
  ConvertNtuple
  SOURCE
    ConvertNtuple.cxx BooNEtoGSimple.cxx
  LIBRARIES
    xml2
    ${LOG4CPP}
    ${GENIE_LIB_LIST}
    ${GSL}
    ${GSLCBLAS}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_GENVECTOR}
    ${ROOT_EGPYTHIA6}
    ${ROOT_EG}
)

