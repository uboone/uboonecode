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
    ${GFWALG}
    ${GFWPARDAT}
    ${GPHDISXS}
    ${GPHDISXS}
    ${GPHCHMXS}
    ${GPHCOHEG} ${GPHCOHXS}
    ${GPHDFRCEG}
    ${GPHDISEG}
    ${GPHXSIG}
    ${GPHDCY}
    ${GPHQELXS}
    ${GPHQELXS}
    ${GFWGHEP}
    ${GFWEG}
    ${GPHCMN}
    ${GFWEG}
    ${GPHRESXS}
    ${GPHHADTRANSP}
    ${GPHHADNZ}
    ${GFWINT}
    ${GPHQELXS}
    ${GPHMNUCEG} ${GPHMNUCXS}
    ${GPHAMNGEG} ${GPHAMNGXS}
    ${GPHNUELEG} ${GPHNUELXS}
    ${GFWNTP}
    ${GPHNUCLST}
    ${GFWNUM}
    ${PHDISXS}
    ${GFWPARDAT}
    ${GPHPDF}
    ${GPHQELEG}
    ${GPHRESEG} ${GPHRESXS}
    ${GFWREG}
    ${GPHRESXS}
    ${GFWUTL}
    ${GTLGEO}
    ${GTLFLX}
    ${GPHMEL}
    ${GPHNDCY}
    ${GPHNNBAROSC}
    ${GFWMSG}
    ${GRWFWK}
    ${GRWIO}
    ${GRWCLC}
    ${GSL}
    ${GSLCBLAS}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_GENVECTOR}
    ${ROOT_EGPYTHIA6}
    ${ROOT_EG}
    ${GENIE_LIB_LIST}
)

