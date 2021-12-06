# project name
project(ConvertNtuple)

find_package(nugen REQUIRED)
find_package(GENIE REQUIRED PUBLIC) # find module in nugen
find_package(log4cpp REQUIRED PUBLIC) # find module in nugen

#include_directories( $ENV{GENIE_INC}/GENIE )

find_library(XML2LIB NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib)

cet_make_exec(ConvertNtuple SOURCE ConvertNtuple.cxx BooNEtoGSimple.cxx
              LIBRARIES ${XML2LIB}
                        GENIE::GFwParDat
                        log4cpp::log4cpp
                        GENIE::GFwAlg
                        GENIE::GFwParDat
                        GENIE::GPhDISXS
                        GENIE::GPhPDF
                        GENIE::GPhXSIg
                        GENIE::GPhDISEG
                        GENIE::GPhHadnz
                        GENIE::GPhDcy
                        GENIE::GFwEG
                        GENIE::GFwEG
                        GENIE::GTlFlx
                        GENIE::GFwGHEP
                        GENIE::GFwInt
                        GENIE::GFwMsg
                        GENIE::GPhNuclSt
                        GENIE::GFwNum
                        GENIE::GFwParDat
                        GENIE::GFwReg
                        GENIE::GPhCmn
                        GENIE::GFwUtl
                        GSL
                        ${GSLCBLAS}
                        ROOT::Physics
                        ROOT::Tree
                )
