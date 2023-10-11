cet_make_exec (
  NAME ConvertNtuple
  SOURCE
  BooNEtoGSimple.cxx
  BooNEtoNuMI.cxx
  ConvertNtuple.cxx
  LIBRARIES
  PRIVATE
  GENIE::GTlFlx
  GENIE::GFwParDat
  ROOT::Tree
  ROOT::Physics
)

