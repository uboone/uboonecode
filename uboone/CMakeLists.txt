

set( UBOONEDAQ_DATATYPES_DIR $ENV{UBOONEDAQ_DATATYPES_DIR} )

add_subdirectory(Geometry)
add_subdirectory(RawData)
if( UBOONEDAQ_DATATYPES_DIR )
	add_subdirectory(BeamData)
	add_subdirectory(BeamDAQ)
endif()
add_subdirectory(Utilities)
add_subdirectory(Database)
add_subdirectory(Simulation)
add_subdirectory(LArG4)
add_subdirectory(PhotonPropagation)
add_subdirectory(DetSim)
add_subdirectory(EventGenerator)
add_subdirectory(CalData)
add_subdirectory(TriggerSim)
add_subdirectory(AnalysisTree)
add_subdirectory(Calibrations)
add_subdirectory(OpticalDetectorSim)
add_subdirectory(OpticalDetectorAna)
add_subdirectory(SpaceCharge)
add_subdirectory(TPCNeutrinoIDFilter)
#add_subdirectory(DataScanner)
#add_subdirectory(LiteMaker)

# Respect the order specified in lar_build, which is:
#
#     "SoftRelTools",
#     "SRT_LAR",
#     "SimpleTypesAndConstants",
#     "Geometry",
#     "SummaryData",
#     "RawData",
#     "RecoBase",
#     "OpticalDetectorData",
#     "AnalysisBase",
#     "Utilities",
#     "Simulation",
#     "MCCheater",
#     "AnalysisAlg",
#     "RecoObjects",
#     "RecoAlg",
#     "LArG4",
#     "PhotonPropagation",
#     "DetSim",
#     "EventGenerator",
#     "Filters",
#     "CalData",
#     "ClusterFinder",
#     "Genfit",
#     "HitFinder",
#     "TrackFinder",
#     "VertexFinder",
#     "EventFinder",
#     "TriggerAlgo",
#     "OpticalDetector",
#     "Calorimetry", 
#     "ShowerFinder",
#     "EventDisplay",
#     "ParticleIdentification",
#     "AnalysisExample",
#     "LArPandoraAlgorithms",
#     "LArPandoraInterface"
