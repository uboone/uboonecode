# include the subfolders which have fhicl files by epoch
add_subdirectory(Feb2016)
add_subdirectory(MCC7)
add_subdirectory(MCC8)
add_subdirectory(MCC9)
add_subdirectory(MCC10)

# Install fcl files in /job subdirectory.
install_fhicl()
