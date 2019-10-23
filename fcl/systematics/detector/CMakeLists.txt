# include the subfolders which have fhicl files by epoch
add_subdirectory(MCC8)
add_subdirectory(MCC9)

# Install fcl files in /job subdirectory.
install_fhicl()
