cet_build_plugin(
  SparseRawDigitTest art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::RawData
  lardataobj::RawData
  lardataobj::RecoBase
)

install_headers()
install_fhicl()
install_source()
