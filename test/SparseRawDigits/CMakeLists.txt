simple_plugin(SparseRawDigitTest "module"
	      ubobj_RawData
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art_root_io::TFileService_service
              art::Persistency_Common
              art::Persistency_Provenance
              art::Utilities
              canvas::canvas
              messagefacility::MF_MessageLogger
              fhiclcpp::fhiclcpp
              cetlib::cetlib
              cetlib_except::cetlib_except
              PostgreSQL::PostgreSQL
              BASENAME_ONLY
             )

install_headers()
install_fhicl()
install_source()
