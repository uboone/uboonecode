art_make( LIBRARY_NAME crt_CRTRawBinarySource
                 LIB_LIBRARIES
                                ${FHICLCPP}
                                art_Persistency_Provenance                              
                 SOURCE_LIBRARIES
#                               artdaq_ArtModules
                                art_Framework_Core
                                art_Framework_IO_Sources
                                art_Framework_Principal
                                art_Persistency_Provenance
                                canvas
                                ${CANVAS_LIB_LIST}
                                MF_MessageLogger
                                ${ROOT_CINTEX}
                                ${ROOT_RIO}
                                ${ROOT_CORE}
                                ${ROOT_NET}
                                ${ROOT_REFLEX}
                                cetlib
                                ${CETLIB_EXCEPT}
                                ${BOOST_LIB}
                                ${FHICLCPP}
                 BASENAME_ONLY
               )



install_headers()
install_source()
