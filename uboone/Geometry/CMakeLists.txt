
art_make( BASENAME_ONLY
	  LIBRARY_NAME uboone_Geometry
	  LIB_LIBRARIES larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        cetlib cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES larcorealg_Geometry
	  		uboone_Geometry
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )							   

add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()


