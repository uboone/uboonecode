art_make( LIB_LIBRARIES
          #  larcore_Geometry
          #  ${MF_MESSAGELOGGER}
          #  ${MF_UTILITIES}
            cetlib cetlib_except
            ${FHICLCPP}
            ${ROOT_CORE}
            ${ROOT_RIO}
            ${ROOT_HIST}
            ${GSL}
        )

install_headers()
install_fhicl()
install_source()
