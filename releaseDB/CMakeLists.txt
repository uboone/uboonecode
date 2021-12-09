
# variables to be used in cmake templates
macro( create_product_variables productname productversion )
   string(TOUPPER ${productname} PRODUCTNAME_UC )
   string(TOLOWER ${productname} PRODUCTNAME_LC )

   # require ${PRODUCTNAME_UC}_VERSION
   set( ${PRODUCTNAME_UC}_VERSION "${productversion}" )
   if ( NOT ${PRODUCTNAME_UC}_VERSION )
     message(FATAL_ERROR "product version is not defined")
   endif ()
   STRING( REGEX REPLACE "_" "." VDOT "${productversion}"  )
   STRING( REGEX REPLACE "^[v]" "" ${PRODUCTNAME_UC}_DOT_VERSION "${VDOT}" )
   set( ${PRODUCTNAME_UC}_STRING ${PRODUCTNAME_LC}.${${PRODUCTNAME_UC}_DOT_VERSION} )
   message(STATUS "${PRODUCTNAME_UC} version is ${${PRODUCTNAME_UC}_VERSION}")
endmacro( create_product_variables)

create_product_variables( uboonecode v09_38_00 )

# define qualifiers
set( BASEQUAL e20 )
set( SQUAL s112 )

create_product_variables( uboone_photon_propagation v01_00_00  )
create_product_variables( swtrigger v02_03_00d )
create_product_variables( larlite v06_71_00d )
create_product_variables( larcv v06_71_00d )
create_product_variables( uboone_oslibs             v1_0_0     )
create_product_variables( uboonedaq_datatypes v6_22_05g )
create_product_variables( uboonedata v09_32_00  )
create_product_variables( ubutil v09_32_00 )
create_product_variables( ublite v09_37_00 )
create_product_variables( ubana v09_37_00 )
create_product_variables( ubreco v09_37_00 )
create_product_variables( ubsim v09_37_00 )
create_product_variables( ubcv v09_37_00 )
create_product_variables( ubevt v09_37_00 )
create_product_variables( ubraw v09_37_00 )
create_product_variables( ubcore v09_37_00 )
create_product_variables( ubcrt v09_37_00 )
create_product_variables( ubobj v09_37_00 )
create_product_variables( larsoft v09_38_00 )

# temporary
create_product_variables( cetpkgsupport v1_14_01 )

# these files should be copied to scisoft for each release
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone-cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-cfg-${UBOONECODE_DOT_VERSION} @ONLY )
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-${version}.html @ONLY )

INSTALL( FILES base_dependency_database README
         DESTINATION "releaseDB" )
