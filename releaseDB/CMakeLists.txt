
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

create_product_variables( uboonecode v10_07_00 )

# define qualifiers
set( BASEQUAL e26 )
set( SQUAL s131 )

create_product_variables( uboone_photon_propagation v01_01_00  )
create_product_variables( larbatch v01_61_01 )
create_product_variables( swtrigger v02_03_00i )
#create_product_variables( larlite v06_71_00i )
#create_product_variables( larcv v06_71_00i )
create_product_variables( ubdl v2_me_06_03b )
create_product_variables( glpk                      v4_65     )
create_product_variables( geant4reweight v01_20_11 )
create_product_variables( wcp v00_18_01 )
create_product_variables( scn v01_00_00 )
create_product_variables( uboone_oslibs             v1_0_0     )
create_product_variables( uboonedaq_datatypes v6_22_07 )
create_product_variables( uboonedata v10_04_04 )
create_product_variables( ubutil v10_07_00 )
create_product_variables( ublite v10_07_00 )
create_product_variables( ubana v10_07_00 )
create_product_variables( ubreco v10_07_00 )
create_product_variables( ubsim v10_07_00 )
create_product_variables( ubcv v10_07_00 )
create_product_variables( ubevt v10_06_00 )
create_product_variables( ubraw v10_06_00 )
create_product_variables( ubcore v10_06_00 )
create_product_variables( ubcrt v10_07_00 )
create_product_variables( ubobj v10_06_00 )
create_product_variables( larsoft v10_07_00 )

# temporary
create_product_variables( cetpkgsupport v1_14_01 )

# these files should be copied to scisoft for each release
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone-cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-cfg-${UBOONECODE_DOT_VERSION} @ONLY )
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-${${CETMODULES_CURRENT_PROJECT_NAME}_UPS_PRODUCT_VERSION}.html @ONLY )

INSTALL( FILES base_dependency_database README
         DESTINATION "releaseDB" )
