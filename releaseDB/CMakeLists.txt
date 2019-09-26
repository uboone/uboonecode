
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

create_product_variables( uboonecode ${version} )

# define qualifiers
set( BASEQUAL e17 )
set( SQUAL s78 )

create_product_variables( uboone_photon_propagation v01_00_00  )
create_product_variables( larbatch                  v01_51_04  )
create_product_variables( cetbuildtools             v7_11_00  )
create_product_variables( swtrigger                 v02_02_05i )
create_product_variables( larlite                   v06_70_00c )
create_product_variables( uboone_oslibs             v1_0_0     )
create_product_variables( uboonedaq_datatypes       v6_22_05a  )
create_product_variables( uboonedata                v01_18_22_08  )
create_product_variables( ubutil                    v08_00_00_19  )
create_product_variables( ublite                    v08_00_00_22  )
create_product_variables( ubana                     v08_00_00_22  )
create_product_variables( ubreco                    v08_00_00_22  )
create_product_variables( ubsim                     v08_00_00_22  )
create_product_variables( ubcv                      v08_00_00_22  )
create_product_variables( ubevt                     v08_00_00_22  )
create_product_variables( ubraw                     v08_00_00_22  )
create_product_variables( ubcore                    v08_00_00_22  )
create_product_variables( ubcrt                     v08_00_00_22  )
create_product_variables( ubobj                     v08_00_00_22  )
create_product_variables( larsoft                   v08_05_00_10  )

# temporary
create_product_variables( cetpkgsupport v1_14_01 )

# these files should be copied to scisoft for each release
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone-cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-cfg-${UBOONECODE_DOT_VERSION} @ONLY )
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone-${version}.html @ONLY )

INSTALL( FILES base_dependency_database
               README
         DESTINATION "${product}/${version}/releaseDB" )
