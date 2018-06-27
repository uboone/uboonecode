
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-nested-anon-types )
endif()

# Main library.


# Subdirectories.

add_subdirectory(Overlays)
add_subdirectory(BernZMQBinarySource)

