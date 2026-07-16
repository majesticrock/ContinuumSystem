# Custom architecture target
set(TARGET_ARCH
    "native"
    CACHE STRING
    "Architecture passed to compiler as -march=<arch> (empty disables)"
)

# Joined compile options
add_library(continuum_options INTERFACE)
target_compile_features(continuum_options
    INTERFACE
        cxx_std_20
)

target_compile_options(continuum_options INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
)
if(TARGET_ARCH)
    target_compile_options(continuum_options INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-march=${TARGET_ARCH}>
    )
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(continuum_options INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-ffast-math> 
    )
else()
    target_compile_definitions(continuum_options INTERFACE DEBUG)
    target_compile_options(continuum_options INTERFACE -g)
endif()

#
# Use MKL if available
#
set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
set(MKL_INTERFACE lp64)

find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT})

if (MKL_FOUND)
    message(STATUS "Configuring Continuum to use Intel MKL")

    target_compile_options(continuum_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(continuum_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(continuum_options INTERFACE $<LINK_ONLY:MKL::MKL>)
    # Let Eigen use MKL
    target_compile_definitions(continuum_options INTERFACE EIGEN_USE_MKL_ALL MROCK_IEOM_DO_NOT_PARALLELIZE)
else()
    message(STATUS "MKL not found or not enabled; Continuum will not use Intel MKL")
endif()