# ContinuumSystem/CMakeLists.txt
cmake_minimum_required(VERSION 3.30)

project(ContinuumSystem VERSION 1.0)

# Set local user-specific settings
# By default, this adds ~/usr/local/ and ../../.mrock/ to the cmake and include paths
# If additional configurations/includes are needed, provide them in cmake/user_settings.cmake
# Extra include dirs are stored in the cmake variable EXTRA_INCLUDE_DIRS
include(cmake/include_settings.cmake)

# Build type settings; defaults to Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#
# Library specifics
#

add_executable(ContinuumSystem 
    sources/ContinuumSystem.cpp
)
set(ADD_SOURCES
    sources/Continuum/ModeHelper.cpp
    sources/Continuum/SCModel.cpp
    sources/Continuum/MomentumRanges.cpp
    sources/Continuum/SplineContainer.cpp
    sources/Continuum/ModelInitializer.cpp
    sources/Continuum/Incrementer.cpp
    sources/Continuum/PhononInteraction.cpp
)
target_sources(ContinuumSystem PRIVATE ${ADD_SOURCES})

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
        iEoM 
        utility
)

message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CXX compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OpenMP_CXX_FOUND: ${OpenMP_CXX_FOUND}")
message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
message(STATUS "OpenMP_CXX_INCLUDE_DIRS: ${OpenMP_CXX_INCLUDE_DIRS}")
message(STATUS "OpenMP_CXX_LIBRARIES: ${OpenMP_CXX_LIBRARIES}")

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/continuum/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUT_HEADER MROCK_GENERATED_INFO_HEADER
    OUT_INCLUDE_DIR MROCK_GENERATED_INCLUDE_DIR
)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/continuum_options.cmake)
# Sets include dir for the generated metadata header file
target_include_directories(continuum_options INTERFACE ${MROCK_GENERATED_INCLUDE_DIR})


target_link_libraries(ContinuumSystem 
    PRIVATE 
        OpenMP::OpenMP_CXX 
        ${mrock_LIBRARIES}
        continuum_options
)

# Add debug option
option(ENABLE_DEBUG "Enable debugging" OFF)

if(ENABLE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
endif()

option(CONTINUUM_FULL_DIAG "Fully diagonalize the transformed dynamical matrix." OFF)
if (CONTINUUM_FULL_DIAG)
    target_compile_definitions(ContinuumSystem PRIVATE CONTINUUM_FULL_DIAG)
endif()

option(BUILD_TESTING 
    "Build and configure the testing suit" 
    OFF
)

if(BUILD_TESTING)
    include(tests/test_config.cmake)
endif()