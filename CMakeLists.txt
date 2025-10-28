# ContinuumSystem/CMakeLists.txt

cmake_minimum_required(VERSION 3.16)

project(ContinuumSystem VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
set(MKL_INTERFACE lp64)

find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT})
if (MKL_FOUND)
    set(USE_MKL ON CACHE BOOL "Use Intel MKL" FORCE)
else()
    set(USE_MKL OFF CACHE BOOL "Use Intel MKL" FORCE)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/MKLFlags.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/SelectClusterBuildOptions.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../PhdUtility/cmake/mrock-information.cmake)

include_directories(
    $ENV{HOME}/usr/local/include
)

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

SET_COMPILER_FLAGS(ContinuumSystem)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(mrock REQUIRED)

target_link_libraries(ContinuumSystem PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX ${mrock_LIBRARIES})
target_include_directories(ContinuumSystem PRIVATE ${mrock_INCLUDE_DIRS})

# Add debug option
option(ENABLE_DEBUG "Enable debugging" OFF)

if(ENABLE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    target_compile_definitions(ContinuumSystem PRIVATE DEBUG)
    target_compile_options(ContinuumSystem PRIVATE -g)
endif()