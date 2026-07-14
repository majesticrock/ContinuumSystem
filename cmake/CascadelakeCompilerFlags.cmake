# ClusterCompilerFlags.cmake

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -march=cascadelake -O3 -ffast-math)
        target_compile_definitions(${TARGET} PRIVATE NDEBUG MROCK_CL1_CASCADE)
        SET_MKL_FLAGS(${TARGET})
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only gcc|Clang is supported.")
    endif()
endfunction()