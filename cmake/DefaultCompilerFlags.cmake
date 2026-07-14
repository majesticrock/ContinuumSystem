# DefaultCompilerFlags.cmake

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -march=native -O0 -g)# -fsanitize=undefined,address
            #target_link_options(${TARGET} PRIVATE -fsanitize=undefined,address)
        else()
            target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -march=native -O3 -ffast-math)
            if(CMAKE_BUILD_TYPE STREQUAL "NDEBUG")
                target_compile_definitions(${TARGET} PRIVATE NDEBUG)
            endif()
            SET_MKL_FLAGS(${TARGET})
        endif()
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC/Clang is supported.")
    endif()
endfunction()