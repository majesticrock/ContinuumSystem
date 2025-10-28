# MKLFlags.cmake

function(SET_MKL_FLAGS TARGET)
    if (USE_MKL)
        message(STATUS "Configuring target ${TARGET} to use Intel MKL")

        target_compile_options(${TARGET} PRIVATE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(${TARGET} PRIVATE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(${TARGET} PRIVATE $<LINK_ONLY:MKL::MKL>)

	    # Let Eigen use MKL
	    target_compile_definitions(${TARGET} PRIVATE EIGEN_USE_MKL_ALL MROCK_IEOM_DO_NOT_PARALLELIZE)
    else()
        message(STATUS "MKL not found or not enabled; ${TARGET} will not use Intel MKL")
    endif()
endfunction()
