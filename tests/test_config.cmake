include(CTest)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

file(GLOB TEST_CONFIGS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_LIST_DIR}/*.config"
)

foreach(CONFIG_FILE ${TEST_CONFIGS})
    get_filename_component(TEST_NAME ${CONFIG_FILE} NAME_WE)
    add_test(
        NAME "continuum_${TEST_NAME}"
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_LIST_DIR}/run_test.py"
            --exe "$<TARGET_FILE:ContinuumSystem>"
            --config "${CONFIG_FILE}"
            --plot-script "${CMAKE_CURRENT_LIST_DIR}/plot.py"
    )

    set_tests_properties("continuum_${TEST_NAME}"
        PROPERTIES
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endforeach()
