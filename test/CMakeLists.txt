# This CMake script configures unit tests and regression tests for exciting

# Option to activate unit tests. Default is ON.
option(UNIT_TESTS "Activate the unit tests for the build system" ON)

# Option to activate regression tests. Default is ON.
option(REGRESSION_TESTS "Activate the regression tests for the build system" OFF)

# Set the path to the installed executable for Exciting.
# This is derived from the CMake installation prefix.
set(INSTALLED_EXECUTABLE "${CMAKE_INSTALL_PREFIX}/bin/exciting")

# Python required for the regression tests
if (REGRESSION_TESTS)
  include(${CMAKE_SOURCE_DIR}/cmake/python3.cmake)
  find_python_module(numpy)
  find_python_module(pytest)
endif()


# Add checks for unit tests if the UNIT_TESTS option is enabled.
if (UNIT_TESTS)
    # Define a test named 'exciting_unit_tests'.
    add_test(
        NAME exciting_unit_tests
        # Command to execute the Exciting binary with unit test flags.
        COMMAND "${INSTALLED_EXECUTABLE}" -run-unit-tests -kill-on-failure
    )
    # Set a timeout for the unit tests. They should not exceed 15 minutes.
    # This timeout helps prevent hanging compilations or tests.
    set_tests_properties(exciting_unit_tests PROPERTIES TIMEOUT 900)
endif()

# Add regression tests if the REGRESSION_TESTS option is enabled.
if (REGRESSION_TESTS)

    # Check if Python3 is found, as it is required for running regression tests.
    if (NOT Python3_FOUND)
        message(FATAL_ERROR "Regression tests require Python3")
    endif()

    # List of regression test names to be run.
    set(TESTS "groundstate" "GW" "BSE" "fastBSE" "hybrid" "properties" "RT-TDDFT" "sirius" "TDDFT" "XANES")

    # Iterate over each test name and add it to the test suite.
    foreach(TEST_NAME ${TESTS})
        add_test(
            NAME exciting_${TEST_NAME}_tests
            # Set the working directory for the tests to the test folder in the source directory.
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test/"
            # Command to run the regression tests using Python3.
            # The 'runtest.py' script is executed with necessary parameters.
            COMMAND ${Python3_EXECUTABLE} runtest.py -e ${exciting_name} -r ${CMAKE_INSTALL_PREFIX}/bin/ -t  ${TEST_NAME}
        )
        # Set a timeout for each regression test. They should not exceed 1 hour.
        set_tests_properties(exciting_${TEST_NAME}_tests PROPERTIES TIMEOUT 3600)
    endforeach()

endif()
