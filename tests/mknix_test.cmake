include(CMakeParseArguments)

# MKNIX_TEST - register a simulation test case.
#
# Required variables set by the parent tests/CMakeLists.txt:
#   MKNIX_RUNNER_PATH      - absolute path to the mknixrunner executable
#   MKNIX_RUN_TEST_SCRIPT  - absolute path to run_test.cmake
#   MKNIX_RESOURCES_DIR    - absolute path to the shared resources/ directory
#
# Parameters:
#   NAME          - unique test name (also used as the CTest test identifier)
#   COMMAND_ARGS  - arguments forwarded to mknixrunner
#   INPUT_FILES   - files to copy into the test working directory (use absolute paths)
#   OUTPUT_FILE   - name of the output file produced by the simulation;
#                   the blessed reference must reside in CMAKE_CURRENT_SOURCE_DIR
#                   under this name

macro(MKNIX_TEST)
  set(options)
  set(oneValueArgs NAME OUTPUT_FILE)
  set(multiValueArgs INPUT_FILES COMMAND_ARGS)
  cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  file(COPY ${TEST_INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # need to replace ; with \\; so that args gets passed as semi colon separated list to command
  string(REPLACE ";" "\\;" COMMAND_ARGS "${TEST_COMMAND_ARGS}")

  add_test(NAME ${TEST_NAME}
      COMMAND ${CMAKE_COMMAND}
      -Dtest_cmd=${MKNIX_RUNNER_PATH}
      -Dtest_args=${COMMAND_ARGS}
      -Doutput_blessed=${CMAKE_CURRENT_SOURCE_DIR}/${TEST_OUTPUT_FILE}
      -Doutput_test=${TEST_OUTPUT_FILE}
      -P ${MKNIX_RUN_TEST_SCRIPT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro()
