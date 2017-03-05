include(CMakeParseArguments)

macro(MKNIX_TEST)
  set(options)
  set(oneValueArgs NAME OUTPUT_FILE)
  set(multiValueArgs INPUT_FILES COMMAND_ARGS)
  cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})

  file(COPY ${TEST_INPUT_FILES} DESTINATION ${TEST_DIR})

  # need to replace ; with \\; so that args gets passed as semi colon separated list to command
  string(REPLACE ";" "\\;" COMMAND_ARGS "${TEST_COMMAND_ARGS}")

  add_test(NAME test1
      COMMAND ${CMAKE_COMMAND}
      -Dtest_cmd=${CMAKE_CURRENT_BINARY_DIR}/mknixrunner
      -Dtest_args=${COMMAND_ARGS}
      -Doutput_blessed=${CMAKE_SOURCE_DIR}/tests/${TEST_OUTPUT_FILE}
      -Doutput_test=${TEST_OUTPUT_FILE}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
      WORKING_DIRECTORY ${TEST_DIR})
endmacro()
