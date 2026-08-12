#-------------------------------------------------
# some argument checking:
# test_cmd is the command to run with all its arguments
if (NOT test_cmd)
  message(FATAL_ERROR "Variable test_cmd not defined")
endif (NOT test_cmd)

# output_blessed contains the name of the "blessed" output file
if (NOT output_blessed)
  message(FATAL_ERROR "Variable output_blessed not defined")
endif (NOT output_blessed)

# output_test contains the name of the output file the test_cmd will  
#produce
if (NOT output_test)
  message(FATAL_ERROR "Variable output_test not defined")
endif (NOT output_test)

# convert the space-separated string to a list
#separate_arguments(test_args)

execute_process(
    COMMAND ${test_cmd} ${test_args}
    RESULT_VARIABLE test_had_error
    OUTPUT_VARIABLE test_out
    ERROR_VARIABLE test_err
)

if (test_had_error)
  message(STATUS "test_had_error: ${test_had_error}")
  message(STATUS "output: ${test_out}")
  message(STATUS "error: ${test_err}")
  message(FATAL_ERROR "Test failed")
endif (test_had_error)

# Read files and filter out GIT_* lines
file(READ ${output_blessed} blessed_content)
file(READ ${output_test} test_content)

string(REGEX REPLACE "(^|\n)GIT_[^\n]*" "" blessed_filtered "${blessed_content}")
string(REGEX REPLACE "(^|\n)GIT_[^\n]*" "" test_filtered "${test_content}")

if (NOT blessed_filtered STREQUAL test_filtered)
  message(SEND_ERROR "${output_test} does not match ${output_blessed}!")
endif ()
#------------------------------------------------- 
