MKNIX_TEST(
    NAME test1_thermalfluence
    COMMAND_ARGS input_2triangles.fem.mknix 2 0
    INPUT_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/input_2triangles.fem.mknix
        ${MKNIX_RESOURCES_DIR}/2triangles.dat
    OUTPUT_FILE 2triangfem.mec
)
