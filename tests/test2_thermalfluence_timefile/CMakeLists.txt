MKNIX_TEST(
    NAME test2_thermalfluence_timefile
    COMMAND_ARGS input_2triangles_timefile.fem.mknix 2 0
    INPUT_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/input_2triangles_timefile.fem.mknix
        ${MKNIX_RESOURCES_DIR}/2triangles.dat
        ${MKNIX_RESOURCES_DIR}/timescale.dat
    OUTPUT_FILE 2triangfem_timefile.mec
)
