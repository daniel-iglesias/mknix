
# NOTE: W-Cp.dat and W-K.dat (temperature-dependent material properties) are
# referenced by the input but are not yet present in resources/. Add them to
# resources/ and include them in INPUT_FILES below when available.

MKNIX_TEST(
    NAME test_thermalflux_timefile
    COMMAND_ARGS input_lamella_fem.mknix 2 0
    INPUT_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/input_lamella_fem.mknix
        ${MKNIX_RESOURCES_DIR}/M15-32-278.dat
        ${MKNIX_RESOURCES_DIR}/PIC_load_scaled.txt
        ${MKNIX_RESOURCES_DIR}/M15-32-278-boundary.dat
    OUTPUT_FILE lamella_278.mec
)
