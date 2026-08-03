
# NOTE: W-Cp.dat and W-K.dat (temperature-dependent material properties) are
# referenced by the input but are not yet present in resources/. Add them to
# resources/ and include them in INPUT_FILES below when available.

MKNIX_TEST(
    NAME test6_thermalflux_2Dmap
    COMMAND_ARGS input_lamella_fem_2D.mknix 2 0
    INPUT_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/input_lamella_fem_2D.mknix
        ${MKNIX_RESOURCES_DIR}/M15-32-278.dat
        ${MKNIX_RESOURCES_DIR}/PIC_load_scaled_2D.txt
        ${MKNIX_RESOURCES_DIR}/M15-32-278-boundary.dat
    OUTPUT_FILE lamella_278.mec
)
