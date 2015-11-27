# compile in debug mode
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# FIND_PACKAGE(VTK)
# IF(VTK_FOUND)
#   INCLUDE(${VTK_USE_FILE})
#   INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
# ELSE(VTK_FOUND)
# #  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
#   FIND_PATH(VTK_DIR VTKConfig.cmake /usr/local/lib/vtk-5.4)
#   INCLUDE(${VTK_USE_FILE})
# #  FIND_PATH(VTK_INCLUDE_DIR vtkActor.h /usr/local/include/vtk-5.4)
# ENDIF(VTK_FOUND)
# 
# FIND_PATH(LMX_INCLUDE_DIR HINTS ${PROJECT_SOURCE_DIR}/LMX-1.0/src NAMES LMX/lmx.h)

# Make sure the compiler can find include files.
include_directories (
${MkniX_SOURCE_DIR}/src 
${MkniX_SOURCE_DIR}/libmknix 
${LMX_INCLUDE_DIR}/ 
# ${VTK_INCLUDE_DIR}
${GMM_INCLUDE_DIR}
)

# Make sure the linker can find the Hello library once it is built.
link_directories (${MkniX_BINARY_DIR}/libmknix)

add_executable (MkniX mknix.cpp)


# link_directories ( ${VTK_DIR})
# Link the executable to the VTK libraries.
# target_link_libraries ( mknixlib
#   vtkRendering
#   vtkGraphics
#   vtkIO
#   vtkCommon
# )

link_directories ( ${PROJECT_BINARY_DIR})
# Link the executable to the VTK libraries.
target_link_libraries ( MkniX
   mknixlib
)
