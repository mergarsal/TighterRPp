# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
project(RelPose LANGUAGES C CXX VERSION 1.0.0)
set(CMAKE_BUILD_TYPE Release)



set(LIBRARY_TARGET_NAME "RelPose")
set(LIBRARY_TARGET_NAME_EXPORT "${LIBRARY_TARGET_NAME}Export")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14


# build the tests
set(BUILD_${LIBRARY_TARGET_NAME}_TESTS ON)



# Build type

# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

set(CODE_PROFILING OFF CACHE BOOL "Turn on code profiling?")
if(${CODE_PROFILING})
message(STATUS "Turning on code profiling for Essential Matrix Estimation")
endif()

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" CACHE STRING "The CMake module path used for this project")

# FIND EIGEN3
set( ENV{EIGEN3_ROOT_DIR} ${CMAKE_SOURCE_DIR}/eigen)
find_package(Eigen3 3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION_STRING})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# FIND ADDITIONAL LIBRARIES
# These next operations make use of the .cmake files shipped with Eigen3
find_package(BLAS REQUIRED)


add_definitions(-march=native)  


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# for SDP solver
find_package(SDPA REQUIRED)
include_directories(${SDPA_INCLUDE_DIRS})



set(${LIBRARY_TARGET_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(${LIBRARY_TARGET_NAME}_SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(${LIBRARY_TARGET_NAME}_TESTS_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/test)


# Expose the include directories for this project
set(${LIBRARY_TARGET_NAME}_ADD_INCLUDES ${ESSENTIAL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} )
set(${LIBRARY_TARGET_NAME}_CERT_INCLUDES ${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS} ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES})



# Get the set of Essential header and source files
set(${LIBRARY_TARGET_NAME}_HDRS
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalTypes.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalUtils.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalLeft.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalRight.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalBoth.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/PrimalAdjugate.h
)


set(${LIBRARY_TARGET_NAME}_SRCS
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/PrimalUtils.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/PrimalLeft.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/PrimalRight.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/PrimalBoth.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/PrimalAdjugate.cpp
)




# Build the Essential library
add_library(${LIBRARY_TARGET_NAME} ${${LIBRARY_TARGET_NAME}_HDRS} ${${LIBRARY_TARGET_NAME}_SRCS}  )



target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC
        # only when building from the source tree
        $<BUILD_INTERFACE:${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}>
        # only when using the lib from the install path
        $<INSTALL_INTERFACE:include>
        ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES}
)


target_link_libraries(${LIBRARY_TARGET_NAME}
        PUBLIC ${BLAS_LIBRARIES} ${M} ${LAPACK} )




if(${CODE_PROFILING})
  set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-pg -g" LINK_FLAGS "-pg -g")
endif()


set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS ${LIBRARY_TARGET_NAME}
  EXPORT  ${LIBRARY_TARGET_NAME_EXPORT}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  INCLUDES DESTINATION "include"
  PUBLIC_HEADER DESTINATION "include/${LIBRARY_TARGET_NAME}"
 )



# building the tests
    if(BUILD_${LIBRARY_TARGET_NAME}_TESTS)
    message(STATUS "Adding tests to build")
    add_subdirectory(test)
endif()
