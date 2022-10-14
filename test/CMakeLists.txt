project(RelPose-Tests CXX)

# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})

# required by SDPA
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread") #  -DEIGEN_DONT_VECTORIZE") # goes to core dumped

add_library(exp_help 
                experimentsHelper.h 
                generateCorrespondences.h 
                #
                experimentsHelper.cpp
                generateCorrespondences.cpp
                )
                
                
                

# GFlags. The namespace patch is borrow from Ceres Solver (see license in
# FindGflags.cmake)

find_package(gflags REQUIRED)



# Google Logging
# message("-- Check for Google Log")
# find_package(glog 0.4.0 REQUIRED)



add_executable(generic_test ${CMAKE_CURRENT_SOURCE_DIR}/generic_test.cpp)
target_link_libraries(generic_test
                                SDPRPP 
                                sdpa 
                                lapack 
                                blas 
                                dmumps_seq
                                exp_help 
                                )
     


