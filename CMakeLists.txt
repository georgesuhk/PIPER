cmake_minimum_required(VERSION 3.12)

project(PIPER
    VERSION 1.0
    LANGUAGES CXX
)

add_compile_options(-fopenmp -O3) 
# add_compile_options(-Wall -fopenmp -O2) 


# finding required packages
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# mutationpp -----------------------------------
set (EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories("/home/inferno/Cambridge/Continuum/Project/CFDCode/Mutationpp/install/include/mutation++")
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
find_package(mutation++ REQUIRED)

#-----------------------------------------------

# specifying include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/base)
include_directories(${PROJECT_SOURCE_DIR}/include/base/utils)
include_directories(${PROJECT_SOURCE_DIR}/include/base/init)
include_directories(${PROJECT_SOURCE_DIR}/include/settings)
include_directories(${PROJECT_SOURCE_DIR}/include/classes)
include_directories(${PROJECT_SOURCE_DIR}/include/classes/source)
include_directories(${PROJECT_SOURCE_DIR}/include/other)




# making libraries ---------------------------------------------------

# base lib ----
file (GLOB_RECURSE BASE_FILES ${PROJECT_SOURCE_DIR}/lib/base/*.cpp)
add_library(baseLib OBJECT ${BASE_FILES})

# classes lib ----
file (GLOB_RECURSE CLASSES ${PROJECT_SOURCE_DIR}/lib/classes/*.cpp)
add_library(classesLib OBJECT ${CLASSES})

# others lib ----
file (GLOB_RECURSE OTHER_FILES ${PROJECT_SOURCE_DIR}/lib/other/*.cpp)
add_library(otherLib OBJECT ${OTHER_FILES})

# EXECUTABLES ------------------------------------------------------------

# MHDEoS
add_executable(MHDEoS ${PROJECT_SOURCE_DIR}/src/MHDEoS.cpp)
set_target_properties(MHDEoS PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(MHDEoS PUBLIC baseLib classesLib mutation++)

# orzag
add_executable(orzag ${PROJECT_SOURCE_DIR}/src/orzag.cpp)
set_target_properties(orzag PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(orzag PUBLIC baseLib classesLib mutation++)

# detachment
add_executable(detachment ${PROJECT_SOURCE_DIR}/src/detachment.cpp)
set_target_properties(detachment PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(detachment PUBLIC baseLib classesLib mutation++)

# detachmentScan
add_executable(detachmentScan ${PROJECT_SOURCE_DIR}/src/detachmentScan.cpp)
set_target_properties(detachmentScan PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(detachmentScan PUBLIC baseLib classesLib mutation++)


# MHDMulti
add_executable(MHDMulti ${PROJECT_SOURCE_DIR}/src/MHDMulti.cpp)
set_target_properties(MHDMulti PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(MHDMulti PUBLIC baseLib classesLib mutation++)


# EXECUTABLES FOR CONSTRUCTION --------------------------------------------

# EoSConstruct
add_executable(EoSConstruct ${PROJECT_SOURCE_DIR}/src/construct/EoSConstruct.cpp)
set_target_properties(EoSConstruct PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(EoSConstruct PUBLIC baseLib classesLib otherLib mutation++)

# resisCompare
add_executable(resisCompare ${PROJECT_SOURCE_DIR}/src/construct/study/resisCompare.cpp)
set_target_properties(resisCompare PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(resisCompare PUBLIC baseLib classesLib mutation++)

# collRateStudy (1D)
add_executable(collRateStudy ${PROJECT_SOURCE_DIR}/src/construct/study/collRateStudy.cpp)
set_target_properties(collRateStudy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(collRateStudy PUBLIC baseLib classesLib mutation++)

# collRate2D (2D)
add_executable(collRate2D ${PROJECT_SOURCE_DIR}/src/construct/study/collRate2D.cpp)
set_target_properties(collRate2D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(collRate2D PUBLIC baseLib classesLib otherLib mutation++)

# resisFinder 
add_executable(resisFinder ${PROJECT_SOURCE_DIR}/src/construct/study/resisFinder.cpp)
set_target_properties(resisFinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_link_libraries(resisFinder PUBLIC baseLib classesLib otherLib mutation++)
