cmake_minimum_required(VERSION 3.12)

project(PIPER
    VERSION 1.0
    LANGUAGES CXX
)

add_compile_options(-O3) 
# add_compile_options(-Wall -fopenmp -O2) 


#finding required packages
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#mutationpp -----------------------------------
set (EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories("/home/inferno/Cambridge/Continuum/Project/CFDCode/Mutationpp/install/include/mutation++")
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
find_package(mutation++ REQUIRED)

#-----------------------------------------------

#specifying include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/base)
include_directories(${PROJECT_SOURCE_DIR}/include/base/utils)
include_directories(${PROJECT_SOURCE_DIR}/include/base/init)
include_directories(${PROJECT_SOURCE_DIR}/include/settings)
include_directories(${PROJECT_SOURCE_DIR}/include/classes)
include_directories(${PROJECT_SOURCE_DIR}/include/classes/source)



#making libraries ---------------------------------------------------

#base lib ----
file (GLOB_RECURSE BASE_FILES ${PROJECT_SOURCE_DIR}/lib/base/*.cpp)
add_library(baseLib OBJECT ${BASE_FILES})

#classes ----
file (GLOB_RECURSE CLASSES ${PROJECT_SOURCE_DIR}/lib/classes/*.cpp)
add_library(classesLib OBJECT ${CLASSES})

#executables ------------------------------------------------------------

#MHDEoS
add_executable(MHDEoS ${PROJECT_SOURCE_DIR}/src/MHDEoS.cpp)
set_target_properties(MHDEoS PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
target_compile_definitions(MHDEoS PRIVATE MHD_COMPILE)
target_link_libraries(MHDEoS PUBLIC baseLib classesLib mutation++)
target_sources(MHDEoS PUBLIC ${MHD_SRC})


