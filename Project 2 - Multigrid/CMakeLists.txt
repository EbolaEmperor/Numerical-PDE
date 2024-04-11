cmake_minimum_required(VERSION 3.5)
project(Multigrid)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE RELEASE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-result -Wno-non-template-friend -Wno-unused-variable -DDIM=${DIM}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DBGLEVEL "The verbosity. " -1)

if(DBGLEVEL GREATER -1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDBGLEVEL=${DBGLEVEL}")
endif()

option(USE_OPENMP "Whether to use the OpenMP parallelism." OFF)

if(USE_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP -fopenmp")
endif()

option(PROFILE "Whether to enable profiling. " OFF)

if(PROFILE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
endif()

message("System : ${CMAKE_SYSTEM_NAME}")
message("Build type : ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

set(ROOT_DIR "${CMAKE_SOURCE_DIR}")
set(MAIN_LIB "Multigrid")
set(COMMON_LIB "${MAIN_LIB}_common")

find_package(jsoncpp)

add_subdirectory(src)
add_subdirectory(examples)
