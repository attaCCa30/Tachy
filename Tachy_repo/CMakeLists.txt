cmake_minimum_required(VERSION 3.16)
project(Tachy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_MPI "Enable MPI parallelism" ON)
option(ENABLE_OPENMP "Enable OpenMP threading" ON)

add_executable(tachy
  src/main.cpp
)

target_include_directories(tachy PRIVATE include)

if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(tachy PRIVATE MPI::MPI_CXX)
  target_compile_definitions(tachy PRIVATE ENABLE_MPI)
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(tachy PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
