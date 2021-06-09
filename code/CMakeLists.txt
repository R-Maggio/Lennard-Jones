# cmake version
cmake_minimum_required(VERSION 3.9)

# project name
project(lennard_jones VERSION 0.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_CXX_FLAGS "-O3 -Werror -Wall -Wextra -Wpedantic")
# set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS "-O3 -Wall")
set(OpenMP_CXX_FLAGS "-fopenmp -lpthread")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# include(CTest)
enable_testing()
add_subdirectory(tests)

# flags

# files

# include
include_directories(./include ./src)

# target
add_executable(lennard_jones ./src/main.cpp ./src/Vector2D.cpp ./src/Particle.cpp ./src/Cell.cpp ./src/Grid.cpp ./src/LJSimulation.cpp)
# add_executable( binary ./src/main.cpp ./src/game_engine.cpp ./src/game_interface.cpp )

# # external libs
# add_subdirectory( ./libs/Logger )
# target_link_libraries( binary logger )

# # add_subdirectory( ./libs/Randomize )
# target_link_libraries( binary randomize )

# # add_library(randomize STATIC IMPORTED)
# set_property(TARGET randomize PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/Randomize/librandomize.a)
# target_link_libraries( binary randomize )


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)