# cmake version
cmake_minimum_required(VERSION 3.0.0)

# project name
project(lennard_jones VERSION 0.1.0)

# include(CTest)
enable_testing()
add_subdirectory(tests)

# flags

# files

# include
include_directories(./include ./src)

# target
add_executable(lennard_jones ./src/main.cpp ./src/Vector2D.cpp ./src/Particle.cpp ./src/Cell.cpp ./src/Grid.cpp)
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