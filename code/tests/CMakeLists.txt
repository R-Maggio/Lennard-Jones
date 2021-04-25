
set(INCLUDE_DIR, "./include")

add_executable(test_grid Grid_TEST.cpp ../src/Vector2D.cpp ../src/Particle.cpp ../src/Cell.cpp ../src/Grid.cpp)
target_include_directories(test_grid PRIVATE ../include)
set_target_properties(test_grid PROPERTIES LINKER_LANGUAGE CXX)

add_test(NAME test_grid COMMAND test_grid)

add_executable(test_particle Particle_TEST.cpp ../src/Vector2D.cpp ../src/Particle.cpp ../src/Cell.cpp ../src/Grid.cpp)
target_include_directories(test_particle PRIVATE ../include)
set_target_properties(test_particle PROPERTIES LINKER_LANGUAGE CXX)

add_test(NAME test_particle COMMAND test_particle)

add_executable(test_LJSimulation LJSimulation_TEST.cpp ../src/Vector2D.cpp ../src/Particle.cpp ../src/Cell.cpp ../src/Grid.cpp ../src/LJSimulation.cpp)
target_include_directories(test_LJSimulation PRIVATE ../include)
set_target_properties(test_LJSimulation PROPERTIES LINKER_LANGUAGE CXX)

add_test(NAME test_LJSimulation COMMAND test_LJSimulation)