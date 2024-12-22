cmake_minimum_required(VERSION 3.10)
project(WaveEquationSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(wave_solver main.cpp)

install(TARGETS wave_solver DESTINATION bin)
install(FILES wave_solver.hpp DESTINATION include)
