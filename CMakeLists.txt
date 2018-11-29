cmake_minimum_required(VERSION 3.12)
project(typeDsingle)

set(CMAKE_CXX_STANDARD 14)

add_executable(typeDsingle main.cpp Coords.h ray.cpp ray.h)