cmake_minimum_required(VERSION 3.16)
project(temporal-abcore)

set(CMAKE_CXX_STANDARD 20)

add_executable(tabcore main.cpp abcore.h abcore.cpp utility.h utility.cpp bigraph.h bigraph.cpp baseline.cpp baseline.h config.h tabcore_baseline.cpp tabcore_baseline.h)