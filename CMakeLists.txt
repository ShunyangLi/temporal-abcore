cmake_minimum_required(VERSION 3.16)
project(temporal-abcore)

set(CMAKE_CXX_STANDARD 20)

add_executable(tabcore main.cpp abcore/abcore.h abcore/abcore.cpp utility/utility.h utility/utility.cpp bigraph/bigraph.h bigraph/bigraph.cpp baseline/baseline.cpp baseline/baseline.h config/config.h baseline/tabcore_baseline.cpp baseline/tabcore_baseline.h online/online.cpp online/online.h)