cmake_minimum_required(VERSION 3.9)
project(quant)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)
include_directories(tests)

add_executable(quant
        fd.cc
        fd.h
        payoff.cc
        payoff.h
        spline.cc
        spline.h)


find_package(Boost)
IF (Boost_FOUND)
    include_directories(${/Users/cedricwang/boost_1_67_0/boost})
endif()


#set (Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
#set (Boost_USE_MULTITHREAD ON)  # enable multithreading
#find_package (Boost COMPONENTS REQUIRED chrono filesystem)


#target_link_libraries (my_target ${Boost_LIBRARIES})