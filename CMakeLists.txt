cmake_minimum_required(VERSION 3.6)
project(coinfinder)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Boost_USE_STATIC_LIBS ON)

set(SOURCE_FILES main.cpp id_lookup.cpp id_lookup.h dataset.cpp dataset.h coincidence.cpp coincidence.h connectivity.cpp connectivity.h binomial_test.cpp binomial_test.h significance.cpp significance.h constants.h parameters.cpp parameters.h elements.cpp elements.h test_cases.cpp test_cases.h phylomax.py common_ancestor.py bugfix.h bugfix.cpp lineage.cpp lineage.h gexf.cpp gexf.h network.cpp network.h)
add_executable(coinfinder ${SOURCE_FILES})

target_compile_features(coinfinder PRIVATE cxx_std_14)

#SET(PYTHON_INCLUDE_PATH /Library/Frameworks/Python.framework/Versions/3.6/include/python3.6m/)
#SET(PYTHON_LIBRARIES /Library/Frameworks/Python.framework/Versions/3.6/lib/libpython3.6m.dylib)
#SET(PYTHON_EXECUTABLE /Library/Frameworks/Python.framework/Versions/3.6/bin/python3)
#SET(PYTHON_INCLUDE_DIR /Library/Frameworks/Python.framework/Versions/3.6/include/python3.6m/)
find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonLibs 3.6 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(coinfinder ${PYTHON_LIBRARIES})

set(Boost_INCLUDE_DIR /Users/s23410fw/Documents/software/boost_1_66_0)
set(BOOST_LIBRARYDIR /Users/s23410fw/Documents/software/boost_1_66_0/libs/)
set(BOOST_ROOT /Users/s23410fw/Documents/software/boost_1_66_0)
find_package(Boost COMPONENTS filesystem system REQUIRED)
#find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
target_link_libraries(coinfinder ${Boost_LIBRARIES})
