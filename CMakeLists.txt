cmake_minimum_required(VERSION 3.6)
project(coinfinder)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Boost_USE_STATIC_LIBS ON)

set(SOURCE_FILES coinfind-code/main.cpp coinfind-code/id_lookup.cpp coinfind-code/id_lookup.h coinfind-code/dataset.cpp coinfind-code/dataset.h coinfind-code/coincidence.cpp coinfind-code/coincidence.h coinfind-code/connectivity.cpp coinfind-code/connectivity.h coinfind-code/binomial_test.cpp coinfind-code/binomial_test.h coinfind-code/significance.cpp coinfind-code/significance.h coinfind-code/constants.h coinfind-code/parameters.cpp coinfind-code/parameters.h coinfind-code/elements.cpp coinfind-code/elements.h coinfind-code/test_cases.cpp coinfind-code/test_cases.h coinfind-code/phylomax.py coinfind-code/common_ancestor.py coinfind-code/bugfix.h coinfind-code/bugfix.cpp coinfind-code/lineage.cpp coinfind-code/lineage.h coinfind-code/gexf.cpp coinfind-code/gexf.h coinfind-code/network.cpp coinfind-code/network.h)
add_executable(coinfinder ${SOURCE_FILES})

target_compile_features(coinfinder PRIVATE cxx_std_14)

find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonLibs 3.6 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(coinfinder ${PYTHON_LIBRARIES})

find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
target_link_libraries(coinfinder ${Boost_LIBRARIES})

#dlopen, dlerror, dlsym undefined reference error fix
#forkpty undefined reference error fix
target_link_libraries(coinfinder ${CMAKE_DL_LIBS})
target_link_libraries(coinfinder "-lutil")
set(CMAKE_MODULE_LINKER_FLAGS "-lutil")
find_package(Threads)
target_link_libraries(coinfinder ${CMAKE_THREAD_LIBS_INIT})

# Find OpenMP
if(APPLE)
#    SET(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
#    SET(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
#    SET(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
#    SET(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

find_package(OPENMP)

if (OPENMP_FOUND)
    include_directories("${OPENMP_INCLUDES}")
    #include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
    link_directories("${OPENMP_LIBRARIES}")
    #link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
