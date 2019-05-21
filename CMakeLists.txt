cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -I/usr/lib/python2.7/site-packages/numpy/core/include/ -g -w")
add_compile_options(-m64)
project(spatialtree)



add_executable(
    spatialtree
    spatialtree.cpp
)

# # find python
# execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE PYTHON_SEARCH_PATH)
# string(REGEX REPLACE "\n$" "" PYTHON_SEARCH_PATH "${PYTHON_SEARCH_PATH}")
# file(GLOB_RECURSE PYTHON_DY_LIBS ${PYTHON_SEARCH_PATH}/lib/libpython*.dylib ${PYTHON_SEARCH_PATH}/lib/libpython*.so)
# if (PYTHON_DY_LIBS)
#     list(GET PYTHON_DY_LIBS 0 PYTHON_LIBRARY)
#     message("-- Find shared libpython: ${PYTHON_LIBRARY}")
# else()
#     message(WARNING "Cannot find shared libpython, try find_package")
# endif()


# find_package(PythonLibs ${PYTHON_VERSION_STRING} 2.7)
find_package(PythonLibs 2.7)
target_include_directories(spatialtree PRIVATE ${PYTHON_INCLUDE_DIRS})

find_package( Boost 1.40 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# find_package(nlohmann_json 3.2.0 REQUIRED)

target_link_libraries(spatialtree ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

