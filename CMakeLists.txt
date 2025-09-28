cmake_minimum_required(VERSION 3.10)
project(spatialtree)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Architecture
add_compile_options(-m64)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Main executable
add_executable(
    spatialtree
    spatialtree.cpp
)

# Optional: Python support for visualization
option(ENABLE_PYTHON "Enable Python support for visualization" OFF)
if(ENABLE_PYTHON)
    find_package(PythonLibs 2.7)
    if(PYTHONLIBS_FOUND)
        target_include_directories(spatialtree PRIVATE ${PYTHON_INCLUDE_DIRS})
        target_link_libraries(spatialtree ${PYTHON_LIBRARIES})
        target_compile_definitions(spatialtree PRIVATE PYTHON_ENABLED)

        # NumPy include path (may need adjustment based on system)
        execute_process(
            COMMAND python -c "import numpy; print(numpy.get_include())"
            OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NUMPY_INCLUDE_DIR)
            target_include_directories(spatialtree PRIVATE ${NUMPY_INCLUDE_DIR})
        endif()
    else()
        message(WARNING "Python libraries not found, visualization features disabled")
    endif()
endif()

# Optional: Boost support for JSON parsing
option(ENABLE_BOOST "Enable Boost for JSON parsing" OFF)
if(ENABLE_BOOST)
    find_package(Boost 1.40 COMPONENTS program_options)
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIR})
        target_link_libraries(spatialtree ${Boost_LIBRARIES})
        target_compile_definitions(spatialtree PRIVATE BOOST_ENABLED)
    else()
        message(WARNING "Boost not found, JSON parsing features may be limited")
    endif()
endif()

# Test executable (optional)
option(BUILD_TESTS "Build test programs" OFF)
if(BUILD_TESTS)
    add_executable(test_spatial tests/test.cpp)
    target_include_directories(test_spatial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

