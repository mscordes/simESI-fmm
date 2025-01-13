cmake_minimum_required(VERSION 3.20)

project(
    simESI-fmm 	
    VERSION 2020.1
    LANGUAGES CXX
    HOMEPAGE_URL "https://github.com/mscordes/simESI-fmm")

set(CMAKE_CXX_STANDARD 20)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(UNIX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
    message(FATAL_ERROR "Unsupported platform")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

file(GLOB_RECURSE SRC_FILES "src/*.cpp")

message(STATUS "Source files: ${SRC_FILES}")

add_executable(simESI-fmm src/main.cpp ${SRC_FILES})

target_include_directories(simESI-fmm PRIVATE src)

set_target_properties(simESI-fmm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}
)
