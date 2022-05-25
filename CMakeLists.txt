# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

set(PROJECT_NAME sim)
project(${PROJECT_NAME})

find_package(ROOT CONFIG REQUIRED)
set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(include)
include_directories(vendor/jsoncpp)

link_directories(${ROOT_LIBRARY_DIR})

file(GLOB SOURCES 
    ${PROJECT_SOURCE_DIR}/source/*.cpp 
    vendor/jsoncpp/jsoncpp.cpp)

# Create the main program using the library.
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} )