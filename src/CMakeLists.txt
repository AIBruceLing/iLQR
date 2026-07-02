cmake_minimum_required(VERSION 3.8)
project(iLQR_standalone)

set(CMAKE_CXX_STANDARD 17)

# Eigen3
find_package(Eigen3 REQUIRED)

# Python3 (matplotlibcpp 依赖)
find_package(Python3 REQUIRED COMPONENTS Development NumPy)

add_executable(iLQR iLQR.cpp)

target_include_directories(iLQR PRIVATE
    ${EIGEN3_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}  # matplotlibcpp.h 所在目录
)

target_link_libraries(iLQR PRIVATE
    ${Python3_LIBRARIES}
)
