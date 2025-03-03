cmake_minimum_required(VERSION 3.10)

project(cpp_project)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(ex01 src/ex01.cpp)
target_link_libraries(ex01 Eigen3::Eigen)

add_executable(ex02 src/ex02.cpp)
target_link_libraries(ex02 Eigen3::Eigen)

add_executable(ex03 src/ex03.cpp)
target_link_libraries(ex03 Eigen3::Eigen)

add_executable(ex04 src/ex04.cpp)
target_link_libraries(ex04 Eigen3::Eigen)

add_executable(ex05 src/ex05.cpp)
target_link_libraries(ex05 Eigen3::Eigen)

add_executable(ex06 src/ex06.cpp)
target_link_libraries(ex06 Eigen3::Eigen)

add_executable(ex07 src/ex07.cpp)
target_link_libraries(ex07 Eigen3::Eigen)

add_executable(ex08 src/ex08.cpp)
target_link_libraries(ex08 Eigen3::Eigen)

add_executable(ex09 src/ex09.cpp)
target_link_libraries(ex09 Eigen3::Eigen)

add_executable(ex10 src/ex10.cpp)
target_link_libraries(ex10 Eigen3::Eigen)

add_executable(ex11 src/ex11.cpp)
target_link_libraries(ex11 Eigen3::Eigen)

add_executable(ex12 src/ex12.cpp)
target_link_libraries(ex12 Eigen3::Eigen)

add_executable(ex13 src/ex13.cpp)
target_link_libraries(ex13 Eigen3::Eigen)

add_executable(ex14 src/ex14.cpp)
target_link_libraries(ex14 Eigen3::Eigen)

add_executable(ex15 src/ex15.cpp)
target_link_libraries(ex15 Eigen3::Eigen)