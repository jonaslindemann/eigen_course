cmake_minimum_required(VERSION 3.10)

# Auto-detect vcpkg toolchain unless the user already set one.
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
	set(_VCPKG_TOOLCHAIN "")

	# 1) Local vcpkg checkout (e.g. repo submodule)
	set(_VCPKG_LOCAL "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
	if(EXISTS "${_VCPKG_LOCAL}")
		set(_VCPKG_TOOLCHAIN "${_VCPKG_LOCAL}")
	endif()

	# 2) vcpkg on PATH
	if(NOT _VCPKG_TOOLCHAIN)
		find_program(VCPKG_EXECUTABLE vcpkg)
		if(VCPKG_EXECUTABLE)
			get_filename_component(_VCPKG_ROOT "${VCPKG_EXECUTABLE}" DIRECTORY)
			set(_VCPKG_FROM_PATH "${_VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
			if(EXISTS "${_VCPKG_FROM_PATH}")
				set(_VCPKG_TOOLCHAIN "${_VCPKG_FROM_PATH}")
			endif()
		endif()
	endif()

	# 3) Environment variables (optional fallback)
	if(NOT _VCPKG_TOOLCHAIN)
		if(DEFINED ENV{VCPKG_ROOT})
			set(_VCPKG_TOOLCHAIN "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
		elseif(DEFINED ENV{VCPKG_INSTALLATION_ROOT})
			set(_VCPKG_TOOLCHAIN "$ENV{VCPKG_INSTALLATION_ROOT}/scripts/buildsystems/vcpkg.cmake")
		endif()
	endif()

	if(_VCPKG_TOOLCHAIN AND EXISTS "${_VCPKG_TOOLCHAIN}")
		set(CMAKE_TOOLCHAIN_FILE "${_VCPKG_TOOLCHAIN}" CACHE STRING "Vcpkg toolchain file")
	endif()
endif()

project(cpp_project)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

function(add_eigen_exe target)
	add_executable(${target} ${ARGN})
	target_link_libraries(${target} Eigen3::Eigen)
	target_include_directories(${target} PRIVATE ${PROJECT_SOURCE_DIR}/include)
endfunction()

set(SIMPLE_EXES
	ex01
	ex02
	ex03
	ex04
	ex05
	ex06
	ex07
	ex08
	ex09
	ex10
	ex11
	ex12
	ex13
	ex14
	ex15
	matrix2_eig
	eigen_c_arrays1
)

foreach(exe IN LISTS SIMPLE_EXES)
	add_eigen_exe(${exe} src/${exe}.cpp)
endforeach()

set(CALFEM_EXES
	live_coding
	matrix1_eig
	matrix3_eig
	matrix4_eig
	matrix5_eig
)

foreach(exe IN LISTS CALFEM_EXES)
	add_eigen_exe(${exe} src/${exe}.cpp src/calfem_eig.cpp)
endforeach()

