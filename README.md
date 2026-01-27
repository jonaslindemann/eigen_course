# eigen_course
Course material for the course Array Computing with Eigen.

## Build instructions

### Prerequisites
- CMake 3.10 or newer
- A C++20 compiler (MSVC, clang, or GCC)
- vcpkg with Eigen installed

### Configure vcpkg
The build auto-detects vcpkg using one of these methods (in order):
1. A local vcpkg checkout at ./vcpkg (e.g., a submodule)
2. A vcpkg executable available on your PATH
3. The VCPKG_ROOT or VCPKG_INSTALLATION_ROOT environment variable

### Configure and build
1. Configure: cmake -S . -B build
2. Build: cmake --build build

Executables are placed under the build directory.
