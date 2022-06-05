cmake_minimum_required (VERSION 3.20.0)
project(alluvol
  LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenVDB REQUIRED)
find_package(pybind11 REQUIRED)

add_subdirectory(lib)

add_executable(allumesh
  src/main.cpp
)

target_link_libraries(allumesh
  PRIVATE
    lib_alluvol
)

add_executable(lscomp
  src/lscomp.cpp
)

target_link_libraries(lscomp
  PRIVATE
    lib_alluvol
)

include(GNUInstallDirs)
install(TARGETS allumesh
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kenny Chu")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "openvdb")
include(CPack)
