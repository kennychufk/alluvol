cmake_minimum_required(VERSION 3.20.0)
project(allumesh
  LANGUAGES CXX C
)

find_package(OpenVDB REQUIRED)

add_subdirectory(src)

include(GNUInstallDirs)
install(TARGETS allumesh
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kenny Chu")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "openvdb")
include(CPack)
