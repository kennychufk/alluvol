list(APPEND lib_alluvol_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/io.cpp
  ${CMAKE_CURRENT_LIST_DIR}/particle_list.cpp
  ${CMAKE_CURRENT_LIST_DIR}/level_set.cpp
  ${CMAKE_CURRENT_LIST_DIR}/geometry.cpp
)

add_library(lib_alluvol
  ${lib_alluvol_SOURCES}
)

target_link_libraries(lib_alluvol
  PUBLIC
    OpenVDB::openvdb
)


target_include_directories(lib_alluvol
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

add_library(_alluvol
  MODULE
    pyalluvol.cpp
)

set_target_properties(_alluvol
  PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    CUDA_VISIBILITY_PRESET "hidden"
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)

target_include_directories(_alluvol
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

target_link_libraries(_alluvol
  PUBLIC
    lib_alluvol
    pybind11::module
)
