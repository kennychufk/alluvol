list(APPEND lib_alluvol_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/io.cpp
  ${CMAKE_CURRENT_LIST_DIR}/particle_list.cpp
  ${CMAKE_CURRENT_LIST_DIR}/level_set.cpp
  ${CMAKE_CURRENT_LIST_DIR}/geometry.cpp
)

add_library(lib_alluvol SHARED ${lib_alluvol_SOURCES})

target_link_libraries(lib_alluvol PUBLIC openvdb)


target_include_directories(lib_alluvol
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

pybind11_add_module(_alluvol pyalluvol.cpp)
target_compile_definitions(_alluvol
  PRIVATE VERSION_INFO=${ALLUVOL_VERSION_INFO})

target_link_libraries(_alluvol
  PUBLIC
    lib_alluvol
    pybind11::module
)
