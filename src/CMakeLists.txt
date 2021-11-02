add_executable(allumesh
  main.cpp
)

target_link_libraries(allumesh
  PRIVATE
    OpenVDB::openvdb
)
