add_executable(allumesh
  main.cpp
)

target_link_libraries(allumesh
  PRIVATE
    lib_alluvol
)

add_executable(lscomp
  lscomp.cpp
)

target_link_libraries(lscomp
  PRIVATE
    lib_alluvol
)
