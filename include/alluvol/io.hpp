#ifndef ALLUVOL_IO_HPP
#define ALLUVOL_IO_HPP

#include <openvdb/math/Mat4.h>
#include <openvdb/math/Math.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/tools/VolumeToMesh.h>

#include <fstream>

#include "alluvol/data_type.hpp"

namespace alluvol {
namespace io {
template <U D, typename M>
std::vector<M> read_alu(const char *filename) {
  std::ifstream stream(filename, std::ios::binary);
  if (!stream.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    abort();
  }
  U linear_shape = 1;
  for (U i = 0; i < D; ++i) {
    U shape_item;
    stream.read(reinterpret_cast<char *>(&shape_item), sizeof(U));
    if (shape_item == 0) {
      std::cerr << "Shape mismatch when reading " << filename << std::endl;
      abort();
    }
    linear_shape *= shape_item;
  }
  U end_of_shape;
  stream.read(reinterpret_cast<char *>(&end_of_shape), sizeof(U));
  if (end_of_shape != 0) {
    std::cerr << "Dimension mismatch when reading " << filename << std::endl;
    abort();
  }
  U num_primitives_per_unit;
  stream.read(reinterpret_cast<char *>(&num_primitives_per_unit), sizeof(U));
  char type_label;
  stream.read(reinterpret_cast<char *>(&type_label), sizeof(char));
  U num_bytes = linear_shape * sizeof(M);
  std::vector<M> host_buffer(linear_shape);

  stream.read(reinterpret_cast<char *>(host_buffer.data()), num_bytes);
  return host_buffer;
}

void write_obj(openvdb::tools::VolumeToMesh const &mesher,
               char const *filename);

void read_obj(const char *filename, std::vector<openvdb::Vec3s> &points,
              std::vector<openvdb::Vec3I> &triangles, openvdb::Real resolution);

std::vector<openvdb::math::Mat4d> read_pile(const char *filename,
                                            openvdb::Real voxel_size);
}  // namespace io
}  // namespace alluvol
#endif
