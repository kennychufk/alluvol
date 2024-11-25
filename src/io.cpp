#include "alluvol/io.hpp"

#include <cstring>
#include <sstream>
#include <string>

#include "alluvol/geometry.hpp"

namespace alluvol {
namespace io {
void write_obj(openvdb::tools::VolumeToMesh const &mesher,
               char const *filename) {
  std::ofstream outfile(filename);
  for (unsigned int i = 0; i < mesher.pointListSize(); ++i) {
    outfile << "v"
            << " " << mesher.pointList()[i][0] << " "
            << mesher.pointList()[i][1] << " " << mesher.pointList()[i][2]
            << std::endl;
  }

  for (int pool_id = 0; pool_id < mesher.polygonPoolListSize(); ++pool_id) {
    openvdb::tools::PolygonPool const &pool = mesher.polygonPoolList()[pool_id];
    // write quad face
    for (unsigned int i = 0; i < pool.numQuads(); ++i)
      outfile << "f"
              << " " << pool.quad(i)[0] + 1 << " " << pool.quad(i)[1] + 1 << " "
              << pool.quad(i)[2] + 1 << " " << pool.quad(i)[3] + 1 << std::endl;
    // write triangle faces
    for (unsigned int i = 0; i < pool.numTriangles(); ++i)
      outfile << "f"
              << " " << pool.triangle(i)[0] + 1 << " "
              << pool.triangle(i)[1] + 1 << " " << pool.triangle(i)[2] + 1
              << std::endl;
  }

  outfile.close();
}

void read_obj(const char *filename, std::vector<openvdb::Vec3s> &points,
              std::vector<openvdb::Vec3I> &triangles,
              openvdb::Real resolution) {
  points.clear();
  triangles.clear();
  std::ifstream file_stream(filename);
  std::stringstream line_stream;
  std::string line;
  std::array<std::string, 4> tokens;
  std::stringstream face_entry_stream;
  std::array<std::string, 3> face_entry_tokens;
  openvdb::Vec3I face;
  int num_tokens;
  int face_token_id;
  file_stream.exceptions(std::ios_base::badbit);
  if (!file_stream.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    abort();
  }
  while (std::getline(file_stream, line)) {
    num_tokens = 0;
    line_stream.clear();
    line_stream.str(line);
    while (num_tokens < 4 &&
           std::getline(line_stream, tokens[num_tokens], ' ')) {
      ++num_tokens;
    }
    if (num_tokens == 0) continue;
    if (tokens[0] == "v") {
      points.emplace_back(std::stof(tokens[1]) * resolution,
                          std::stof(tokens[2]) * resolution,
                          std::stof(tokens[3]) * resolution);
    } else if (tokens[0] == "f") {
      for (U face_entry_id = 1; face_entry_id <= 3; ++face_entry_id) {
        face_token_id = 0;
        face_entry_stream.clear();
        face_entry_stream.str(tokens[face_entry_id]);
        while (face_token_id < 3 &&
               std::getline(face_entry_stream, face_entry_tokens[face_token_id],
                            '/')) {
          if (face_token_id == 0) {
            if (face_entry_id == 1)
              face.x() = std::stoul(face_entry_tokens[face_token_id]);
            if (face_entry_id == 2)
              face.y() = std::stoul(face_entry_tokens[face_token_id]);
            if (face_entry_id == 3)
              face.z() = std::stoul(face_entry_tokens[face_token_id]);
          }
          face_token_id += 1;
        }
      }
      triangles.push_back(face - 1);  // OBJ vertex: one-based indexing
    }
  }
}

std::vector<openvdb::math::Mat4d> read_pile(const char *filename,
                                            openvdb::Real voxel_size) {
  std::ifstream stream(filename, std::ios::binary);
  F3 x;
  F3 v;
  F4 q;
  F3 omega;
  std::vector<openvdb::math::Mat4d> matrices;
  if (!stream.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    abort();
  }
  while (true) {
    stream.read(reinterpret_cast<char *>(&x), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&v), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&q), sizeof(F4));
    stream.read(reinterpret_cast<char *>(&omega), sizeof(F3));
    matrices.push_back(get_transform_matrix(x, q, voxel_size));
    if (stream.peek() == std::ifstream::traits_type::eof()) break;
  }
  return matrices;
}

}  // namespace io
}  // namespace alluvol
