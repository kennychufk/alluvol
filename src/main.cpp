#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Statistics.h>
#include <openvdb/tree/LeafNode.h>

#include <chrono>
#include <iomanip>
#include <ratio>
#include <string>

#include "alluvol/io.hpp"
#include "alluvol/level_set.hpp"
#include "alluvol/particle_list.hpp"

using namespace alluvol;

bool is_file_exist(const char *fileName) {
  std::ifstream infile(fileName);
  return infile.good();
}
// VDB will take the particles from my simulation code and build a mesh from
// them.
int main(int argc, char *argv[]) {
  auto t0 = std::chrono::high_resolution_clock::now();
  std::vector<F3> particle_x = io::read_alu<1, F3>(argv[2]);
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> fp_ms = t1 - t0;
  std::cout << "read alu: " << fp_ms.count() << std::endl;
  std::string radius_str(argv[1]);
  openvdb::Real radius = std::stod(radius_str);
  t0 = std::chrono::high_resolution_clock::now();

  openvdb::Real voxel_size = radius / 1.001 * 2.0 / sqrt(3.0) / 2.0;
  std::string alu_filename(argv[2]);
  std::string vdb_filename =
      alu_filename.substr(0, alu_filename.size() - 4).append(".vdb");
  openvdb::FloatGrid::Ptr filtered_grid = nullptr;
  if (is_file_exist(vdb_filename.c_str())) {
    openvdb::initialize();
    std::cout << "VDB file " << vdb_filename << " exists" << std::endl;
    openvdb::io::File file(vdb_filename);
    file.open();
    openvdb::GridBase::Ptr base_grid =
        file.readGrid(file.beginName().gridName());
    file.close();
    filtered_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
  } else {
    filtered_grid =
        create_liquid_level_set(particle_x, radius, voxel_size, 0.2, 0.3);
  }

  openvdb::FloatGrid::Ptr buoy_ls;
  bool has_buoy = (std::strlen(argv[5]) > 1 || argv[5][0] != '0');
  if (has_buoy) {
    buoy_ls = create_mesh_level_set(argv[5], voxel_size);
  }

  openvdb::FloatGrid::Ptr container_ls =
      create_mesh_level_set(argv[4], voxel_size);

  openvdb::FloatGrid::Ptr agitator_ls;
  bool has_agitator = (std::strlen(argv[6]) > 1 || argv[6][0] != '0');
  if (has_agitator) {
    std::string agitator_scale_str(argv[7]);
    double scale = std::stod(agitator_scale_str);
    agitator_ls = create_mesh_level_set(argv[6], voxel_size, scale);
  }

  std::vector<openvdb::math::Mat4d> matrices =
      io::read_pile(argv[3], voxel_size);
  openvdb::math::Mat4d &agitator_matrix = matrices[matrices.size() - 1];
  int num_buoys = matrices.size() - 2;
  for (int i = 0; i < num_buoys + 1 /* also subtract agitator*/; ++i) {
    openvdb::FloatGrid::Ptr transformed_obj_ls = transform_level_set(
        (i == num_buoys ? agitator_ls : buoy_ls), matrices[i + 1]);

    t0 = std::chrono::high_resolution_clock::now();
    openvdb::tools::csgDifference(*filtered_grid, *transformed_obj_ls);
    t1 = std::chrono::high_resolution_clock::now();
    fp_ms = t1 - t0;
    std::cout << "csgDifference: " << fp_ms.count() << std::endl;
  }

  {
    t0 = std::chrono::high_resolution_clock::now();
    openvdb::tools::csgIntersection(*filtered_grid, *container_ls);
    t1 = std::chrono::high_resolution_clock::now();
    fp_ms = t1 - t0;
    std::cout << "csgIntersection: " << fp_ms.count() << std::endl;
  }

  openvdb::tools::VolumeToMesh mesher;
  t0 = std::chrono::high_resolution_clock::now();
  mesher.operator()<openvdb::FloatGrid>(*filtered_grid);
  t1 = std::chrono::high_resolution_clock::now();
  fp_ms = t1 - t0;
  std::cout << "mesher: " << fp_ms.count() << std::endl;

  // write vertices
  t0 = std::chrono::high_resolution_clock::now();
  io::write_obj(mesher, argv[8]);
  t1 = std::chrono::high_resolution_clock::now();
  fp_ms = t1 - t0;
  std::cout << "export: " << fp_ms.count() << std::endl;
}
