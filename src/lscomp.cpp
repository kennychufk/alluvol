#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/math/Mat4.h>
#include <openvdb/math/Math.h>
#include <openvdb/math/Stats.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Statistics.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tree/LeafNode.h>

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <ratio>
#include <sstream>
#include <string>

#include "alluvol/inside_squared_accumulation.hpp"
#include "alluvol/io.hpp"
#include "alluvol/level_set.hpp"
#include "alluvol/particle_list.hpp"

using namespace alluvol;

int main(int argc, char *argv[]) {
  std::string radius_str1(argv[1]);
  std::vector<F3> particle_x1 = io::read_alu<1, F3>(argv[2]);
  openvdb::Real radius1 = std::stod(radius_str1);

  std::string radius_str2(argv[3]);
  std::vector<F3> particle_x2 = io::read_alu<1, F3>(argv[4]);
  openvdb::Real radius2 = std::stod(radius_str2);

  openvdb::Real voxel_size1 = radius1 / 1.001 * 2.0 / sqrt(3.0) / 2.0;
  openvdb::Real voxel_size2 = radius2 / 1.001 * 2.0 / sqrt(3.0) / 2.0;
  openvdb::FloatGrid::Ptr liquid_ls1 =
      create_liquid_level_set(particle_x1, radius1, voxel_size1);
  openvdb::FloatGrid::Ptr liquid_ls2 =
      create_liquid_level_set(particle_x2, radius2, voxel_size2);
  std::cout << "active voxel count 1 " << liquid_ls1->activeVoxelCount()
            << std::endl;
  std::cout << "active voxel count 2 " << liquid_ls2->activeVoxelCount()
            << std::endl;

  openvdb::FloatGrid::Ptr liquid_ls_intersection =
      openvdb::tools::csgIntersectionCopy(*liquid_ls1, *liquid_ls2);
  openvdb::tools::csgUnion(*liquid_ls1, *liquid_ls2, false);
  openvdb::tools::csgDifference(*liquid_ls1, *liquid_ls_intersection, true);
  std::cout << "active voxel count after csg " << liquid_ls1->activeVoxelCount()
            << std::endl;
  openvdb::math::Stats stat1 =
      openvdb::tools::statistics(liquid_ls1->cbeginValueOn());
  std::cout << "min: " << stat1.min() << " max: " << stat1.max()
            << " avg: " << stat1.avg() << " var: " << stat1.var() << std::endl;
  InsideSquaredAccumulation acc;
  openvdb::tools::accumulate(liquid_ls1->cbeginValueOn(), acc);
  std::cout << "sum: " << acc.sum << std::endl;

  openvdb::tools::VolumeToMesh mesher;
  mesher.operator()<openvdb::FloatGrid>(*liquid_ls1);
  io::write_obj(mesher, argv[5]);
}
