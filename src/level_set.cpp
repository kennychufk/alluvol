#include "alluvol/level_set.hpp"

#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

#include "alluvol/geometry.hpp"
#include "alluvol/io.hpp"
#include "alluvol/particle_list.hpp"

namespace alluvol {
openvdb::FloatGrid::Ptr create_liquid_level_set(
    std::vector<F3> const& particle_x, openvdb::Real radius,
    openvdb::Real voxel_size, openvdb::Real mask_min, openvdb::Real mask_max,
    openvdb::Real half_width, U num_dilation, U num_mean_curvature,
    U num_erosion) {
  ParticleList particle_list(particle_x.size(), 1, 1);
  particle_list.set(particle_x, radius);
  openvdb::FloatGrid::Ptr ls =
      openvdb::createLevelSet<openvdb::FloatGrid>(voxel_size, half_width);
  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid, openvdb::Vec3s>
      raster(*ls);

  raster.setGrainSize(1);  // a value of zero disables threading
  raster.rasterizeSpheres(particle_list);
  raster.finalize(true);
  openvdb::tools::LevelSetFilter<openvdb::FloatGrid, openvdb::FloatGrid> filter(
      *ls);
  filter.dilate(num_dilation);
  if (mask_min != 0 || mask_max != 0) {
    filter.setMaskRange(0.2, 0.3);
    filter.invertMask();
  }
  for (U i = 0; i < num_mean_curvature; ++i) {
    filter.meanCurvature();
  }
  filter.erode(num_erosion);
  // filter.meanCurvature();
  filter.gaussian();
  return ls;
}

openvdb::FloatGrid::Ptr create_mesh_level_set(const char* filename,
                                              openvdb::Real voxel_size,
                                              openvdb::Real scale,
                                              openvdb::Real ex_band,
                                              openvdb::Real in_band) {
  const int kConversionFlags = 0;
  openvdb::math::Transform::Ptr transform =
      openvdb::math::Transform::createLinearTransform(voxel_size);
  std::vector<openvdb::Vec3s> point_list;
  std::vector<openvdb::Vec3I> prim_list;

  io::read_obj(filename, point_list, prim_list, scale / voxel_size);
  openvdb::tools::QuadAndTriangleDataAdapter<openvdb::Vec3s, openvdb::Vec3I>
      mesh(point_list, prim_list);
  return openvdb::tools::meshToVolume<openvdb::FloatGrid>(
      mesh, *transform, ex_band, in_band, kConversionFlags, nullptr);
}

openvdb::FloatGrid::Ptr transform_level_set(openvdb::FloatGrid::Ptr grid,
                                            openvdb::math::Mat4d const& m) {
  openvdb::FloatGrid::Ptr transformed_grid = openvdb::FloatGrid::create();
  openvdb::tools::GridTransformer transformer(m);
  transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(
      *grid, *transformed_grid);
  return transformed_grid;
}

openvdb::FloatGrid::Ptr transform_level_set(openvdb::FloatGrid::Ptr grid,
                                            F3 const& x, F4 const& q,
                                            openvdb::Real voxel_size) {
  openvdb::math::Mat4d m = get_transform_matrix(x, q, voxel_size);
  return transform_level_set(grid, m);
}
}  // namespace alluvol
