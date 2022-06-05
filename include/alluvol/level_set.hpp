#ifndef ALLUVOL_LEVEL_SET_HPP
#define ALLUVOL_LEVEL_SET_HPP

#include <openvdb/openvdb.h>

#include "alluvol/data_type.hpp"

namespace alluvol {
openvdb::FloatGrid::Ptr create_liquid_level_set(
    std::vector<F3> const& particle_x, openvdb::Real radius,
    openvdb::Real voxel_size, openvdb::Real mask_min = 0,
    openvdb::Real mask_max = 0, openvdb::Real half_width = 3,
    U num_dilation = 3, U num_mean_curvature = 1, U num_erosion = 3);
openvdb::FloatGrid::Ptr create_mesh_level_set(const char* filename,
                                              openvdb::Real voxel_size,
                                              openvdb::Real scale = 1,
                                              openvdb::Real ex_band = 3,
                                              openvdb::Real in_band = 3);
openvdb::FloatGrid::Ptr transform_level_set(openvdb::FloatGrid::Ptr grid,
                                            openvdb::math::Mat4d const& m);
openvdb::FloatGrid::Ptr transform_level_set(openvdb::FloatGrid::Ptr grid,
                                            F3 const& x, F4 const& q, openvdb::Real voxel_size);
}  // namespace alluvol

#endif
