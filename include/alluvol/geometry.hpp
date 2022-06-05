#ifndef ALLUVOL_GEOMETRY_HPP
#define ALLUVOL_GEOMETRY_HPP

#include <openvdb/math/Mat4.h>
#include <openvdb/openvdb.h>

#include "alluvol/data_type.hpp"

namespace alluvol {
openvdb::math::Mat4d get_transform_matrix(const F3& x, const F4& q,
                                          openvdb::Real voxel_size);
}  // namespace alluvol

#endif
