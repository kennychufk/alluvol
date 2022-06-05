#include "alluvol/geometry.hpp"

namespace alluvol {
openvdb::math::Mat4d get_transform_matrix(const F3& x, const F4& q,
                                          openvdb::Real voxel_size) {
  return openvdb::math::Mat4d(
      static_cast<double>(1 - 2 * (q.y * q.y + q.z * q.z)),
      static_cast<double>(2 * (q.x * q.y + q.z * q.w)),
      static_cast<double>(2 * (q.x * q.z - q.y * q.w)), static_cast<double>(0),
      static_cast<double>(2 * (q.x * q.y - q.z * q.w)),
      static_cast<double>(1 - 2 * (q.x * q.x + q.z * q.z)),
      static_cast<double>(2 * (q.y * q.z + q.x * q.w)), static_cast<double>(0),
      static_cast<double>(2 * (q.x * q.z + q.y * q.w)),
      static_cast<double>(2 * (q.y * q.z - q.x * q.w)),
      static_cast<double>(1 - 2 * (q.x * q.x + q.y * q.y)),
      static_cast<double>(0), static_cast<double>(x.x / voxel_size),
      static_cast<double>(x.y / voxel_size),
      static_cast<double>(x.z / voxel_size), static_cast<double>(1));
}
}  // namespace alluvol
