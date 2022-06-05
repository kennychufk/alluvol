#ifndef ALLUVOL_PARTICLE_LIST_HPP
#define ALLUVOL_PARTICLE_LIST_HPP

#include <openvdb/math/Vec3.h>
#include <openvdb/openvdb.h>

#include "alluvol/data_type.hpp"

namespace alluvol {

struct Particle {
  openvdb::Vec3R p, v;
  openvdb::Real r;
  Particle();
  Particle(const Particle &_p);
};

class ParticleList {
 protected:
  openvdb::Real radius_scale_;
  openvdb::Real velocity_scale_;
  std::vector<Particle> particles_;

 public:
  typedef openvdb::Vec3R PosType;
  ParticleList(size_t size, openvdb::Real radius_scale = 1,
               openvdb::Real velocity_scale = 1);
  ParticleList(openvdb::Real radius_scale = 1,
               openvdb::Real velocity_scale = 1);
  void set(int i, const openvdb::Vec3R &p, const openvdb::Real &r,
           const openvdb::Vec3R &v = openvdb::Vec3R(0, 0, 0));
  void set(const std::vector<F3> &x_list, openvdb::Real r);
  void add(const openvdb::Vec3R &p, const openvdb::Real &r,
           const openvdb::Vec3R &v = openvdb::Vec3R(0, 0, 0));
  /// @return coordinate bbox in the space of the specified transfrom
  openvdb::CoordBBox getBBox(const openvdb::GridBase &grid);
  // typedef int AttributeType;
  // The methods below are only required for the unit-tests
  openvdb::Vec3R pos(int n) const;
  openvdb::Vec3R vel(int n) const;
  openvdb::Real radius(int n) const;

  //////////////////////////////////////////////////////////////////////////////
  /// The methods below are the only ones required by tools::ParticleToLevelSet
  /// @note We return by value since the radius and velocities are modified
  /// by the scaling factors! Also these methods are all assumed to
  /// be thread-safe.

  /// Return the total number of particles in list.
  ///  Always required!
  size_t size() const;

  /// Get the world space position of n'th particle.
  /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
  void getPos(size_t n, openvdb::Vec3R &pos) const;

  void getPosRad(size_t n, openvdb::Vec3R &pos, openvdb::Real &rad) const;
  void getPosRadVel(size_t n, openvdb::Vec3R &pos, openvdb::Real &rad,
                    openvdb::Vec3R &vel) const;
  // The method below is only required for attribute transfer
  void getAtt(size_t n, openvdb::Vec3s &att) const;
};
}  // namespace alluvol

#endif
