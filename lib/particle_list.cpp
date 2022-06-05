#include "alluvol/particle_list.hpp"

namespace alluvol {
Particle::Particle() {
  p = openvdb::Vec3R(0);
  v = openvdb::Vec3R(0);
  r = 0;
}

Particle::Particle(const Particle &_p) {
  p = _p.p;
  v = _p.v;
  r = _p.r;
}

ParticleList::ParticleList(size_t size, openvdb::Real radius_scale,
                           openvdb::Real velocity_scale)
    : radius_scale_(radius_scale), velocity_scale_(velocity_scale) {
  particles_.resize(size);
}
ParticleList::ParticleList(openvdb::Real radius_scale,
                           openvdb::Real velocity_scale)
    : radius_scale_(radius_scale), velocity_scale_(velocity_scale) {
  particles_.resize(0);
}

void ParticleList::set(int i, const openvdb::Vec3R &p, const openvdb::Real &r,
                       const openvdb::Vec3R &v) {
  Particle pa;
  pa.p = p;
  pa.r = r;
  pa.v = v;
  particles_[i] = pa;
}

void ParticleList::set(const std::vector<F3> &x_list, openvdb::Real r) {
  for (U i = 0; i < x_list.size(); ++i) {
    F3 const &pos = x_list[i];
    set(i,
        openvdb::Vec3R(static_cast<openvdb::Real>(pos.x),
                       static_cast<openvdb::Real>(pos.y),
                       static_cast<openvdb::Real>(pos.z)),
        r);
  }
}

void ParticleList::add(const openvdb::Vec3R &p, const openvdb::Real &r,
                       const openvdb::Vec3R &v) {
  Particle pa;
  pa.p = p;
  pa.r = r;
  pa.v = v;
  particles_.push_back(pa);
}

openvdb::CoordBBox ParticleList::getBBox(const openvdb::GridBase &grid) {
  openvdb::CoordBBox bbox;
  openvdb::Coord &min = bbox.min(), &max = bbox.max();
  openvdb::Vec3R pos;
  openvdb::Real rad, invDx = 1.0 / grid.voxelSize()[0];
  for (size_t n = 0, e = this->size(); n < e; ++n) {
    this->getPosRad(n, pos, rad);
    const openvdb::Vec3d xyz = grid.worldToIndex(pos);
    const openvdb::Real r = rad * invDx;
    for (int i = 0; i < 3; ++i) {
      min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
      max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil(xyz[i] + r));
    }
  }
  return bbox;
}
openvdb::Vec3R ParticleList::pos(int n) const { return particles_[n].p; }
openvdb::Vec3R ParticleList::vel(int n) const {
  return velocity_scale_ * particles_[n].v;
}
openvdb::Real ParticleList::radius(int n) const {
  return radius_scale_ * particles_[n].r;
}

size_t ParticleList::size() const { return particles_.size(); }

/// Get the world space position of n'th particle.
/// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
void ParticleList::getPos(size_t n, openvdb::Vec3R &pos) const {
  pos = particles_[n].p;
}

void ParticleList::getPosRad(size_t n, openvdb::Vec3R &pos,
                             openvdb::Real &rad) const {
  pos = particles_[n].p;
  rad = radius_scale_ * particles_[n].r;
}
void ParticleList::getPosRadVel(size_t n, openvdb::Vec3R &pos,
                                openvdb::Real &rad, openvdb::Vec3R &vel) const {
  pos = particles_[n].p;
  rad = radius_scale_ * particles_[n].r;
  vel = velocity_scale_ * particles_[n].v;
}
// The method below is only required for attribute transfer
void ParticleList::getAtt(size_t n, openvdb::Vec3s &att) const {
  att = particles_[n].v;
}

}  // namespace alluvol
