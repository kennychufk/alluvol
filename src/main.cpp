#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/math/Math.h>
#include <openvdb/openvdb.h>
#include <openvdb/tree/LeafNode.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "openvdb/tools/LevelSetFilter.h"
#include "openvdb/tools/ParticlesToLevelSet.h"
#include "openvdb/tools/VolumeToMesh.h"

struct float3 {
  float x;
  float y;
  float z;
};
using U = unsigned int;
using F = float;
using F3 = float3;

template <U D, typename M>
std::vector<M> read_file(const char *filename) {
  std::ifstream stream(filename, std::ios::binary);
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

class MyParticleList {
 protected:
  struct MyParticle {
    openvdb::Vec3R p, v;
    openvdb::Real r;
    MyParticle() {
      p = openvdb::Vec3R(0);
      v = openvdb::Vec3R(0);
      r = 0;
    }
    MyParticle(const MyParticle &_p) {
      p = _p.p;
      v = _p.v;
      r = _p.r;
    }
  };
  openvdb::Real mRadiusScale;
  openvdb::Real mVelocityScale;
  std::vector<MyParticle> mParticleList;

 public:
  typedef openvdb::Vec3R PosType;

  MyParticleList(size_t size, openvdb::Real rScale = 1,
                 openvdb::Real vScale = 1)
      : mRadiusScale(rScale), mVelocityScale(vScale) {
    mParticleList.resize(size);
  }
  MyParticleList(openvdb::Real rScale = 1, openvdb::Real vScale = 1)
      : mRadiusScale(rScale), mVelocityScale(vScale) {
    mParticleList.resize(0);
  }
  // void free() { mParticleList.swap(std::vector<MyParticle>()); }
  void set(int i, const openvdb::Vec3R &p, const openvdb::Real &r,
           const openvdb::Vec3R &v = openvdb::Vec3R(0, 0, 0)) {
    MyParticle pa;
    pa.p = p;
    pa.r = r;
    pa.v = v;
    mParticleList[i] = pa;
  }
  void add(const openvdb::Vec3R &p, const openvdb::Real &r,
           const openvdb::Vec3R &v = openvdb::Vec3R(0, 0, 0)) {
    MyParticle pa;
    pa.p = p;
    pa.r = r;
    pa.v = v;
    mParticleList.push_back(pa);
  }
  /// @return coordinate bbox in the space of the specified transfrom
  openvdb::CoordBBox getBBox(const openvdb::GridBase &grid) {
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
  // typedef int AttributeType;
  // The methods below are only required for the unit-tests
  openvdb::Vec3R pos(int n) const { return mParticleList[n].p; }
  openvdb::Vec3R vel(int n) const {
    return mVelocityScale * mParticleList[n].v;
  }
  openvdb::Real radius(int n) const {
    return mRadiusScale * mParticleList[n].r;
  }

  //////////////////////////////////////////////////////////////////////////////
  /// The methods below are the only ones required by tools::ParticleToLevelSet
  /// @note We return by value since the radius and velocities are modified
  /// by the scaling factors! Also these methods are all assumed to
  /// be thread-safe.

  /// Return the total number of particles in list.
  ///  Always required!
  size_t size() const { return mParticleList.size(); }

  /// Get the world space position of n'th particle.
  /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
  void getPos(size_t n, openvdb::Vec3R &pos) const { pos = mParticleList[n].p; }

  void getPosRad(size_t n, openvdb::Vec3R &pos, openvdb::Real &rad) const {
    pos = mParticleList[n].p;
    rad = mRadiusScale * mParticleList[n].r;
  }
  void getPosRadVel(size_t n, openvdb::Vec3R &pos, openvdb::Real &rad,
                    openvdb::Vec3R &vel) const {
    pos = mParticleList[n].p;
    rad = mRadiusScale * mParticleList[n].r;
    vel = mVelocityScale * mParticleList[n].v;
  }
  // The method below is only required for attribute transfer
  void getAtt(size_t n, openvdb::Index32 &att) const {
    att = openvdb::Index32(n);
  }
};
// VDB will take the particles from my simulation code and build a mesh from
// them.
int main(int argc, char *argv[]) {
  std::vector<F3> particle_x = read_file<1, F3>(argv[1]);
  int frame = 0;
  openvdb::Real radius = 1e-3;
  MyParticleList pa(particle_x.size(), 1, 1);
  for (U i = 0; i < particle_x.size(); ++i) {
    F3 const &pos = particle_x[i];
    pa.set(i,
           openvdb::Vec3R(static_cast<openvdb::Real>(pos.x),
                          static_cast<openvdb::Real>(pos.y),
                          static_cast<openvdb::Real>(pos.z)),
           radius);
  }
  openvdb::Real voxelSize = radius / 1.001 * 2.0 / sqrt(3.0) / 2.0,
                halfWidth = 2.0;
  openvdb::FloatGrid::Ptr ls =
      openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid, openvdb::Index32>
      raster(*ls);

  raster.setGrainSize(1);  // a value of zero disables threading
  raster.rasterizeSpheres(pa);
  raster.finalize(true);
  openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*ls);
  std::cout << "filtering" << std::endl;
  filter.dilate(5);
  std::cout << "dilated" << std::endl;
  for (int i = 0; i < 1; ++i) {
    filter.meanCurvature();
  }
  std::cout << "curved" << std::endl;
  filter.erode();
  std::cout << "eroded" << std::endl;
  filter.meanCurvature();
  std::cout << "final smooth" << std::endl;
  // filter.gaussian(radius*8);
  // std::cout<<"gaussian"<<std::endl;

  openvdb::CoordBBox bbox = pa.getBBox(*ls);
  std::cout << bbox.min() << std::endl;
  std::cout << bbox.max() << std::endl;

  // openvdb::tools::volumeToMesh(*ls, points, triangles, quads, 0.0, 0.5);

  openvdb::tools::VolumeToMesh mesher;
  mesher.operator()<openvdb::FloatGrid>(filter.grid());
  std::cout << "No. of polygon pool = " << mesher.polygonPoolListSize()
            << std::endl;
  printf("meshing done\n");

  std::ofstream outfile(argv[2]);

  // write vertices
  for (unsigned int i = 0; i < mesher.pointListSize(); ++i) {
    outfile << "v"
            << " " << mesher.pointList()[i][0] << " "
            << mesher.pointList()[i][1] << " " << mesher.pointList()[i][2]
            << std::endl;
  }

  for (int pool_id = 0; pool_id < mesher.polygonPoolListSize(); ++pool_id) {
    openvdb::tools::PolygonPool const &pool = mesher.polygonPoolList()[pool_id];
    std::cout << "num of triangles = " << pool.numTriangles() << "< "
              << pool.numQuads() << std::endl;
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

  // // write vertices
  // for (unsigned int i = 0; i < points.size(); ++i)
  //   outfile << "v"
  //           << " " << points[i][0] << " " << points[i][1] << " " <<
  //           points[i][2]
  //           << std::endl;
  // // write quad face
  // for (unsigned int i = 0; i < quads.size(); ++i)
  //   outfile << "f"
  //           << " " << quads[i][0] + 1 << " " << quads[i][1] + 1 << " "
  //           << quads[i][2] + 1 << " " << quads[i][3] + 1 << std::endl;
  // // write triangle faces
  // for (unsigned int i = 0; i < triangles.size(); ++i)
  //   outfile << "f"
  //           << " " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << "
  //           "
  //           << triangles[i][2] + 1 << std::endl;
  outfile.close();
  // // points.swap(std::vector<openvdb::Vec3s>());
  // // quads.swap(std::vector<openvdb::Vec4I>());
  // // triangles.swap(std::vector<openvdb::Vec3I>());
}
