#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/math/Mat4.h>
#include <openvdb/math/Math.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tree/LeafNode.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "openvdb/tools/Composite.h"
#include "openvdb/tools/LevelSetFilter.h"
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/ParticlesToLevelSet.h"
#include "openvdb/tools/VolumeToMesh.h"

struct float3 {
  float x;
  float y;
  float z;
};
struct float4 {
  float x;
  float y;
  float z;
  float w;
};
using U = unsigned int;
using F = float;
using F3 = float3;
using F4 = float4;

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

void load_obj(const char *filename, std::vector<openvdb::Vec3s> &points,
              std::vector<openvdb::Vec3I> &triangles, float resolution) {
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
                                            float voxel_size) {
  std::ifstream stream(filename, std::ios::binary);
  F3 x;
  F3 v;
  F4 q;
  F3 omega;
  std::vector<openvdb::math::Mat4d> matrices;
  while (true) {
    stream.read(reinterpret_cast<char *>(&x), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&v), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&q), sizeof(F4));
    stream.read(reinterpret_cast<char *>(&omega), sizeof(F3));
    std::cout << x.x << ", " << x.y << ", " << x.z << std::endl;
    std::cout << q.x << ", " << q.y << ", " << q.z << ", " << q.w << std::endl;

    // openvdb::math::Mat4d mat(
    //     static_cast<double>(1 - 2 * (q.y * q.y + q.z * q.z)),
    //     static_cast<double>(2 * (q.x * q.y + q.z * q.w)),
    //     static_cast<double>(2 * (q.x * q.z - q.y * q.w)),
    //     static_cast<double>(0),
    //     static_cast<double>(2 * (q.x * q.y - q.z * q.w)),
    //     static_cast<double>(1 - 2 * (q.x * q.x + q.z * q.z)),
    //     static_cast<double>(2 * (q.y * q.z + q.x * q.w)),
    //     static_cast<double>(0),
    //     static_cast<double>(2 * (q.x * q.z + q.y * q.w)),
    //     static_cast<double>(2 * (q.y * q.z - q.x * q.w)),
    //     static_cast<double>(1 - 2 * (q.x * q.x + q.y * q.y)),
    //     static_cast<double>(0), static_cast<double>(0),
    //     static_cast<double>(0),
    //     static_cast<double>(0), static_cast<double>(1));
    // openvdb::math::Mat4d mat = openvdb::math::Mat4d::identity();
    // mat.postTranslate(openvdb::math::Vec3d(x.x/voxel_size, x.y/voxel_size,
    // x.z/voxel_size));
    // mat.preScale(openvdb::math::Vec3d(1.0/voxel_size, 1.0/voxel_size, 1.0/voxel_size));
    openvdb::math::Mat4d mat(
        static_cast<double>(1 - 2 * (q.y * q.y + q.z * q.z)),
        static_cast<double>(2 * (q.x * q.y + q.z * q.w)),
        static_cast<double>(2 * (q.x * q.z - q.y * q.w)),
        static_cast<double>(0),
        static_cast<double>(2 * (q.x * q.y - q.z * q.w)),
        static_cast<double>(1 - 2 * (q.x * q.x + q.z * q.z)),
        static_cast<double>(2 * (q.y * q.z + q.x * q.w)),
        static_cast<double>(0),
        static_cast<double>(2 * (q.x * q.z + q.y * q.w)),
        static_cast<double>(2 * (q.y * q.z - q.x * q.w)),
        static_cast<double>(1 - 2 * (q.x * q.x + q.y * q.y)),
        static_cast<double>(0), static_cast<double>(x.x / voxel_size),
        static_cast<double>(x.y / voxel_size),
        static_cast<double>(x.z / voxel_size), static_cast<double>(1));
    // mat.postScale(openvdb::math::Vec3d(1.0/voxel_size, 1.0/voxel_size, 1.0/voxel_size));

    matrices.push_back(mat);
    if (stream.peek() == std::ifstream::traits_type::eof()) break;
  }
  return matrices;
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
  filter.dilate(3);
  std::cout << "dilated" << std::endl;
  // for (int i = 0; i < 1; ++i) {
  //   filter.meanCurvature();
  // }
  std::cout << "curved" << std::endl;
  filter.erode();
  std::cout << "eroded" << std::endl;
  // filter.meanCurvature();
  // std::cout << "final smooth" << std::endl;
  // filter.gaussian();
  std::cout << "gaussian" << std::endl;
  //
  std::vector<openvdb::Vec3s> point_list;
  std::vector<openvdb::Vec3I> prim_list;
  load_obj("/home/kennychufk/workspace/pythonWs/alluvion-optim/buoy.obj",
           point_list, prim_list, 1 / voxelSize);
  // openvdb::FloatGrid::Ptr collision_obj_ls =
  // openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*linearTransform,
  // points, triangles);
  openvdb::tools::QuadAndTriangleDataAdapter<openvdb::Vec3s, openvdb::Vec3I>
      mesh(point_list, prim_list);

  int conversionFlags = 0;  // unsignedDistanceFieldConversion ?
                            // openvdb::tools::UNSIGNED_DISTANCE_FIELD : 0;
  openvdb::math::Mat4d mat = openvdb::math::Mat4d::identity();
  std::vector<openvdb::math::Mat4d> matrices = read_pile(
      "/home/kennychufk/workspace/pythonWs/alluvion-optim/rl-truth-f2caa5/"
      "73.pile",
      voxelSize);
  // openvdb::math::Transform::Ptr transform =
  //     openvdb::math::Transform::createLinearTransform(matrices[1]);
  openvdb::math::Transform::Ptr transform =
      openvdb::math::Transform::createLinearTransform(voxelSize);
  // transform->preScale(voxel_size);
  // transform->postTranslate(openvdb::math::Vec3d(voxel_size * 0.5, voxel_size
  // * 0.5, voxel_size * 0.5));

  float exBand = 3.0f;
  float inBand = 3.0f;
  openvdb::FloatGrid::Ptr collision_obj_ls =
      openvdb::tools::meshToVolume<openvdb::FloatGrid>(
          mesh, *transform, exBand, inBand, conversionFlags, nullptr);

  openvdb::FloatGrid::Ptr filtered_grid = filter.grid().deepCopy();
  for (int i = 1; i < 8; ++i) {
    openvdb::FloatGrid::Ptr transformed_collision_obj_ls =
        openvdb::FloatGrid::create();
    openvdb::tools::GridTransformer transformer(matrices[i]);
    transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(
        *collision_obj_ls, *transformed_collision_obj_ls);

    // openvdb::tools::csgDifference(*filtered_grid, *collision_obj_ls);
    openvdb::tools::csgDifference(*filtered_grid,
                                  *transformed_collision_obj_ls);
  }

  openvdb::tools::VolumeToMesh mesher;
  mesher.operator()<openvdb::FloatGrid>(*filtered_grid);
  // mesher.operator()<openvdb::FloatGrid>(*collision_obj_ls);
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
