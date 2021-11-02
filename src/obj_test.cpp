#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/math/Math.h>
#include <openvdb/openvdb.h>
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

std::vector<openvdb::math::Transform::Ptr> read_pile(const char *filename) {
  std::ifstream stream(filename, std::ios::binary);
  stream.exceptions(std::ios_base::failbit);
  std::cout<<"good: "<<stream.good()<<std::endl;
  std::cout<<"eof: "<<stream.eof()<<std::endl;
  std::cout<<"fail: "<<stream.fail()<<std::endl;
  std::cout<<"bad: "<<stream.bad()<<std::endl;
  F3 x;
  F3 v;
  F4 q;
  F3 omega;
  std::vector<openvdb::math::Transform::Ptr> transforms;
  stream.seekg (0, stream.end);
  int length = stream.tellg();
  std::cout<<"length = "<<length<<std::endl;
  stream.seekg (0, stream.beg);
  while (true) {
    stream.read(reinterpret_cast<char *>(&x), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&v), sizeof(F3));
    stream.read(reinterpret_cast<char *>(&q), sizeof(F4));
    stream.read(reinterpret_cast<char *>(&omega), sizeof(F3));
    std::cout << x.x << ", " << x.y << ", " << x.z << std::endl;
    // openvdb::math::Mat4d mat(
    //     static_cast<double>(1 - 2 * (q.y * q.y + q.z * q.z)),
    //     static_cast<double>(2 * (q.x * q.y - q.z * q.w)),
    //     static_cast<double>(2 * (q.x * q.z + q.y * q.w)),
    //     static_cast<double>(x.x),
    //     static_cast<double>(2 * (q.x * q.y + q.z * q.w)),
    //     static_cast<double>(1 - 2 * (q.x * q.x + q.z * q.z)),
    //     static_cast<double>(2 * (q.y * q.z - q.x * q.w)),
    //     static_cast<double>(x.y),
    //     static_cast<double>(2 * (q.x * q.z - q.y * q.w)),
    //     static_cast<double>(2 * (q.y * q.z + q.x * q.w)),
    //     static_cast<double>(1 - 2 * (q.x * q.x + q.y * q.y)),
    //     static_cast<double>(x.z), static_cast<double>(0),
    //     static_cast<double>(0), static_cast<double>(0),
    //     static_cast<double>(1));

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
        static_cast<double>(0), static_cast<double>(x.x),
        static_cast<double>(x.y), static_cast<double>(x.z),
        static_cast<double>(1));

    transforms.push_back(openvdb::math::Transform::createLinearTransform(mat));
    if (stream.peek() == std::ifstream::traits_type::eof()) break;
  }
  return transforms;
}

int main(int argc, char *argv[]) {
  std::vector<openvdb::math::Transform::Ptr> transforms = read_pile("/home/kennychufk/workspace/pythonWs/alluvion-optim/rl-truth-f2caa5/73.pile");
}
