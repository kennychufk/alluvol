#ifndef ALLUVOL_DATA_TYPE_HPP
#define ALLUVOL_DATA_TYPE_HPP

namespace alluvol {
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
}  // namespace alluvol

#endif
