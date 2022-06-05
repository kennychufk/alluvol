#ifndef ALLUVOL_SQUARED_STAT_HPP
#define ALLUVOL_SQUARED_STAT_HPP

#include <openvdb/tools/Statistics.h>

namespace alluvol {

struct InsideSquaredAccumulation {
  double sum;
  inline void operator()(const openvdb::FloatGrid::ValueOnCIter &iter) {
    double inside_dist =
        (*iter) < 0 ? static_cast<double>(*iter) : static_cast<double>(0);
    if (iter.isVoxelValue())
      sum += inside_dist * inside_dist;
    else
      sum += inside_dist * inside_dist * iter.getVoxelCount();
  }
  // Accumulate another functor's Stats object into this functor's.
  inline void join(InsideSquaredAccumulation &other) { sum += other.sum; }
};
}  // namespace alluvol

#endif
