#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Statistics.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <stdexcept>

#include "alluvol/inside_squared_accumulation.hpp"
#include "alluvol/io.hpp"
#include "alluvol/level_set.hpp"

using namespace alluvol;
namespace py = pybind11;

PYBIND11_MODULE(_alluvol, m) {
  m.doc() = "OpenVDB for SPH";
  using namespace pybind11;

  py::class_<openvdb::FloatGrid, openvdb::FloatGrid::Ptr>(m, "FloatGrid")
      .def("activeVoxelCount", &openvdb::FloatGrid::activeVoxelCount)
      .def("statistics",
           [](const openvdb::FloatGrid& grid) {
             return openvdb::tools::statistics(grid.cbeginValueOn());
           })
      .def("calculate_inside_squared_distance_sum",
           [](const openvdb::FloatGrid& grid) {
             InsideSquaredAccumulation acc;
             openvdb::tools::accumulate(grid.cbeginValueOn(), acc);
             return acc.sum;
           })
      .def("deepCopy", &openvdb::FloatGrid::deepCopy);

  py::class_<openvdb::tools::VolumeToMesh>(m, "VolumeToMesh")
      .def(py::init<double, double, bool>(), py::arg("isovalue") = 0,
           py::arg("adaptivity") = 0,
           py::arg("relaxDisorientedTriangles") = true)
      .def("__call__", [](openvdb::tools::VolumeToMesh& volume_to_mesh,
                          const openvdb::FloatGrid& grid) {
        volume_to_mesh.operator()<openvdb::FloatGrid>(grid);
      });

  py::class_<openvdb::math::Stats>(m, "Stats")
      .def_property_readonly("min", &openvdb::math::Stats::min)
      .def_property_readonly("max", &openvdb::math::Stats::max)
      .def_property_readonly("range", &openvdb::math::Stats::range)
      .def_property_readonly("avg", &openvdb::math::Stats::avg)
      .def_property_readonly("mean", &openvdb::math::Stats::mean)
      .def_property_readonly("var", &openvdb::math::Stats::var)
      .def_property_readonly("std", &openvdb::math::Stats::std);

  m.def(
      "create_liquid_level_set",
      [](py::array_t<F> particle_x_py, openvdb::Real radius,
         openvdb::Real voxel_size, openvdb::Real mask_min,
         openvdb::Real mask_max, openvdb::Real half_width, U num_dilation,
         U num_mean_curvature, U num_erosion) {
        if (particle_x_py.ndim() != 2) {
          throw std::runtime_error("Number of dimensions must be 2");
        }
        if (particle_x_py.shape(1) != 3) {
          throw std::runtime_error(
              "The size of the second dimension must be 3");
        }
        std::vector<F3> particle_x =
            std::vector<F3>(reinterpret_cast<F3 const*>(particle_x_py.data()),
                            reinterpret_cast<F3 const*>(particle_x_py.data() +
                                                        particle_x_py.size()));
        return create_liquid_level_set(particle_x, radius, voxel_size, mask_min,
                                       mask_max, half_width, num_dilation,
                                       num_mean_curvature, num_erosion);
      },
      py::arg("particle_x"), py::arg("radius"), py::arg("voxel_size"),
      py::arg("mask_min") = 0, py::arg("mask_max") = 0,
      py::arg("half_width") = 3, py::arg("num_dilation") = 3,
      py::arg("num_mean_curvature") = 1, py::arg("num_erosion") = 3);
  m.def("create_mesh_level_set", &create_mesh_level_set, py::arg("filename"),
        py::arg("voxel_size"), py::arg("scale") = 1, py::arg("ex_band") = 3,
        py::arg("in_band") = 3);
  m.def("transform_level_set",
        py::overload_cast<openvdb::FloatGrid::Ptr, F3 const&, F4 const&,
                          openvdb::Real>(&transform_level_set),
        py::arg("grid"), py::arg("x"), py::arg("q"), py::arg("voxel_size"));
  m.def("csgDifference", &openvdb::tools::csgDifference<openvdb::FloatGrid>,
        py::arg("a"), py::arg("b"), py::arg("prune") = true);
  m.def("csgUnion", &openvdb::tools::csgUnion<openvdb::FloatGrid>, py::arg("a"),
        py::arg("b"), py::arg("prune") = true);
  m.def("csgIntersection", &openvdb::tools::csgIntersection<openvdb::FloatGrid>,
        py::arg("a"), py::arg("b"), py::arg("prune") = true);

  m.def("csgDifferenceCopy",
        &openvdb::tools::csgDifferenceCopy<openvdb::FloatGrid>, py::arg("a"),
        py::arg("b"));
  m.def("csgUnionCopy", &openvdb::tools::csgUnionCopy<openvdb::FloatGrid>,
        py::arg("a"), py::arg("b"));
  m.def("csgIntersectionCopy",
        &openvdb::tools::csgIntersectionCopy<openvdb::FloatGrid>, py::arg("a"),
        py::arg("b"));

  py::module m_io = m.def_submodule("io", "I/O operations");
  m_io.def("read_alu", &io::read_alu<1, F3>, py::arg("filename"));
  m_io.def("write_obj", &io::write_obj, py::arg("mesher"), py::arg("filename"));
  m_io.def("read_pile", &io::read_pile, py::arg("filename"),
           py::arg("voxel_size"));
}
