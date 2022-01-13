//
//  pybind11_sparse_matrix.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/09.
//

#ifndef COMPNAL_PYBIND11_SPARSE_MATRIX_HPP_
#define COMPNAL_PYBIND11_SPARSE_MATRIX_HPP_

#include "../../src/sparse_matrix/all.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<typename RealType>
void pybind11SparseMatrixCRS(py::module &m) {
   
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   auto c = py::class_<CRS>(m, "CRS", py::module_local());
   c.def(py::init<const std::int64_t, const std::int64_t>(), "row_dim"_a = 0, "col_dim"_a = 0);
   c.def(py::init<const std::vector<std::vector<RealType>>&>(), "mat_vec");
   c.def_readonly("row_dim", &CRS::row_dim);
   c.def_readonly("col_dim", &CRS::col_dim);
   c.def_readonly("row", &CRS::row);
   c.def_readonly("col", &CRS::col);
   c.def_readonly("val", &CRS::val);
   c.def("assign", &CRS::Assign);
   c.def("multiply_by_scalar", py::overload_cast<const RealType>(&CRS::MultiplyByScalar), "scalar"_a);
   c.def("free", &CRS::Free);
   c.def("Clear", &CRS::Clear);
   c.def("is_symmetric", &CRS::isSymmetric, "threshold"_a = 0.000000000000001/*pow(10,-15)*/);
   c.def("print_matrix", &CRS::PrintMatrix, "display_name"_a = "matrix");
   c.def("print_info", &CRS::PrintInfo, "display_name"_a = "matrix");
   c.def("__mul__", [](const CRS &lhs, const CRS &rhs) {
      CRS mat;
      compnal::sparse_matrix::CalculateMatrixMatrixProduct(&mat, 1.0, lhs, 1.0, rhs);
      return mat;
   }, py::is_operator());
   c.def("__add__", [](const CRS &lhs, const CRS &rhs) {
      CRS mat;
      compnal::sparse_matrix::CalculateMatrixMatrixSum(&mat, 1.0, lhs, 1.0, rhs);
      return mat;
   }, py::is_operator());
   c.def("__sub__", [](const CRS &lhs, const CRS &rhs) {
      CRS mat;
      compnal::sparse_matrix::CalculateMatrixMatrixSum(&mat, 1.0, lhs, -1.0, rhs);
      return mat;
   }, py::is_operator());
   c.def("__iadd__", [](CRS &self, const CRS &rhs) {
      const CRS lhs = self;
      compnal::sparse_matrix::CalculateMatrixMatrixSum(&self, 1.0, lhs, 1.0, rhs);
      return self;
   }, py::is_operator());
   c.def("__isub__", [](CRS &self, const CRS &rhs) {
      const CRS lhs = self;
      compnal::sparse_matrix::CalculateMatrixMatrixSum(&self, 1.0, lhs, -1.0, rhs);
      return self;
   }, py::is_operator());
   c.def("__imul__", [](CRS &self, const CRS &rhs) {
      const CRS lhs = self;
      compnal::sparse_matrix::CalculateMatrixMatrixProduct(&self, 1.0, lhs, 1.0, rhs);
      return self;
   }, py::is_operator());
   c.def("__eq__", [](CRS &self, const CRS &rhs) {
      return self == rhs;
   }, py::is_operator());
   c.def("__ne__", [](CRS &self, const CRS &rhs) {
      return self != rhs;
   }, py::is_operator());
   
}

template<typename RealType>
void pybind11SparseMatrixBraketVector(py::module &m) {
   
   using BKV = compnal::sparse_matrix::BraketVector<RealType>;
   py::class_<BKV>(m, "BraketVector", py::module_local())
      .def(py::init<const std::int64_t>(), "dim"_a = 0)
      .def(py::init<const std::vector<RealType>&>(), "vector"_a)
      .def_readonly("val", &BKV::val)
      .def("fill", &BKV::Fill, "val"_a)
      .def("free", &BKV::Fill)
      .def("clear", &BKV::Clear)
      .def("assign", py::overload_cast<const BKV&>(&BKV::Assign), "vector"_a)
      .def("assign", py::overload_cast<const std::vector<RealType>&>(&BKV::Assign), "vector"_a)
      .def("normalize", &BKV::Normalize, "scalar"_a = 1.0)
      .def("multiply_by_scalar", &BKV::MultiplyByScalar, "scalar"_a)
      .def("L2_norm", &BKV::L2Norm)
      .def("__repr__", [](const BKV& self) {
         std::ostringstream out;
         out << "[";
         for (std::size_t i = 0; i < self.val.size(); ++i) {
            if (i + 1 == self.val.size()) {
               out << self.val.at(i);
            }
            else {
               out << self.val.at(i) << ", ";
            }
         }
         out << "]";
         return out.str();
      });
   
}

void pybind11SparseMatrixParameters(py::module &m) {
   
   using PML = compnal::sparse_matrix::ParametersLanczos;
   py::class_<PML>(m, "ParamsDiag", py::module_local())
      .def(py::init<>())
      .def_readwrite("min_step", &PML::min_step)
      .def_readwrite("max_step", &PML::max_step)
      .def_readwrite("threshold", &PML::acc)
      .def_readwrite("flag_use_initial_vec", &PML::flag_use_initial_vec)
      .def_readwrite("flag_store_vec", &PML::flag_store_vec)
      .def_readwrite("flag_output_info", &PML::flag_output_info)
      .def_readonly("flag_symmetric_crs", &PML::flag_symmetric_crs);
   
   using PMCG = compnal::sparse_matrix::ParametersCG;
   py::class_<PMCG>(m, "ParamsCG", py::module_local())
      .def(py::init<>())
      .def_readwrite("max_step", &PMCG::max_step)
      .def_readwrite("threshold", &PMCG::acc)
      .def_readwrite("flag_use_initial_vec", &PMCG::flag_use_initial_vec)
      .def_readwrite("flag_output_info", &PMCG::flag_output_info)
      .def_readwrite("flag_symmetric_crs", &PMCG::flag_symmetric_crs);
   
   using PMII = compnal::sparse_matrix::ParametersII;
   py::class_<PMII>(m, "ParamsII", py::module_local())
      .def(py::init<>())
      .def_readwrite("max_step", &PMII::max_step)
      .def_readwrite("threshold", &PMII::acc)
      .def_readwrite("diag_add", &PMII::diag_add)
      .def_readwrite("flag_output_info", &PMII::flag_output_info)
      .def_readwrite("cg", &PMII::cg);
   
   
   using PMA = compnal::sparse_matrix::ParametersAll;
   py::class_<PMA>(m, "ParamsAll", py::module_local())
      .def(py::init<>())
      .def_readwrite("lanczos", &PMA::lanczos)
      .def_readwrite("ii", &PMA::ii);
   
}



#endif /* COMPNAL_PYBIND11_SPARSE_MATRIX_HPP_ */
