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

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<typename RealType>
void pybind11SparseMatrixCRS(py::module &m) {
   
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   py::class_<CRS>(m, "CRS", py::module_local())
   .def(py::init<const std::size_t, const std::size_t>(), "row_dim"_a = 0, "col_dim"_a = 0)
   .def(py::init<const std::vector<std::vector<RealType>>&>(), "mat_vec")
   .def_readonly("row_dim", &CRS::row_dim)
   .def_readonly("col_dim", &CRS::col_dim)
   .def_readonly("row", &CRS::row)
   .def_readonly("col", &CRS::col)
   .def_readonly("val", &CRS::val)
   .def("assign", &CRS::Assign)
   .def("multiply_by_scalar", &CRS::MultiplyByScalar, "scalar"_a)
   .def("free", &CRS::Free)
   .def("Clear", &CRS::Clear)
   .def("is_symmetric", &CRS::isSymmetric, "threshold"_a = 0.000000000000001/*pow(10,-15)*/)
   .def("__repr__", [](const CRS& self) {
      std::ostringstream out;
      for (std::size_t i = 0; i < self.row_dim; ++i) {
         for (std::size_t j = self.row.at(i); j < self.row.at(i+1); ++j) {
            out << "M" << "[";
            out << std::noshowpos << std::left << std::setw(3) << i << "][";
            out << std::left << std::setw(3) << self.col[j] << "]=";
            out << std::showpos << self.val[j] << std::endl;
         }
      }
      out << std::noshowpos;
      return out.str();
   });
   
}

template<typename RealType>
void pybind11SparseMatrixBraketVector(py::module &m) {
   
   using BKV = compnal::sparse_matrix::BraketVector<RealType>;
   py::class_<BKV>(m, "BraketVector", py::module_local())
   .def(py::init<const std::size_t>(), "dim"_a = 0)
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
         if (i +1 == self.val.size()) {
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



#endif /* COMPNAL_PYBIND11_SPARSE_MATRIX_HPP_ */
