//
//  pybind11_declare.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/21.
//

#ifndef pybind11_declare_hpp
#define pybind11_declare_hpp

#include "solver.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<typename RealType>
void pybindCRS(py::module &m) {
   
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   py::class_<CRS>(m, "CRS")
   .def(py::init<>())
   .def(py::init<const int, const int>(), "row_dim"_a, "col_dim"_a)
   .def(py::init<const std::vector<std::vector<RealType>> &>(), "vec"_a)
   .def("resize_row", &CRS::ResizeRow, "row_dim"_a)
   .def("resize_col_val", &CRS::ResizeColVal, "num_elements"_a)
   .def("set_row_dim", &CRS::SetRowDim, "row_dim"_a)
   .def("set_col_dim", &CRS::SetColDim, "col_dim"_a)
   .def("get_row_dim", &CRS::GetRowDim)
   .def("get_col_dim", &CRS::GetColDim)
   .def("get_num_elements", &CRS::GetNumElements)
   .def("get_size_row", &CRS::GetSizeRow)
   .def("get_size_col", &CRS::GetSizeCol)
   .def("get_size_val", &CRS::GetSizeVal)
   .def("push_row", &CRS::PushRow, "size"_a)
   .def("push_col", &CRS::PushCol, "col"_a)
   .def("push_val", &CRS::PushVal, "val"_a)
   .def("get_row", &CRS::GetRow)
   .def("get_col", &CRS::GetCol)
   .def("get_val", &CRS::GetVal)
   .def("clear", &CRS::Clear)
   .def("free", &CRS::Free)
   .def("create_copy", &CRS::CreateCopy)
   .def("sort_col", &CRS::SortCol)
   .def("print", &CRS::Print, "display_name"_a = "Matrix")
   .def("print_info", &CRS::PrintInfo)
   .def("check_symmetric", &CRS::CheckSymmetric)
   .def("row", [](CRS& self, const int64_t index, const int64_t size) { return self.Row(index) = size;})
   .def("col", [](CRS& self, const int64_t index, const int64_t  col) { return self.Col(index) = col;} )
   .def("val", [](CRS& self, const int64_t index, const RealType val) { return self.Val(index) = val;} );
   
   
}

template<class ModelClass1D>
void pybindSolverExactDiag1D(py::module &m) {
   
   using ED1D = compnal::solver::ExactDiag1D<ModelClass1D>;
   
   py::class_<ED1D>(m, "ExactDiag1D")
   .def(py::init<const ModelClass1D &>(), "model"_a)
   //.def_readwrite("model"      , &ED1D::model)
   .def_readwrite("diag_params", &ED1D::diag_params)
   //.def_readwrite("diag_info"  , &ED1D::diag_info)
   .def("set_basis"    , &ED1D::SetBasis)
   .def("diagonalize"  , &ED1D::Diagonalize)
   .def("get_gs_vactor", &ED1D::GetGSVector)
   .def("get_gs_value" , &ED1D::GetGSValue)
   .def("get_basis"    , &ED1D::GetBasis)
   .def("get_dim"      , &ED1D::GetDim);
   
   
};



#endif /* pybind11_declare_hpp */
