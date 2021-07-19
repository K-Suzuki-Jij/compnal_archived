//
//  pybind11_declare.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/21.
//

#ifndef pybind11_declare_hpp
#define pybind11_declare_hpp

#include "solver.hpp"
#include "model.hpp"
#include "sparse_matrix.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<typename RealType>
void pybindHeisenberg1D(py::module &m) {
   
   using HB1D = compnal::model::Heisenberg1D<RealType>;
   
   py::class_<HB1D>(m, "Heisenberg1D")
   .def(py::init<const int>(), "system_size"_a)
   .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
   .def(py::init<const int, const BoundaryCondition>(), "system_size"_a, "boundary_condition"_a)
   .def(py::init<const int, const double, const BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a)
   .def("get_boundary_condition"    , &HB1D::GetBoundaryCondition   )
   .def("get_system_size"           , &HB1D::GetSystemSize          )
   .def("get_dim_onsite"            , &HB1D::GetDimOnsite           )
   .def("get_spin"                  , &HB1D::GetMagnitudeSpin       )
   .def("get_total_sz"              , &HB1D::GetTotalSz             )
   .def("get_num_conserved_quantity", &HB1D::GetNumConservedQuantity)
   .def("get_flag_recalc_basis"     , &HB1D::GetFlagRecalcBasis     )
   .def("get_operator_ham"          , &HB1D::GetOperatorHam         )
   .def("get_operator_Sx"           , &HB1D::GetOperatorSx          )
   .def("get_operator_iSy"          , &HB1D::GetOperatoriSy         )
   .def("get_operator_Sz"           , &HB1D::GetOperatorSz          )
   .def("get_operator_Sp"           , &HB1D::GetOperatorSp          )
   .def("get_operator_Sm"           , &HB1D::GetOperatorSm          )
   .def("get_Jz"  , &HB1D::GetJz )
   .def("get_Jxy" , &HB1D::GetJxy)
   .def("get_Hz"  , &HB1D::GetHz )
   .def("get_Dz"  , &HB1D::GetDz )
   .def("set_spin"    , &HB1D::SetMagnitudeSpin, "spin"_a)
   .def("set_total_sz", &HB1D::SetTotalSz      , "total_sz"_a)
   .def("set_system_size", &HB1D::SetSystemSize, "system_size"_a)
   .def("set_boundary_condition", &HB1D::SetBoundaryCondition, "boundary_condition"_a)
   
   
   
}


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
   .def("print_info", &CRS::PrintInfo)
   .def("check_symmetric", &CRS::CheckSymmetric)
   .def("reset_row", [](CRS& self, const int64_t index, const int64_t size) { self.Row(index) = size;})
   .def("reset_col", [](CRS& self, const int64_t index, const int64_t  col) { self.Col(index) = col;} )
   .def("reset_val", [](CRS& self, const int64_t index, const RealType val) { self.Val(index) = val;} )
   .def("row", py::overload_cast<const int64_t>(&CRS::Row), "index"_a)
   .def("col", py::overload_cast<const int64_t>(&CRS::Col), "index"_a)
   .def("val", py::overload_cast<const int64_t>(&CRS::Val), "index"_a)
   .def("print", [](CRS& self, const std::string display_name) {
      for (int64_t i = 0; i < self.GetRowDim(); ++i) {
         for (int64_t j = self.Row(i); j < self.Row(i + 1); ++j) {
            py::print(display_name, "[", i, "][", self.Col(j), "]=", self.Val(j), "sep"_a="");
         }
      }
   }, "display_name"_a = "Matrix");
   
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
