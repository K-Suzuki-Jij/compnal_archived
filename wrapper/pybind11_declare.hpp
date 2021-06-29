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


template<class ModelClass1D>
inline void pybindSolverExactDiag1D(py::module &m) {
   
   using ED1D = compnal::solver::ExactDiag1D<ModelClass1D>;
   
   py::class_<ED1D>(m, "ExactDiag1D")
   .def(py::init<const ModelClass1D &>(), "model")
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
