//
//  pybind11_model.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/09.
//

#ifndef COMPNAL_PYBIND11_MODEL_HPP_
#define COMPNAL_PYBIND11_MODEL_HPP_

#include "../../src/model/all.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<typename RealType>
void pybindModelXXZ1D(py::module &m) {
   
   using XXZ1D = compnal::model::XXZ_1D<RealType>;
      
   py::class_<XXZ1D>(m, "XXZ_1D")
      .def(py::init<>())
      .def(py::init<const int>(), "system_size"_a)
      .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
      .def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a)
      .def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a)
      .def("set_system_size", &XXZ1D::SetSystemSize, "system_size")
      .def("set_spin"       , &XXZ1D::SetMagnitudeSpin, "spin")
      .def("set_total_sz"   , &XXZ1D::SetTotalSz, "total_sz")
      .def("set_boundary_condition", &XXZ1D::SetBoundaryCondition, "boundary_condition"_a)
      .def("set_J_z" , py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_z" , py::overload_cast<const RealType>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_xy", py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("set_J_xy", py::overload_cast<const RealType>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("set_h_z" ,&XXZ1D::SetHz, "h_z"_a)
      .def("set_D_z" ,&XXZ1D::SetDz, "D_z"_a)
      .def("calculate_target_dim"  , &XXZ1D::CalculateTargetDim)
      .def("get_boundary_condition", &XXZ1D::GetBoundaryCondition)
      .def("get_system_size"       , &XXZ1D::GetSystemSize)
      .def("get_dim_onsite"        , &XXZ1D::GetDimOnsite)
      .def("get_num_conserved_quantity", &XXZ1D::GetNumConservedQuantity)
      .def("get_spin"              , &XXZ1D::GetMagnitudeSpin)
      .def("get_total_sz"          , &XXZ1D::GetTotalSz)
      .def("get_onsite_operator_Ham", &XXZ1D::GetOnsiteOperatorHam)
      .def("get_onsite_operator_Sx" , &XXZ1D::GetOnsiteOperatorSx )
      .def("get_onsite_operator_iSy", &XXZ1D::GetOnsiteOperatoriSy)
      .def("get_onsite_operator_Sz" , &XXZ1D::GetOnsiteOperatorSz )
      .def("get_onsite_operator_Sp" , &XXZ1D::GetOnsiteOperatorSp )
      .def("get_onsite_operator_Sm" , &XXZ1D::GetOnsiteOperatorSm )
      .def("get_J_z" , py::overload_cast<>(&XXZ1D::GetJz, py::const_))
      .def("get_J_xy", py::overload_cast<>(&XXZ1D::GetJxy, py::const_))
      .def("get_J_z" , py::overload_cast<const std::size_t>(&XXZ1D::GetJz, py::const_) , "index"_a)
      .def("get_J_xy", py::overload_cast<const std::size_t>(&XXZ1D::GetJxy, py::const_), "index"_a)
      .def("get_h_z", &XXZ1D::GetHz)
      .def("get_D_z", &XXZ1D::GetDz)
      .def_static("create_onsite_operator_Sx"  , &XXZ1D::CreateOnsiteOperatorSx  , "spin"_a)
      .def_static("create_onsite_operator_iSy" , &XXZ1D::CreateOnsiteOperatoriSy , "spin"_a)
      .def_static("create_onsite_operator_Sz"  , &XXZ1D::CreateOnsiteOperatorSz  , "spin"_a)
      .def_static("create_onsite_operator_SzSz", &XXZ1D::CreateOnsiteOperatorSzSz, "spin"_a)
      .def_static("create_onsite_operator_Sp"  , &XXZ1D::CreateOnsiteOperatorSp  , "spin"_a)
      .def_static("create_onsite_operator_Sm"  , &XXZ1D::CreateOnsiteOperatorSm  , "spin"_a)
      .def_static("create_onsite_operator_Ham" , &XXZ1D::CreateOnsiteOperatorHam , "spin"_a, "h_z"_a = 0.0, "D_z"_a = 0.0)
   
   
   ;
   
   
}


#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
