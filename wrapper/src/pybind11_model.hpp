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

void pybind11ModelBoundaryCondition(py::module &m) {
   py::enum_<compnal::utility::BoundaryCondition>(m, "BC")
           .value("OBC", compnal::utility::BoundaryCondition::OBC)
           .value("PBC", compnal::utility::BoundaryCondition::PBC)
           .value("SSD", compnal::utility::BoundaryCondition::SSD)
   ;
}

template<typename RealType>
void pybind11ModelBaseElectron1D(py::module &m) {
   
   using BUE1D = compnal::model::BaseU1Electron_1D<RealType>;
   
   py::class_<BUE1D>(m, "BaseU1Electron_1D", py::module_local())
   .def(py::init<>())
   .def(py::init<const int>(), "system_size"_a)
   .def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a)
   .def("print_basis_onsite"     , &BUE1D::PrintBasisOnsite)
   .def("calculate_target_dim"   , py::overload_cast<>(&BUE1D::CalculateTargetDim, py::const_))
   .def("calculate_target_dim"   , py::overload_cast<const int, const double>(&BUE1D::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a)
   .def_property("system_size"   , &BUE1D::GetSystemSize   , &BUE1D::SetSystemSize)
   .def_property("total_electron", &BUE1D::GetTotalElectron, &BUE1D::SetTotalElectron)
   .def_property("total_sz"      , &BUE1D::GetTotalSz      , &BUE1D::SetTotalSz)
   .def_property_readonly("num_conserved_quantity", &BUE1D::GetNumConservedQuantity)
   .def_property_readonly("CUp"         , &BUE1D::GetOnsiteOperatorCUp        )
   .def_property_readonly("CDown"       , &BUE1D::GetOnsiteOperatorCDown      )
   .def_property_readonly("CUpDagger"   , &BUE1D::GetOnsiteOperatorCUpDagger  )
   .def_property_readonly("CDownDagger" , &BUE1D::GetOnsiteOperatorCDownDagger)
   .def_property_readonly("NCUp"        , &BUE1D::GetOnsiteOperatorNCUp       )
   .def_property_readonly("NCDown"      , &BUE1D::GetOnsiteOperatorNCDown     )
   .def_property_readonly("NC"          , &BUE1D::GetOnsiteOperatorNC         )
   .def_property_readonly("Sx"          , &BUE1D::GetOnsiteOperatorSx         )
   .def_property_readonly("iSy"         , &BUE1D::GetOnsiteOperatoriSy        )
   .def_property_readonly("Sz"          , &BUE1D::GetOnsiteOperatorSz         )
   .def_property_readonly("Sp"          , &BUE1D::GetOnsiteOperatorSp         )
   .def_property_readonly("Sm"          , &BUE1D::GetOnsiteOperatorSm         );
}

template<typename RealType>
void pybind11ModelBaseSpin1D(py::module &m) {
   
   using BUS1D = compnal::model::BaseU1Spin_1D<RealType>;
   
   py::class_<BUS1D>(m, "BaseU1Spin_1D", py::module_local())
      .def(py::init<>())
      .def(py::init<const int>(), "system_size"_a)
      .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
      .def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a)
      .def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a)
      .def("print_basis_onsite"     , &BUS1D::PrintBasisOnsite)
      .def("calculate_target_dim"   , py::overload_cast<>(&BUS1D::CalculateTargetDim, py::const_))
      .def("calculate_target_dim"   , py::overload_cast<const double>(&BUS1D::CalculateTargetDim, py::const_), "total_sz"_a)
      .def_property_readonly("Sx" , &BUS1D::GetOnsiteOperatorSx )
      .def_property_readonly("iSy", &BUS1D::GetOnsiteOperatoriSy)
      .def_property_readonly("Sz" , &BUS1D::GetOnsiteOperatorSz )
      .def_property_readonly("Sp" , &BUS1D::GetOnsiteOperatorSp )
      .def_property_readonly("Sm" , &BUS1D::GetOnsiteOperatorSm )
      .def_property("system_size"       , &BUS1D::GetSystemSize, &BUS1D::SetSystemSize)
      .def_property("spin"              , &BUS1D::GetMagnitudeSpin, &BUS1D::SetMagnitudeSpin)
      .def_property("total_sz"          , &BUS1D::GetTotalSz, &BUS1D::SetTotalSz)
      .def_property("boundary_condition", &BUS1D::GetBoundaryCondition, &BUS1D::SetBoundaryCondition);

}

template<typename RealType>
void pybind11ModelGeneral1D(py::module &m) {
   
   using BUS1D    = compnal::model::BaseU1Spin_1D<RealType>;
   using GM1Dspin = compnal::model::GeneralModel_1D<BUS1D>;
   
   using BUE1D    = compnal::model::BaseU1Electron_1D<RealType>;
   using GM1DElec = compnal::model::GeneralModel_1D<BUE1D>;
   
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   py::class_<GM1Dspin, BUS1D> gm1d_spin(m, "U1Spin_1D"    , py::module_local());
   py::class_<GM1DElec, BUE1D> gm1d_elec(m, "U1Electron_1D", py::module_local());
   
   gm1d_spin.def(py::init<>());
   gm1d_spin.def(py::init<const int>(), "system_size"_a);
   gm1d_spin.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin.def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a);
   gm1d_spin.def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a);
   gm1d_spin.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1Dspin::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_spin.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1Dspin::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_spin.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1Dspin::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_spin.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1Dspin::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
      
   gm1d_elec.def(py::init<>());
   gm1d_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_elec.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   gm1d_elec.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1DElec::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_elec.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1DElec::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_elec.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1DElec::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_elec.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1DElec::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   
}

template<typename RealType>
void pybind11ModelXXZ1D(py::module &m) {
   
   using XXZ1D = compnal::model::XXZ_1D<RealType>;
   
   py::class_<XXZ1D, compnal::model::BaseU1Spin_1D<RealType>>(m, "XXZ_1D", py::module_local())
      .def(py::init<>())
      .def(py::init<const int>(), "system_size"_a)
      .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
      .def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a)
      .def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a)
      .def("set_J_z" , py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_z" , py::overload_cast<const RealType>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_xy", py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("set_J_xy", py::overload_cast<const RealType>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("set_h_z" ,&XXZ1D::SetHz, "h_z"_a)
      .def("set_D_z" ,&XXZ1D::SetDz, "D_z"_a)
      .def("get_J_z" , py::overload_cast<>(&XXZ1D::GetJz, py::const_))
      .def("get_J_xy", py::overload_cast<>(&XXZ1D::GetJxy, py::const_))
      .def("get_J_z" , py::overload_cast<const std::int64_t>(&XXZ1D::GetJz, py::const_) , "index"_a)
      .def("get_J_xy", py::overload_cast<const std::int64_t>(&XXZ1D::GetJxy, py::const_), "index"_a)
      .def("get_h_z", &XXZ1D::GetHz)
      .def("get_D_z", &XXZ1D::GetDz)
      .def("print_info", &XXZ1D::PrintInfo)
      .def_property_readonly("Ham", &XXZ1D::GetOnsiteOperatorHam);
   
}


#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
