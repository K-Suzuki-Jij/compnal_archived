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
void pybind11ModelXXZ1D(py::module &m) {
   
   using XXZ1D = compnal::model::XXZ_1D<RealType>;
      
   py::class_<XXZ1D>(m, "XXZ_1D", py::module_local())
      .def(py::init<>())
      .def(py::init<const int>(), "system_size"_a)
      .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
      .def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a)
      .def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a)
      .def("set_J_z" , py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_z" , py::overload_cast<const RealType>(&XXZ1D::SetJz), "J_z"_a)
      .def("set_J_xy", py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("set_J_xy", py::overload_cast<const RealType>(&XXZ1D::SetJxy), "J_xy"_a)
      .def("calculate_target_dim"  , py::overload_cast<>(&XXZ1D::CalculateTargetDim, py::const_))
      .def("set_h_z" ,&XXZ1D::SetHz, "h_z"_a)
      .def("set_D_z" ,&XXZ1D::SetDz, "D_z"_a)
      .def("get_dim_onsite"        , &XXZ1D::GetDimOnsite)
      .def("get_num_conserved_quantity", &XXZ1D::GetNumConservedQuantity)
      .def("get_J_z" , py::overload_cast<>(&XXZ1D::GetJz, py::const_))
      .def("get_J_xy", py::overload_cast<>(&XXZ1D::GetJxy, py::const_))
      .def("get_J_z" , py::overload_cast<const std::int64_t>(&XXZ1D::GetJz, py::const_) , "index"_a)
      .def("get_J_xy", py::overload_cast<const std::int64_t>(&XXZ1D::GetJxy, py::const_), "index"_a)
      .def("get_h_z", &XXZ1D::GetHz)
      .def("get_D_z", &XXZ1D::GetDz)
      .def_property("system_size", &XXZ1D::GetSystemSize, &XXZ1D::SetSystemSize)
      .def_property("spin", &XXZ1D::GetMagnitudeSpin, &XXZ1D::SetMagnitudeSpin)
      .def_property("total_sz", &XXZ1D::GetTotalSz, &XXZ1D::SetTotalSz)
      .def_property("boundary_condition", &XXZ1D::GetBoundaryCondition, &XXZ1D::SetBoundaryCondition)
      .def_property_readonly("Ham", &XXZ1D::GetOnsiteOperatorHam)
      .def_property_readonly("Sx" , &XXZ1D::GetOnsiteOperatorSx )
      .def_property_readonly("iSy", &XXZ1D::GetOnsiteOperatoriSy)
      .def_property_readonly("Sz" , &XXZ1D::GetOnsiteOperatorSz )
      .def_property_readonly("Sp" , &XXZ1D::GetOnsiteOperatorSp )
      .def_property_readonly("Sm" , &XXZ1D::GetOnsiteOperatorSm )
      .def_static("create_onsite_operator_Sx"  , &XXZ1D::CreateOnsiteOperatorSx  , "spin"_a)
      .def_static("create_onsite_operator_iSy" , &XXZ1D::CreateOnsiteOperatoriSy , "spin"_a)
      .def_static("create_onsite_operator_Sz"  , &XXZ1D::CreateOnsiteOperatorSz  , "spin"_a)
      .def_static("create_onsite_operator_Sp"  , &XXZ1D::CreateOnsiteOperatorSp  , "spin"_a)
      .def_static("create_onsite_operator_Sm"  , &XXZ1D::CreateOnsiteOperatorSm  , "spin"_a)
      .def_static("create_onsite_operator_Ham" , &XXZ1D::CreateOnsiteOperatorHam , "spin"_a, "h_z"_a = 0.0, "D_z"_a = 0.0)
      .def("__repr__", [](const XXZ1D &self) {
         std::string bc = "None";
         if (self.GetBoundaryCondition() == compnal::utility::BoundaryCondition::OBC) {
            bc = "OBC";
         }
         else if (self.GetBoundaryCondition() == compnal::utility::BoundaryCondition::PBC) {
            bc = "PBC";
         }
         else if (self.GetBoundaryCondition() == compnal::utility::BoundaryCondition::SSD) {
            bc = "SSD";
         }
         
         std::ostringstream out;
         out << "////////////////////////" << std::endl;
         out << "Heisenberg Model Infomation:" << std::endl;
         out << "boundary_condition     = " << bc                             << std::endl;
         out << "system_size            = " << self.GetSystemSize()           << std::endl;
         out << "spin                   = " << self.GetMagnitudeSpin()        << std::endl;
         out << "total_2sz              = " << self.GetTotalSz()              << std::endl;
         out << "dim_target             = " << self.CalculateTargetDim()      << std::endl;
         out << "dim_onsite             = " << self.GetDimOnsite()            << std::endl;
         out << "num_conserved_quantity = " << self.GetNumConservedQuantity() << std::endl;
         out << "////////////////////////" << std::endl;
         out << "Heisenberg Model Interaction" << std::endl;
         out << "Sz-Sz Interaction: J_z ="     << std::endl;
         for (std::size_t i = 0; i < self.GetJz().size(); ++i) {
            out << i + 1 << "-th neighber: " << self.GetJz().at(i) << std::endl;
         }
         out << "Sx-Sx, Sy-Sy Interactions: J_xy =" << std::endl;
         for (std::size_t i = 0; i < self.GetJxy().size(); ++i) {
            out << i + 1 << "-th neighber: " << self.GetJxy().at(i) << std::endl;
         }
         out << "External Magnetic Fields for the z-direction: h_z = " << self.GetHz() << std::endl;
         out << "Uniaxial Anisotropy for the z-direction: D_z = "      << self.GetDz();
         return out.str();
      });
   
}

template<typename RealType>
void pybind11ModelU1Spin1D(py::module &m) {
   
   using U1Spin1D = compnal::model::U1Spin_1D<RealType>;
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   py::class_<U1Spin1D>(m, "U1Spin_1D", py::module_local())
      .def(py::init<>())
      .def(py::init<const int>(), "system_size"_a)
      .def(py::init<const int, const double>(), "system_size"_a, "spin"_a)
      .def("calculate_target_dim"  , py::overload_cast<>(&U1Spin1D::CalculateTargetDim, py::const_))
      .def("get_dim_onsite"        , &U1Spin1D::GetDimOnsite)
      .def("get_num_conserved_quantity", &U1Spin1D::GetNumConservedQuantity)
      .def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&U1Spin1D::AddOnsitePotential), "val"_a, "m"_a, "site"_a)
      .def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&U1Spin1D::AddOnsitePotential), "m"_a, "site"_a)
      .def("add_interaction", py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&U1Spin1D::AddInteraction), "val"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a)
      .def("add_interaction", py::overload_cast<const CRS&, const int, const CRS&, const int>(&U1Spin1D::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a)
      .def_property("system_size", &U1Spin1D::GetSystemSize, &U1Spin1D::SetSystemSize)
      .def_property("spin", &U1Spin1D::GetMagnitudeSpin, &U1Spin1D::SetMagnitudeSpin)
      .def_property("total_sz", &U1Spin1D::GetTotalSz, &U1Spin1D::SetTotalSz)
      .def_property_readonly("Sx" , &U1Spin1D::GetOnsiteOperatorSx )
      .def_property_readonly("iSy", &U1Spin1D::GetOnsiteOperatoriSy)
      .def_property_readonly("Sz" , &U1Spin1D::GetOnsiteOperatorSz )
      .def_property_readonly("Sp" , &U1Spin1D::GetOnsiteOperatorSp )
      .def_property_readonly("Sm" , &U1Spin1D::GetOnsiteOperatorSm )
      .def_static("create_onsite_operator_Sx"  , &U1Spin1D::CreateOnsiteOperatorSx  , "spin"_a)
      .def_static("create_onsite_operator_iSy" , &U1Spin1D::CreateOnsiteOperatoriSy , "spin"_a)
      .def_static("create_onsite_operator_Sz"  , &U1Spin1D::CreateOnsiteOperatorSz  , "spin"_a)
      .def_static("create_onsite_operator_Sp"  , &U1Spin1D::CreateOnsiteOperatorSp  , "spin"_a)
      .def_static("create_onsite_operator_Sm"  , &U1Spin1D::CreateOnsiteOperatorSm  , "spin"_a)
      .def("__repr__", [](const U1Spin1D &self) {
         std::ostringstream out;
         out << "U(1) Spin Model Infomation:" << std::endl;
         out << "system_size            = " << self.GetSystemSize()           << std::endl;
         out << "spin                   = " << self.GetMagnitudeSpin()        << std::endl;
         out << "total_2sz              = " << self.GetTotalSz()              << std::endl;
         out << "dim_target             = " << self.CalculateTargetDim()      << std::endl;
         out << "dim_onsite             = " << self.GetDimOnsite()            << std::endl;
         out << "num_conserved_quantity = " << self.GetNumConservedQuantity() << std::endl;
         return out.str();
      });
   
}


#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
