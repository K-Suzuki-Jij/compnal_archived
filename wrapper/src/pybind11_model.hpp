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
   py::enum_<compnal::utility::BoundaryCondition>(m, "bc")
           .value("OBC", compnal::utility::BoundaryCondition::OBC)
           .value("PBC", compnal::utility::BoundaryCondition::PBC)
           .value("SSD", compnal::utility::BoundaryCondition::SSD)
   ;
}

template<typename RealType>
void pybind11ModelXXZ1D(py::module &m) {
   
   using XXZ1D = compnal::model::XXZ_1D<RealType>;
      
   py::class_<XXZ1D>(m, "XXZ_1D", py::module_local())
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
      .def("__repr__", [](const XXZ1D& self) {
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
      })
   
   ;
   
   
}


#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
