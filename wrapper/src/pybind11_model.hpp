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
#include <pybind11/iostream.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

void pybind11ModelBoundaryCondition(py::module &m) {
   py::enum_<compnal::utility::BoundaryCondition>(m, "BoundaryCondition")
      .value("OBC", compnal::utility::BoundaryCondition::OBC)
      .value("PBC", compnal::utility::BoundaryCondition::PBC)
      .value("SSD", compnal::utility::BoundaryCondition::SSD)
   ;
}

template<typename RealType>
void pybind11ModelBaseElectron1D(py::module &m) {
   
   using BUE1D = compnal::model::BaseU1Electron_1D<RealType>;
   
   auto c = py::class_<BUE1D>(m, "_BaseU1Electron_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "total_electron"_a, "total_sz"_a);
   
   //Public Member Functions
   c.def("calculate_target_dim", py::overload_cast<>(&BUE1D::CalculateTargetDim, py::const_));
   c.def("calculate_target_dim", py::overload_cast<const int, const double>(&BUE1D::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   c.def("print_basis_onsite"  , [](BUE1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintBasisOnsite();
   });
   
   //Static Member Functions
   c.def_static("make_onsite_operator_c_up"  , &BUE1D::CreateOnsiteOperatorCUp);
   c.def_static("make_onsite_operator_c_down", &BUE1D::CreateOnsiteOperatorCDown);
   c.def_static("make_onsite_operator_c_up_dagger"  , &BUE1D::CreateOnsiteOperatorCUpDagger);
   c.def_static("make_onsite_operator_c_down_dagger", &BUE1D::CreateOnsiteOperatorCDownDagger);
   c.def_static("make_onsite_operator_nc_up", &BUE1D::CreateOnsiteOperatorNCUp);
   c.def_static("make_onsite_operator_nc_down", &BUE1D::CreateOnsiteOperatorNCDown);
   c.def_static("make_onsite_operator_nc", &BUE1D::CreateOnsiteOperatorNC);
   c.def_static("make_onsite_operator_sx", &BUE1D::CreateOnsiteOperatorSx);
   c.def_static("make_onsite_operator_isy", &BUE1D::CreateOnsiteOperatoriSy);
   c.def_static("make_onsite_operator_sz", &BUE1D::CreateOnsiteOperatorSz);
   c.def_static("make_onsite_operator_sp", &BUE1D::CreateOnsiteOperatorSp);
   c.def_static("make_onsite_operator_sm", &BUE1D::CreateOnsiteOperatorSm);
   
   //Properties
   c.def_property("system_size"   , &BUE1D::GetSystemSize   , &BUE1D::SetSystemSize);
   c.def_property("total_electron", &BUE1D::GetTotalElectron, &BUE1D::SetTotalElectron);
   c.def_property("total_sz"      , &BUE1D::GetTotalSz      , &BUE1D::SetTotalSz);
   
   //Read Only Properties
   c.def_property_readonly("c_up"         , &BUE1D::GetOnsiteOperatorCUp        );
   c.def_property_readonly("c_down"       , &BUE1D::GetOnsiteOperatorCDown      );
   c.def_property_readonly("c_up_dagger"  , &BUE1D::GetOnsiteOperatorCUpDagger  );
   c.def_property_readonly("c_down_dagger", &BUE1D::GetOnsiteOperatorCDownDagger);
   c.def_property_readonly("nc_up"        , &BUE1D::GetOnsiteOperatorNCUp       );
   c.def_property_readonly("nc_down"      , &BUE1D::GetOnsiteOperatorNCDown     );
   c.def_property_readonly("nc" , &BUE1D::GetOnsiteOperatorNC );
   c.def_property_readonly("sx" , &BUE1D::GetOnsiteOperatorSx );
   c.def_property_readonly("isy", &BUE1D::GetOnsiteOperatoriSy);
   c.def_property_readonly("sz" , &BUE1D::GetOnsiteOperatorSz );
   c.def_property_readonly("sp" , &BUE1D::GetOnsiteOperatorSp );
   c.def_property_readonly("sm" , &BUE1D::GetOnsiteOperatorSm );
   
}

template<typename RealType>
void pybind11ModelBaseSpin1D(py::module &m) {
   
   using BUS1D = compnal::model::BaseU1Spin_1D<RealType>;
   
   auto c = py::class_<BUS1D>(m, "_BaseU1Spin_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   
   //Public Member Functions
   c.def("calculate_target_dim", py::overload_cast<>(&BUS1D::CalculateTargetDim, py::const_));
   c.def("calculate_target_dim", py::overload_cast<const double>(&BUS1D::CalculateTargetDim, py::const_), "total_sz"_a);
   c.def("print_basis_onsite"  , [](BUS1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintBasisOnsite();
   });
   
   //Static Member Functions
   c.def_static("make_onsite_operator_sx" , &BUS1D::CreateOnsiteOperatorSx , "spin"_a);
   c.def_static("make_onsite_operator_isy", &BUS1D::CreateOnsiteOperatoriSy, "spin"_a);
   c.def_static("make_onsite_operator_sz" , &BUS1D::CreateOnsiteOperatorSz , "spin"_a);
   c.def_static("make_onsite_operator_sp" , &BUS1D::CreateOnsiteOperatorSp , "spin"_a);
   c.def_static("make_onsite_operator_sm" , &BUS1D::CreateOnsiteOperatorSm , "spin"_a);
   
   //Properties
   c.def_property("system_size" , &BUS1D::GetSystemSize, &BUS1D::SetSystemSize);
   c.def_property("spin"        , &BUS1D::GetMagnitudeSpin, &BUS1D::SetMagnitudeSpin);
   c.def_property("total_sz"    , &BUS1D::GetTotalSz, &BUS1D::SetTotalSz);
   
   //Read Only Properties
   c.def_property_readonly("sx" , &BUS1D::GetOnsiteOperatorSx );
   c.def_property_readonly("isy", &BUS1D::GetOnsiteOperatoriSy);
   c.def_property_readonly("sz" , &BUS1D::GetOnsiteOperatorSz );
   c.def_property_readonly("sp" , &BUS1D::GetOnsiteOperatorSp );
   c.def_property_readonly("sm" , &BUS1D::GetOnsiteOperatorSm );

}

template<typename RealType>
void pybind11ModelBaseSpinElectron1D(py::module &m) {
   
   using BUSE1D = compnal::model::BaseU1SpinElectron_1D<RealType>;
   
   auto c = py::class_<BUSE1D>(m, "_BaseU1SpinElectron_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   c.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   c.def(py::init<const int, const double, const int>(), "system_size"_a, "spin"_a, "total_electron"_a);
   
   //Public Member Functions
   c.def("calculate_target_dim", py::overload_cast<>(&BUSE1D::CalculateTargetDim, py::const_));
   c.def("calculate_target_dim", py::overload_cast<const int, const double>(&BUSE1D::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   c.def("print_basis_onsite"  , [](BUSE1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintBasisOnsite();
   });
   
   //Static Member Functions
   c.def_static("make_onsite_operator_c_up"         , &BUSE1D::CreateOnsiteOperatorCUp        , "spin"_a);
   c.def_static("make_onsite_operator_c_down"       , &BUSE1D::CreateOnsiteOperatorCDown      , "spin"_a);
   c.def_static("make_onsite_operator_c_up_dagger"  , &BUSE1D::CreateOnsiteOperatorCUpDagger  , "spin"_a);
   c.def_static("make_onsite_operator_c_down_dagger", &BUSE1D::CreateOnsiteOperatorCDownDagger, "spin"_a);
   c.def_static("make_onsite_operator_nc_up"        , &BUSE1D::CreateOnsiteOperatorNCUp       , "spin"_a);
   c.def_static("make_onsite_operator_nc_down"      , &BUSE1D::CreateOnsiteOperatorNCDown     , "spin"_a);
   c.def_static("make_onsite_operator_nc"           , &BUSE1D::CreateOnsiteOperatorNC         , "spin"_a);
   c.def_static("make_onsite_operator_sx_c"         , &BUSE1D::CreateOnsiteOperatorSxC        , "spin"_a);
   c.def_static("make_onsite_operator_isy_c"        , &BUSE1D::CreateOnsiteOperatoriSyC       , "spin"_a);
   c.def_static("make_onsite_operator_sz_c"         , &BUSE1D::CreateOnsiteOperatorSzC        , "spin"_a);
   c.def_static("make_onsite_operator_sp_c"         , &BUSE1D::CreateOnsiteOperatorSpC        , "spin"_a);
   c.def_static("make_onsite_operator_sm_c"         , &BUSE1D::CreateOnsiteOperatorSmC        , "spin"_a);
   c.def_static("make_onsite_operator_sx_l"         , &BUSE1D::CreateOnsiteOperatorSxL        , "spin"_a);
   c.def_static("make_onsite_operator_isy_l"        , &BUSE1D::CreateOnsiteOperatoriSyL       , "spin"_a);
   c.def_static("make_onsite_operator_sz_l"         , &BUSE1D::CreateOnsiteOperatorSzL        , "spin"_a);
   c.def_static("make_onsite_operator_sp_l"         , &BUSE1D::CreateOnsiteOperatorSpL        , "spin"_a);
   c.def_static("make_onsite_operator_sm_l"         , &BUSE1D::CreateOnsiteOperatorSmL        , "spin"_a);
   
   //Properties
   c.def_property("system_size"   , &BUSE1D::GetSystemSize, &BUSE1D::SetSystemSize);
   c.def_property("spin"          , &BUSE1D::GetMagnitudeLSpin, &BUSE1D::SetMagnitudeLSpin);
   c.def_property("total_sz"      , &BUSE1D::GetTotalSz, &BUSE1D::SetTotalSz);
   c.def_property("total_electron", &BUSE1D::GetTotalElectron, &BUSE1D::SetTotalElectron);
   
   //Read Only Properties
   c.def_property_readonly("c_up"         , &BUSE1D::GetOnsiteOperatorCUp);
   c.def_property_readonly("c_down"       , &BUSE1D::GetOnsiteOperatorCDown);
   c.def_property_readonly("c_up_dagger"  , &BUSE1D::GetOnsiteOperatorCUpDagger);
   c.def_property_readonly("c_down_dagger", &BUSE1D::GetOnsiteOperatorCDownDagger);
   c.def_property_readonly("nc_up"        , &BUSE1D::GetOnsiteOperatorNCUp);
   c.def_property_readonly("nc_down"      , &BUSE1D::GetOnsiteOperatorNCDown);
   c.def_property_readonly("nc"           , &BUSE1D::GetOnsiteOperatorNC);
   c.def_property_readonly("sx_c" , &BUSE1D::GetOnsiteOperatorSxC );
   c.def_property_readonly("isy_c", &BUSE1D::GetOnsiteOperatoriSyC);
   c.def_property_readonly("sz_c" , &BUSE1D::GetOnsiteOperatorSzC );
   c.def_property_readonly("sp_c" , &BUSE1D::GetOnsiteOperatorSpC );
   c.def_property_readonly("sm_c" , &BUSE1D::GetOnsiteOperatorSmC );
   c.def_property_readonly("sx_l" , &BUSE1D::GetOnsiteOperatorSxL );
   c.def_property_readonly("isy_l", &BUSE1D::GetOnsiteOperatoriSyL);
   c.def_property_readonly("sz_l" , &BUSE1D::GetOnsiteOperatorSzL );
   c.def_property_readonly("sp_l" , &BUSE1D::GetOnsiteOperatorSpL );
   c.def_property_readonly("sm_l" , &BUSE1D::GetOnsiteOperatorSmL );
   c.def_property_readonly("sc_sl", &BUSE1D::GetOnsiteOperatorSCSL);
   
}

template<typename RealType>
void pybind11ModelBaseSpinMultiElectrons1D(py::module &m) {
   
   using BUSME1D = compnal::model::BaseU1SpinMultiElectrons_1D<RealType>;
   
   auto c = py::class_<BUSME1D>(m, "_BaseU1SpinMultiElectrons_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   c.def(py::init<const int, const std::vector<int>>(), "system_size"_a, "total_electron_list"_a);
   c.def(py::init<const int, const double, const std::vector<int>>(), "system_size"_a, "spin"_a, "total_electron_list"_a);
   
   //Public Member Functions
   c.def("calculate_target_dim", py::overload_cast<>(&BUSME1D::CalculateTargetDim, py::const_));
   c.def("calculate_target_dim", py::overload_cast<const double>(&BUSME1D::CalculateTargetDim, py::const_), "total_sz"_a);
   c.def("calculate_target_dim", py::overload_cast<const std::vector<int>&, const double>(&BUSME1D::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   c.def("print_basis_onsite"  , [](BUSME1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintBasisOnsite();
   });
   
   //Static Member Functions
   c.def_static("make_onsite_operator_c_up", &BUSME1D::CreateOnsiteOperatorCUp, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_c_down", &BUSME1D::CreateOnsiteOperatorCDown, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_c_up_dagger", &BUSME1D::CreateOnsiteOperatorCUpDagger, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_c_down_dagger", &BUSME1D::CreateOnsiteOperatorCDownDagger, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_nc_up"  , &BUSME1D::CreateOnsiteOperatorNCUp, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_nc_down", &BUSME1D::CreateOnsiteOperatorNCDown, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_nc"     , py::overload_cast<const double, const int, const int>(&BUSME1D::CreateOnsiteOperatorNC), "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_nc"     , py::overload_cast<const double, const int>(&BUSME1D::CreateOnsiteOperatorNC), "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sx_c"   , &BUSME1D::CreateOnsiteOperatorSxC , "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_isy_c"  , &BUSME1D::CreateOnsiteOperatoriSyC, "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sz_c"   , &BUSME1D::CreateOnsiteOperatorSzC , "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sp_c"   , &BUSME1D::CreateOnsiteOperatorSpC , "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sm_c"   , &BUSME1D::CreateOnsiteOperatorSmC , "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sx_l"   , &BUSME1D::CreateOnsiteOperatorSxL , "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_isy_l"  , &BUSME1D::CreateOnsiteOperatoriSyL, "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sz_l"   , &BUSME1D::CreateOnsiteOperatorSzL , "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sp_l"   , &BUSME1D::CreateOnsiteOperatorSpL , "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_sm_l"   , &BUSME1D::CreateOnsiteOperatorSmL , "spin"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_scsl"   , py::overload_cast<const double, const int, const int>(&BUSME1D::CreateOnsiteOperatorSCSL), "spin"_a, "orbital"_a, "num_orbital"_a);
   c.def_static("make_onsite_operator_scsl"   , py::overload_cast<const double, const int>(&BUSME1D::CreateOnsiteOperatorSCSL), "spin"_a, "num_orbital"_a);
   
   //Properties
   c.def_property("system_size"   , &BUSME1D::GetSystemSize, &BUSME1D::SetSystemSize);
   c.def_property("spin"          , &BUSME1D::GetMagnitudeLSpin, &BUSME1D::SetMagnitudeLSpin);
   c.def_property("total_sz"      , &BUSME1D::GetTotalSz, &BUSME1D::SetTotalSz);
   c.def_property("total_electron", py::overload_cast<>(&BUSME1D::GetTotalElectron, py::const_), &BUSME1D::SetTotalElectron);
   
   //Read Only Properties
   c.def_property_readonly("dim_onsite_electron", &BUSME1D::GetDimOnsiteElectron);
   c.def_property_readonly("dim_onsite_all_electrons", &BUSME1D::GetDimOnsiteAllElectrons);
   c.def_property_readonly("dim_onsite_lspin", &BUSME1D::GetDimOnsiteLSpin);
   c.def_property_readonly("dim_onsite", &BUSME1D::GetDimOnsite);
   c.def_property_readonly("num_electron_orbital", &BUSME1D::GetNumElectronOrbital);
   
   c.def_property_readonly("c_up"         , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorCUp, py::const_));
   c.def_property_readonly("c_down"       , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorCDown, py::const_));
   c.def_property_readonly("c_up_dagger"  , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorCUpDagger, py::const_));
   c.def_property_readonly("c_down_dagger", py::overload_cast<>(&BUSME1D::GetOnsiteOperatorCDownDagger, py::const_));
   c.def_property_readonly("nc_up"        , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorNCUp, py::const_));
   c.def_property_readonly("nc_down"      , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorNCDown, py::const_));
   c.def_property_readonly("nc"           , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorNC, py::const_));
   c.def_property_readonly("nc_tot"       , &BUSME1D::GetOnsiteOperatorNCTot);
   c.def_property_readonly("sx_c" , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorSxC , py::const_));
   c.def_property_readonly("isy_c", py::overload_cast<>(&BUSME1D::GetOnsiteOperatoriSyC, py::const_));
   c.def_property_readonly("sz_c" , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorSzC , py::const_));
   c.def_property_readonly("sp_c" , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorSpC , py::const_));
   c.def_property_readonly("sm_c" , py::overload_cast<>(&BUSME1D::GetOnsiteOperatorSmC , py::const_));
   c.def_property_readonly("sx_l" , &BUSME1D::GetOnsiteOperatorSxL );
   c.def_property_readonly("isy_l", &BUSME1D::GetOnsiteOperatoriSyL);
   c.def_property_readonly("sz_l" , &BUSME1D::GetOnsiteOperatorSzL );
   c.def_property_readonly("sp_l" , &BUSME1D::GetOnsiteOperatorSpL );
   c.def_property_readonly("sm_l" , &BUSME1D::GetOnsiteOperatorSmL );
   c.def_property_readonly("sc_sl", py::overload_cast<>(&BUSME1D::GetOnsiteOperatorSCSL, py::const_));
   
}

template<typename RealType>
void pybind11ModelGeneral1D(py::module &m) {
   
   using BUS1D    = compnal::model::BaseU1Spin_1D<RealType>;
   using GM1Dspin = compnal::model::GeneralModel_1D<BUS1D>;
   
   using BUE1D    = compnal::model::BaseU1Electron_1D<RealType>;
   using GM1DElec = compnal::model::GeneralModel_1D<BUE1D>;
   
   using BUSE1D   = compnal::model::BaseU1SpinElectron_1D<RealType>;
   using GM1DSE   = compnal::model::GeneralModel_1D<BUSE1D>;
   
   using BUSME1D  = compnal::model::BaseU1SpinMultiElectrons_1D<RealType>;
   using GM1DSME  = compnal::model::GeneralModel_1D<BUSME1D>;
   
   using CRS = compnal::sparse_matrix::CRS<RealType>;
   
   py::class_<GM1Dspin, BUS1D> gm1d_spin(m, "U1Spin_1D"    , py::module_local());
   gm1d_spin.def(py::init<>());
   gm1d_spin.def(py::init<const int>(), "system_size"_a);
   gm1d_spin.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1Dspin::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_spin.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1Dspin::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_spin.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1Dspin::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_spin.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1Dspin::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   
   py::class_<GM1DElec, BUE1D> gm1d_elec(m, "U1Electron_1D", py::module_local());
   gm1d_elec.def(py::init<>());
   gm1d_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_elec.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   gm1d_elec.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1DElec::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_elec.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1DElec::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_elec.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1DElec::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_elec.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1DElec::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   
   py::class_<GM1DSE, BUSE1D> gm1d_spin_elec(m, "U1SpinElectron_1D", py::module_local());
   gm1d_spin_elec.def(py::init<>());
   gm1d_spin_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_spin_elec.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin_elec.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   gm1d_spin_elec.def(py::init<const int, const double, const int>(), "system_size"_a, "spin"_a, "total_electron"_a);
   gm1d_spin_elec.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1DSE::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_spin_elec.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1DSE::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_spin_elec.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1DSE::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_spin_elec.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1DSE::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   
   py::class_<GM1DSME, BUSME1D> gm1d_spin_multi_elec(m, "U1SpinMultiElectrons_1D", py::module_local());
   gm1d_spin_multi_elec.def(py::init<>());
   gm1d_spin_multi_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const std::vector<int>&>(), "system_size"_a, "total_electron"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const double, const std::vector<int>&>(), "system_size"_a, "spin"_a, "total_electron"_a);
   gm1d_spin_multi_elec.def("add_onsite_potential", py::overload_cast<const RealType, const CRS&, const int>(&GM1DSME::AddOnsitePotential), "value"_a, "m"_a, "site"_a);
   gm1d_spin_multi_elec.def("add_onsite_potential", py::overload_cast<const CRS&, const int>(&GM1DSME::AddOnsitePotential), "m"_a, "site"_a);
   gm1d_spin_multi_elec.def("add_interaction"     , py::overload_cast<const RealType, const CRS&, const int, const CRS&, const int>(&GM1DSME::AddInteraction), "value"_a, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   gm1d_spin_multi_elec.def("add_interaction"     , py::overload_cast<const CRS&, const int, const CRS&, const int>(&GM1DSME::AddInteraction), "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a);
   
}

template<typename RealType>
void pybind11ModelXXZ1D(py::module &m) {
   
   using XXZ1D = compnal::model::XXZ_1D<RealType>;
   
   auto c = py::class_<XXZ1D, compnal::model::BaseU1Spin_1D<RealType>>(m, "XXZ_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   c.def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a);
   c.def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a);
   
   //Public Member Functions
   c.def("set_J_z" , py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJz), "J_z"_a);
   c.def("set_J_z" , py::overload_cast<const RealType>(&XXZ1D::SetJz), "J_z"_a);
   c.def("set_J_xy", py::overload_cast<const std::vector<RealType>&>(&XXZ1D::SetJxy), "J_xy"_a);
   c.def("set_J_xy", py::overload_cast<const RealType>(&XXZ1D::SetJxy), "J_xy"_a);
   c.def("set_h_z" , &XXZ1D::SetHz, "h_z"_a);
   c.def("set_D_z" , &XXZ1D::SetDz, "D_z"_a);
   c.def("get_J_z" , py::overload_cast<>(&XXZ1D::GetJz, py::const_));
   c.def("get_J_xy", py::overload_cast<>(&XXZ1D::GetJxy, py::const_));
   c.def("get_J_z" , py::overload_cast<const std::int64_t>(&XXZ1D::GetJz, py::const_) , "index"_a);
   c.def("get_J_xy", py::overload_cast<const std::int64_t>(&XXZ1D::GetJxy, py::const_), "index"_a);
   c.def("get_h_z" , &XXZ1D::GetHz);
   c.def("get_D_z" , &XXZ1D::GetDz);
   c.def("print_info", [](XXZ1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintInfo();
   });
   
   //Properties
   c.def_property("boundary_condition", &XXZ1D::GetBoundaryCondition, &XXZ1D::SetBoundaryCondition);
   
   //Read Only Properties
   c.def_property_readonly("Ham", &XXZ1D::GetOnsiteOperatorHam);
   
}

template<typename RealType>
void pybind11ModelHubbard1D(py::module &m) {
   
   using HBM1D = compnal::model::Hubbard_1D<RealType>;
   
   auto c = py::class_<HBM1D, compnal::model::BaseU1Electron_1D<RealType>>(m, "Hubbard_1D", py::module_local());
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   c.def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a);
   c.def(py::init<const int, const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "total_electron"_a, "boundary_condition"_a);
   
   //Public Member Functions
   c.def("set_t"  , py::overload_cast<const RealType>(&HBM1D::SetHopping), "t"_a);
   c.def("set_t"  , py::overload_cast<const std::vector<RealType>&>(&HBM1D::SetHopping), "t"_a);
   c.def("set_V"  , py::overload_cast<const RealType>(&HBM1D::SetIntersiteCoulomb), "V"_a);
   c.def("set_V"  , py::overload_cast<const std::vector<RealType>&>(&HBM1D::SetIntersiteCoulomb), "V"_a);
   c.def("set_U"  , &HBM1D::SetOnsiteCoulomb, "U"_a);
   c.def("set_h_z", &HBM1D::SetMagneticField, "h_z"_a);
   c.def("print_info", [](HBM1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintInfo();
   });
   
   //Properties
   c.def_property("boundary_condition", &HBM1D::GetBoundaryCondition, &HBM1D::SetBoundaryCondition);
   
   //Read Only Properties
   c.def_property_readonly("Ham", &HBM1D::GetOnsiteOperatorHam);
   
}

template<typename RealType>
void pybind11ModelKondoLattice1D(py::module &m) {
   
   using KLM1D = compnal::model::KondoLattice_1D<RealType>;
   
   auto c = py::class_<KLM1D, compnal::model::BaseU1SpinElectron_1D<RealType>>(m, "KondoLattice1D");
   
   //Constructors
   c.def(py::init<>());
   c.def(py::init<const int>(), "system_size"_a);
   c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   c.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   c.def(py::init<const int, const double, const int>(), "system_size"_a, "spin"_a, "total_electron"_a);
   c.def(py::init<const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "boundary_condition"_a);
   c.def(py::init<const int, const double, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "boundary_condition"_a);
   c.def(py::init<const int, const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "total_electron"_a, "boundary_condition"_a);
   c.def(py::init<const int, const double, const int, const compnal::utility::BoundaryCondition>(), "system_size"_a, "spin"_a, "total_electron"_a, "boundary_condition"_a);
   
   //Member Functions
   c.def("set_hopping", py::overload_cast<const RealType>(&KLM1D::SetHopping), "t"_a);
   c.def("set_hopping", py::overload_cast<const std::vector<RealType>&>(&KLM1D::SetHopping), "t"_a);
   c.def("set_J"   , &KLM1D::SetJ  , "J"_a   );
   c.def("set_J_z" , &KLM1D::SetJz , "J_z"_a );
   c.def("set_J_xy", &KLM1D::SetJxy, "J_xy"_a);
   c.def("set_D_z" , &KLM1D::SetDz , "D_z"_a );
   c.def("set_h_z" , &KLM1D::SetHz , "D_z"_a );
   c.def("get_hopping", py::overload_cast<const std::int64_t>(&KLM1D::GetHopping, py::const_), "i"_a);
   c.def("get_hopping", py::overload_cast<>(&KLM1D::GetHopping, py::const_));
   c.def("get_J_z" , &KLM1D::GetJz );
   c.def("get_J_xy", &KLM1D::GetJxy);
   c.def("get_D_z" , &KLM1D::GetDz );
   c.def("get_h_z" , &KLM1D::GetHz );
   c.def("print_info", [](KLM1D &self) {
      py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
      self.PrintInfo();
   });

   //Properties
   c.def_property("boundary_condition", &KLM1D::GetBoundaryCondition, &KLM1D::SetBoundaryCondition);
   c.def_property_readonly("Ham", &KLM1D::GetOnsiteOperatorHam);
   
   
}

#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
