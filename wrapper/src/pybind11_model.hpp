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
void pybind11ModelBaseSpin1D(py::module &m) {
   
   using CppC = compnal::model::BaseU1Spin_1D<RealType>;
   auto  py_c = py::class_<CppC>(m, "_BaseU1Spin_1D", py::module_local());
   
   //Constructors
   py_c.def(py::init<>());
   py_c.def(py::init<const int>(), "system_size"_a);
   py_c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   py_c.def(py::init<const int, const double, const double>(), "system_size"_a, "spin"_a, "total_sz"_a);

   //Public Member Functions
   py_c.def("_calculate_target_dim", py::overload_cast<>(&CppC::CalculateTargetDim, py::const_));
   py_c.def("_calculate_target_dim", py::overload_cast<const double>(&CppC::CalculateTargetDim, py::const_), "total_sz"_a);
   py_c.def("_set_system_size"   , &CppC::SetSystemSize   , "system_size"_a);
   py_c.def("_set_magnitude_spin", &CppC::SetMagnitudeSpin, "spin"_a       );
   py_c.def("_set_total_sz"      , &CppC::SetTotalSz      , "total_sz"_a   );
   py_c.def("_get_system_size"   , &CppC::GetSystemSize   );
   py_c.def("_get_dim_onsite"    , &CppC::GetDimOnsite    );
   py_c.def("_get_total_sz"      , &CppC::GetTotalSz      );
   py_c.def("_get_magnitude_spin", &CppC::GetMagnitudeSpin);
   py_c.def("_get_onsite_operator_sx" , &CppC::GetOnsiteOperatorSx );
   py_c.def("_get_onsite_operator_isy", &CppC::GetOnsiteOperatoriSy);
   py_c.def("_get_onsite_operator_sz" , &CppC::GetOnsiteOperatorSz );
   py_c.def("_get_onsite_operator_sp" , &CppC::GetOnsiteOperatorSp );
   py_c.def("_get_onsite_operator_sm" , &CppC::GetOnsiteOperatorSm );

   //Static Member Functions
   py_c.def_static("_make_onsite_operator_sx" , &CppC::CreateOnsiteOperatorSx , "spin"_a);
   py_c.def_static("_make_onsite_operator_isy", &CppC::CreateOnsiteOperatoriSy, "spin"_a);
   py_c.def_static("_make_onsite_operator_sz" , &CppC::CreateOnsiteOperatorSz , "spin"_a);
   py_c.def_static("_make_onsite_operator_sp" , &CppC::CreateOnsiteOperatorSp , "spin"_a);
   py_c.def_static("_make_onsite_operator_sm" , &CppC::CreateOnsiteOperatorSm , "spin"_a);
   py_c.def_static("_is_valid_q_number", py::overload_cast<const int, const double, const double>(&CppC::isValidQNumber), "system_size"_a, "spin"_a, "total_sz"_a);
   
}

template<typename RealType>
void pybind11ModelBaseElectron1D(py::module &m) {
   
   using CppC = compnal::model::BaseU1Electron_1D<RealType>;
   auto  py_c = py::class_<CppC>(m, "_BaseU1Electron_1D", py::module_local());
   
   //Constructors
   py_c.def(py::init<>());
   py_c.def(py::init<const int>(), "system_size"_a);
   py_c.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   py_c.def(py::init<const int, const int, const double>(), "system_size"_a, "total_electron"_a, "total_sz"_a);

   //Public Member Functions
   py_c.def("_calculate_target_dim", py::overload_cast<>(&CppC::CalculateTargetDim, py::const_));
   py_c.def("_calculate_target_dim", py::overload_cast<const int, const double>(&CppC::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   py_c.def("_set_system_size"   , &CppC::SetSystemSize   , "system_size"_a   );
   py_c.def("_set_total_sz"      , &CppC::SetTotalSz      , "total_sz"_a      );
   py_c.def("_set_total_electron", &CppC::SetTotalElectron, "total_electron"_a);
   py_c.def("_get_system_size"   , &CppC::GetSystemSize   );
   py_c.def("_get_total_sz"      , &CppC::GetTotalSz      );
   py_c.def("_get_total_electron", &CppC::GetTotalElectron);
   py_c.def("_get_dim_onsite"    , &CppC::GetDimOnsite    );
   py_c.def("_get_onsite_operator_c_up"         , &CppC::GetOnsiteOperatorCUp        );
   py_c.def("_get_onsite_operator_c_down"       , &CppC::GetOnsiteOperatorCDown      );
   py_c.def("_get_onsite_operator_c_up_dagger"  , &CppC::GetOnsiteOperatorCUpDagger  );
   py_c.def("_get_onsite_operator_c_down_dagger", &CppC::GetOnsiteOperatorCDownDagger);
   py_c.def("_get_onsite_operator_nc_up"        , &CppC::GetOnsiteOperatorNCUp       );
   py_c.def("_get_onsite_operator_nc_down"      , &CppC::GetOnsiteOperatorNCDown     );
   py_c.def("_get_onsite_operator_nc"           , &CppC::GetOnsiteOperatorNC         );
   py_c.def("_get_onsite_operator_sx"           , &CppC::GetOnsiteOperatorSx         );
   py_c.def("_get_onsite_operator_isy"          , &CppC::GetOnsiteOperatoriSy        );
   py_c.def("_get_onsite_operator_sz"           , &CppC::GetOnsiteOperatorSz         );
   py_c.def("_get_onsite_operator_sp"           , &CppC::GetOnsiteOperatorSp         );
   py_c.def("_get_onsite_operator_sm"           , &CppC::GetOnsiteOperatorSm         );
   
   //Static Member Functions
   py_c.def_static("_make_onsite_operator_c_up"         , &CppC::CreateOnsiteOperatorCUp        );
   py_c.def_static("_make_onsite_operator_c_down"       , &CppC::CreateOnsiteOperatorCDown      );
   py_c.def_static("_make_onsite_operator_c_up_dagger"  , &CppC::CreateOnsiteOperatorCUpDagger  );
   py_c.def_static("_make_onsite_operator_c_down_dagger", &CppC::CreateOnsiteOperatorCDownDagger);
   py_c.def_static("_make_onsite_operator_nc_up"        , &CppC::CreateOnsiteOperatorNCUp       );
   py_c.def_static("_make_onsite_operator_nc_down"      , &CppC::CreateOnsiteOperatorNCDown     );
   py_c.def_static("_make_onsite_operator_nc"           , &CppC::CreateOnsiteOperatorNC         );
   py_c.def_static("_make_onsite_operator_sx"           , &CppC::CreateOnsiteOperatorSx         );
   py_c.def_static("_make_onsite_operator_isy"          , &CppC::CreateOnsiteOperatoriSy        );
   py_c.def_static("_make_onsite_operator_sz"           , &CppC::CreateOnsiteOperatorSz         );
   py_c.def_static("_make_onsite_operator_sp"           , &CppC::CreateOnsiteOperatorSp         );
   py_c.def_static("_make_onsite_operator_sm"           , &CppC::CreateOnsiteOperatorSm         );
   py_c.def_static("_is_valid_q_number", py::overload_cast<const int, const int, const double>(&CppC::isValidQNumber), "system_size"_a, "total_electron"_a, "total_sz"_a);
   
}

template<typename RealType>
void pybind11ModelBaseSpinElectron1D(py::module &m) {
   
   using CppC = compnal::model::BaseU1SpinElectron_1D<RealType>;
   auto  py_c = py::class_<CppC>(m, "_BaseU1SpinElectron_1D", py::module_local());
   
   //Constructors
   py_c.def(py::init<>());
   py_c.def(py::init<const int>(), "system_size"_a);
   py_c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   py_c.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   py_c.def(py::init<const int, const double, const int>(), "system_size"_a, "spin"_a, "total_electron"_a);
   py_c.def(py::init<const int, const double, const int, const double>(), "system_size"_a, "spin"_a, "total_electron"_a, "total_sz"_a);

   //Public Member Functions
   py_c.def("_calculate_target_dim", py::overload_cast<>(&CppC::CalculateTargetDim, py::const_));
   py_c.def("_calculate_target_dim", py::overload_cast<const int, const double>(&CppC::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   py_c.def("_set_system_size"     , &CppC::SetSystemSize    , "system_size"_a   );
   py_c.def("_set_magnitude_spin"  , &CppC::SetMagnitudeLSpin, "spin"_a          );
   py_c.def("_set_total_sz"        , &CppC::SetTotalSz       , "total_sz"_a      );
   py_c.def("_set_total_electron"  , &CppC::SetTotalElectron , "total_electron"_a);
   py_c.def("_get_system_size"     , &CppC::GetSystemSize    );
   py_c.def("_get_magnitude_spin"  , &CppC::GetMagnitudeLSpin);
   py_c.def("_get_total_sz"        , &CppC::GetTotalSz       );
   py_c.def("_get_total_electron"  , &CppC::GetTotalElectron );
   py_c.def("_get_dim_onsite"      , &CppC::GetDimOnsite     );
   py_c.def("_get_onsite_operator_c_up"         , &CppC::GetOnsiteOperatorCUp        );
   py_c.def("_get_onsite_operator_c_down"       , &CppC::GetOnsiteOperatorCDown      );
   py_c.def("_get_onsite_operator_c_up_dagger"  , &CppC::GetOnsiteOperatorCUpDagger  );
   py_c.def("_get_onsite_operator_c_down_dagger", &CppC::GetOnsiteOperatorCDownDagger);
   py_c.def("_get_onsite_operator_nc_up"        , &CppC::GetOnsiteOperatorNCUp       );
   py_c.def("_get_onsite_operator_nc_down"      , &CppC::GetOnsiteOperatorNCDown     );
   py_c.def("_get_onsite_operator_nc"           , &CppC::GetOnsiteOperatorNC         );
   py_c.def("_get_onsite_operator_sxc" , &CppC::GetOnsiteOperatorSxC );
   py_c.def("_get_onsite_operator_isyc", &CppC::GetOnsiteOperatoriSyC);
   py_c.def("_get_onsite_operator_szc" , &CppC::GetOnsiteOperatorSzC );
   py_c.def("_get_onsite_operator_spc" , &CppC::GetOnsiteOperatorSpC );
   py_c.def("_get_onsite_operator_smc" , &CppC::GetOnsiteOperatorSmC );
   py_c.def("_get_onsite_operator_sxl" , &CppC::GetOnsiteOperatorSxL );
   py_c.def("_get_onsite_operator_isyl", &CppC::GetOnsiteOperatoriSyL);
   py_c.def("_get_onsite_operator_szl" , &CppC::GetOnsiteOperatorSzL );
   py_c.def("_get_onsite_operator_spl" , &CppC::GetOnsiteOperatorSpL );
   py_c.def("_get_onsite_operator_sml" , &CppC::GetOnsiteOperatorSmL );
   py_c.def("_get_onsite_operator_scsl", &CppC::GetOnsiteOperatorSCSL);
   
   //Static Member Functions
   py_c.def_static("_make_onsite_operator_c_up"         , &CppC::CreateOnsiteOperatorCUp        , "spin"_a);
   py_c.def_static("_make_onsite_operator_c_down"       , &CppC::CreateOnsiteOperatorCDown      , "spin"_a);
   py_c.def_static("_make_onsite_operator_c_up_dagger"  , &CppC::CreateOnsiteOperatorCUpDagger  , "spin"_a);
   py_c.def_static("_make_onsite_operator_c_down_dagger", &CppC::CreateOnsiteOperatorCDownDagger, "spin"_a);
   py_c.def_static("_make_onsite_operator_nc_up"        , &CppC::CreateOnsiteOperatorNCUp       , "spin"_a);
   py_c.def_static("_make_onsite_operator_nc_down"      , &CppC::CreateOnsiteOperatorNCDown     , "spin"_a);
   py_c.def_static("_make_onsite_operator_nc"           , &CppC::CreateOnsiteOperatorNC         , "spin"_a);
   py_c.def_static("_make_onsite_operator_sxc"          , &CppC::CreateOnsiteOperatorSxC        , "spin"_a);
   py_c.def_static("_make_onsite_operator_isyc"         , &CppC::CreateOnsiteOperatoriSyC       , "spin"_a);
   py_c.def_static("_make_onsite_operator_szc"          , &CppC::CreateOnsiteOperatorSzC        , "spin"_a);
   py_c.def_static("_make_onsite_operator_spc"          , &CppC::CreateOnsiteOperatorSpC        , "spin"_a);
   py_c.def_static("_make_onsite_operator_smc"          , &CppC::CreateOnsiteOperatorSmC        , "spin"_a);
   py_c.def_static("_make_onsite_operator_sxl"          , &CppC::CreateOnsiteOperatorSxL        , "spin"_a);
   py_c.def_static("_make_onsite_operator_isyl"         , &CppC::CreateOnsiteOperatoriSyL       , "spin"_a);
   py_c.def_static("_make_onsite_operator_szl"          , &CppC::CreateOnsiteOperatorSzL        , "spin"_a);
   py_c.def_static("_make_onsite_operator_spl"          , &CppC::CreateOnsiteOperatorSpL        , "spin"_a);
   py_c.def_static("_make_onsite_operator_sml"          , &CppC::CreateOnsiteOperatorSmL        , "spin"_a);

}

template<typename RealType>
void pybind11ModelBaseSpinMultiElectrons1D(py::module &m) {
   
   using CppC = compnal::model::BaseU1SpinMultiElectrons_1D<RealType>;
   auto py_c = py::class_<CppC>(m, "_BaseU1SpinMultiElectrons_1D", py::module_local());
   
   //Constructors
   py_c.def(py::init<>());
   py_c.def(py::init<const int>(), "system_size"_a);
   py_c.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   py_c.def(py::init<const int, const std::vector<int>>(), "system_size"_a, "total_electron"_a);
   py_c.def(py::init<const int, const double, const std::vector<int>>(), "system_size"_a, "spin"_a, "total_electron"_a);
   py_c.def(py::init<const int, const double, const std::vector<int>, const double>(), "system_size"_a, "spin"_a, "total_electron"_a, "total_sz"_a);

   //Public Member Functions
   py_c.def("_calculate_target_dim", py::overload_cast<>(&CppC::CalculateTargetDim, py::const_));
   py_c.def("_calculate_target_dim", py::overload_cast<const double>(&CppC::CalculateTargetDim, py::const_), "total_sz"_a);
   py_c.def("_calculate_target_dim", py::overload_cast<const std::vector<int>&, const double>(&CppC::CalculateTargetDim, py::const_), "total_electron"_a, "total_sz"_a);
   py_c.def("_set_system_size"         , &CppC::SetSystemSize    , "system_size"_a   );
   py_c.def("_set_magnitude_spin"      , &CppC::SetMagnitudeLSpin, "spin"_a          );
   py_c.def("_set_total_sz"            , &CppC::SetTotalSz       , "total_sz"_a      );
   py_c.def("_set_total_electron"      , &CppC::SetTotalElectron , "total_electron"_a);
   py_c.def("_get_system_size"         , &CppC::GetSystemSize        );
   py_c.def("_get_magnitude_spin"      , &CppC::GetMagnitudeLSpin    );
   py_c.def("_get_total_sz"            , &CppC::GetTotalSz           );
   py_c.def("_get_dim_onsite"          , &CppC::GetDimOnsite         );
   py_c.def("_get_num_electron_orbital", &CppC::GetNumElectronOrbital);
   py_c.def("_get_total_electron"      , py::overload_cast<>(&CppC::GetTotalElectron, py::const_));
   py_c.def("_get_total_electron"      , py::overload_cast<const int>(&CppC::GetTotalElectron, py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_c_up"         , py::overload_cast<>(&CppC::GetOnsiteOperatorCUp        , py::const_));
   py_c.def("_get_onsite_operator_c_down"       , py::overload_cast<>(&CppC::GetOnsiteOperatorCDown      , py::const_));
   py_c.def("_get_onsite_operator_c_up_dagger"  , py::overload_cast<>(&CppC::GetOnsiteOperatorCUpDagger  , py::const_));
   py_c.def("_get_onsite_operator_c_down_dagger", py::overload_cast<>(&CppC::GetOnsiteOperatorCDownDagger, py::const_));
   py_c.def("_get_onsite_operator_nc_up"        , py::overload_cast<>(&CppC::GetOnsiteOperatorNCUp       , py::const_));
   py_c.def("_get_onsite_operator_nc_down"      , py::overload_cast<>(&CppC::GetOnsiteOperatorNCDown     , py::const_));
   py_c.def("_get_onsite_operator_nc"           , py::overload_cast<>(&CppC::GetOnsiteOperatorNC         , py::const_));
   py_c.def("_get_onsite_operator_sx_c"         , py::overload_cast<>(&CppC::GetOnsiteOperatorSxC        , py::const_));
   py_c.def("_get_onsite_operator_isy_c"        , py::overload_cast<>(&CppC::GetOnsiteOperatoriSyC       , py::const_));
   py_c.def("_get_onsite_operator_sz_c"         , py::overload_cast<>(&CppC::GetOnsiteOperatorSzC        , py::const_));
   py_c.def("_get_onsite_operator_sp_c"         , py::overload_cast<>(&CppC::GetOnsiteOperatorSpC        , py::const_));
   py_c.def("_get_onsite_operator_sm_c"         , py::overload_cast<>(&CppC::GetOnsiteOperatorSmC        , py::const_));
   py_c.def("_get_onsite_operator_sc_sl"        , py::overload_cast<>(&CppC::GetOnsiteOperatorSCSL       , py::const_));
   py_c.def("_get_onsite_operator_c_up"         , py::overload_cast<const int>(&CppC::GetOnsiteOperatorCUp        , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_c_down"       , py::overload_cast<const int>(&CppC::GetOnsiteOperatorCDown      , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_c_up_dagger"  , py::overload_cast<const int>(&CppC::GetOnsiteOperatorCUpDagger  , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_c_down_dagger", py::overload_cast<const int>(&CppC::GetOnsiteOperatorCDownDagger, py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_nc_up"        , py::overload_cast<const int>(&CppC::GetOnsiteOperatorNCUp       , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_nc_down"      , py::overload_cast<const int>(&CppC::GetOnsiteOperatorNCDown     , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_nc"           , py::overload_cast<const int>(&CppC::GetOnsiteOperatorNC         , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_sx_c"         , py::overload_cast<const int>(&CppC::GetOnsiteOperatorSxC        , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_isy_c"        , py::overload_cast<const int>(&CppC::GetOnsiteOperatoriSyC       , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_sz_c"         , py::overload_cast<const int>(&CppC::GetOnsiteOperatorSzC        , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_sp_c"         , py::overload_cast<const int>(&CppC::GetOnsiteOperatorSpC        , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_sm_c"         , py::overload_cast<const int>(&CppC::GetOnsiteOperatorSmC        , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_sc_sl"        , py::overload_cast<const int>(&CppC::GetOnsiteOperatorSCSL       , py::const_), "orbital"_a);
   py_c.def("_get_onsite_operator_nc_tot", &CppC::GetOnsiteOperatorNCTot);
   py_c.def("_get_onsite_operator_sx_l"  , &CppC::GetOnsiteOperatorSxL );
   py_c.def("_get_onsite_operator_isy_l" , &CppC::GetOnsiteOperatoriSyL);
   py_c.def("_get_onsite_operator_sz_l"  , &CppC::GetOnsiteOperatorSzL );
   py_c.def("_get_onsite_operator_sp_l"  , &CppC::GetOnsiteOperatorSpL );
   py_c.def("_get_onsite_operator_sm_l"  , &CppC::GetOnsiteOperatorSmL );
   
   //Static Member Functions
   py_c.def_static("_make_onsite_operator_c_up"         , &CppC::CreateOnsiteOperatorCUp        , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_c_down"       , &CppC::CreateOnsiteOperatorCDown      , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_c_up_dagger"  , &CppC::CreateOnsiteOperatorCUpDagger  , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_c_down_dagger", &CppC::CreateOnsiteOperatorCDownDagger, "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_nc_up"        , &CppC::CreateOnsiteOperatorNCUp       , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_nc_down"      , &CppC::CreateOnsiteOperatorNCDown     , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_nc"           , &CppC::CreateOnsiteOperatorNC         , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sx_c"         , &CppC::CreateOnsiteOperatorSxC        , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_isy_c"        , &CppC::CreateOnsiteOperatoriSyC       , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sz_c"         , &CppC::CreateOnsiteOperatorSzC        , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sp_c"         , &CppC::CreateOnsiteOperatorSpC        , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sm_c"         , &CppC::CreateOnsiteOperatorSmC        , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_scsl"         , &CppC::CreateOnsiteOperatorSCSL       , "spin"_a, "orbital"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_nc_tot"       , &CppC::CreateOnsiteOperatorNCTot      , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sx_l"         , &CppC::CreateOnsiteOperatorSxL        , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_isy_l"        , &CppC::CreateOnsiteOperatoriSyL       , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sz_l"         , &CppC::CreateOnsiteOperatorSzL        , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sp_l"         , &CppC::CreateOnsiteOperatorSpL        , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_sm_l"         , &CppC::CreateOnsiteOperatorSmL        , "spin"_a, "num_orbital"_a);
   py_c.def_static("_make_onsite_operator_scsl_tot"     , &CppC::CreateOnsiteOperatorSCSLTot    , "spin"_a, "num_orbital"_a);
   
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
      
   py::class_<GM1Dspin, BUS1D> gm1d_spin(m, "_U1Spin_1D", py::module_local());
   gm1d_spin.def(py::init<>());
   gm1d_spin.def(py::init<const int>(), "system_size"_a);
   gm1d_spin.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin.def(py::init<const int, const double, const double>(), "system_size"_a, "spin"_a, "total_sz"_a);
   gm1d_spin.def("_add_onsite_potential", &GM1Dspin::AddOnsitePotential, "m"_a, "site"_a, "value"_a = 1.0);
   gm1d_spin.def("_add_interaction"     , &GM1Dspin::AddInteraction, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a, "value"_a = 1.0);
   
   py::class_<GM1DElec, BUE1D> gm1d_elec(m, "_U1Electron_1D", py::module_local());
   gm1d_elec.def(py::init<>());
   gm1d_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_elec.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   gm1d_elec.def(py::init<const int, const int, const double>(), "system_size"_a, "total_electron"_a, "total_sz"_a);
   gm1d_elec.def("_add_onsite_potential", &GM1DElec::AddOnsitePotential, "m"_a, "site"_a, "value"_a = 1.0);
   gm1d_elec.def("_add_interaction"     , &GM1DElec::AddInteraction, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a, "value"_a = 1.0);
   
   py::class_<GM1DSE, BUSE1D> gm1d_spin_elec(m, "_U1SpinElectron_1D", py::module_local());
   gm1d_spin_elec.def(py::init<>());
   gm1d_spin_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_spin_elec.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin_elec.def(py::init<const int, const int>(), "system_size"_a, "total_electron"_a);
   gm1d_spin_elec.def(py::init<const int, const double, const int>(), "system_size"_a, "spin"_a, "total_electron"_a);
   gm1d_spin_elec.def(py::init<const int, const double, const int, const double>(), "system_size"_a, "spin"_a, "total_electron"_a, "total_sz"_a);
   gm1d_spin_elec.def("_add_onsite_potential", &GM1DSE::AddOnsitePotential, "m"_a, "site"_a, "value"_a = 1.0);
   gm1d_spin_elec.def("_add_interaction"     , &GM1DSE::AddInteraction, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a, "value"_a = 1.0);
   
   py::class_<GM1DSME, BUSME1D> gm1d_spin_multi_elec(m, "_U1SpinMultiElectrons_1D", py::module_local());
   gm1d_spin_multi_elec.def(py::init<>());
   gm1d_spin_multi_elec.def(py::init<const int>(), "system_size"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const double>(), "system_size"_a, "spin"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const std::vector<int>&>(), "system_size"_a, "total_electron"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const double, const std::vector<int>&>(), "system_size"_a, "spin"_a, "total_electron"_a);
   gm1d_spin_multi_elec.def(py::init<const int, const double, const std::vector<int>, const double>(), "system_size"_a, "spin"_a, "total_electron"_a, "total_sz"_a);
   gm1d_spin_multi_elec.def("_add_onsite_potential", &GM1DSME::AddOnsitePotential, "m"_a, "site"_a, "value"_a = 1.0);
   gm1d_spin_multi_elec.def("_add_interaction"     , &GM1DSME::AddInteraction, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a, "value"_a = 1.0);
   
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
