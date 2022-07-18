//  Copyright 2022 Kohei Suzuki
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Created by Kohei Suzuki on 2022/07/13.
//

#ifndef COMPNAL_PYBIND11_SOLVER_HPP_
#define COMPNAL_PYBIND11_SOLVER_HPP_

#include "../../../cpp/src/solver/all.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>


namespace compnal {
namespace wrapper {

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

void pybind11SolverUpdater(py::module &m) {
   
   py::enum_<solver::Updater>(m, "Updater")
      .value("METROPOLIS", solver::Updater::METROPOLIS)
      .value("HEAT_BATH" , solver::Updater::HEAT_BATH);
   
}


template<class ModelType>
void pybind11SolverClassicalMonteCarlo(py::module &m) {
   
   using ClassicalMC = solver::ClassicalMonteCarlo<ModelType>;
   
   auto py_class = py::class_<ClassicalMC>(m, "ClassicalMonteCarlo", py::module_local());
   
   //Constructors
   py_class.def(py::init<const ModelType, const solver::Updater>(), "model"_a, "updater"_a);
   
   //Public Member Functions
   py_class.def("set_num_sweeps", &ClassicalMC::SetNumSweeps  , "num_sweeps"_a );
   py_class.def("set_num_samples", &ClassicalMC::SetNumSamples, "num_samples"_a);
   py_class.def("set_temperature", &ClassicalMC::SetTemperature, "T"_a);
   py_class.def("set_inverse_temperature", &ClassicalMC::SetInverseTemperature, "beta"_a);
   py_class.def("get_num_sweeps", &ClassicalMC::GetNumSweeps);
   py_class.def("get_num_samples", &ClassicalMC::GetNumSamples);
   py_class.def("get_samples", &ClassicalMC::GetSamples);
   py_class.def("get_energies", &ClassicalMC::GetEnergies);
   py_class.def("get_magnetizations", &ClassicalMC::GetMagnetizations);
   py_class.def("clear_samples", &ClassicalMC::ClearSamples);
   py_class.def("run", &ClassicalMC::Run);
   py_class.def_readonly("model", &ClassicalMC::model);
   py_class.def_readonly("updater", &ClassicalMC::updater);
   
}




} // namespace wrapper
} // namespace compnal


#endif /* COMPNAL_PYBIND11_SOLVER_HPP_ */
