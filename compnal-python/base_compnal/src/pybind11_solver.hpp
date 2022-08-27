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

#include "../../../compnal-cpp/src/solver/all.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>


namespace compnal {
namespace wrapper {

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

void pybind11SolverCMCUpdater(py::module &m) {
   
   py::enum_<solver::CMCUpdater>(m, "CMCUpdater")
      .value("METROPOLIS", solver::CMCUpdater::METROPOLIS)
      .value("HEAT_BATH" , solver::CMCUpdater::HEAT_BATH);
   
}


template<class ModelType>
void pybind11SolverClassicalMonteCarlo(py::module &m, const std::string &post_name = "") {
   
   using CMC = solver::ClassicalMonteCarlo<ModelType>;
   using IndexType = typename CMC::IndexType;
   std::string name = std::string("ClassicalMonteCarlo") + post_name;

   auto py_class = py::class_<CMC>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const ModelType&, const solver::CMCUpdater>(), "model"_a, "updater"_a=solver::CMCUpdater::METROPOLIS);
   
   //Public Member Functions
   py_class.def("set_num_sweeps", &CMC::SetNumSweeps  , "num_sweeps"_a );
   py_class.def("set_num_samples", &CMC::SetNumSamples, "num_samples"_a);
   py_class.def("set_num_threads", &CMC::SetNumThreads, "num_threads"_a);
   py_class.def("set_inverse_temperature", &CMC::SetInverseTemperature, "inverse_temperature"_a);
   py_class.def("set_temperature", &CMC::SetTemperature, "temperature"_a);
   py_class.def("set_cmc_updater", &CMC::SetCMCUpdater, "cmc_updater"_a);
   py_class.def("get_num_sweeps", &CMC::GetNumSweeps);
   py_class.def("get_num_samples", &CMC::GetNumSamples);
   py_class.def("get_num_threads", &CMC::GetNumThreads);
   py_class.def("get_samples", &CMC::GetSamples);
   py_class.def("get_temperature", &CMC::GetTemperature);
   py_class.def("get_inverse_temperature", &CMC::GetInverseTemperature);
   py_class.def("get_seed", &CMC::GetSeed);
   py_class.def("get_cmc_updater", &CMC::GetCMCUpdater);
   py_class.def("run", py::overload_cast<>(&CMC::Run));
   py_class.def("run", py::overload_cast<const std::uint64_t>(&CMC::Run), "seed"_a);
   py_class.def("calculate_sample_average", &CMC::CalculateSampleAverage);
   py_class.def("calculate_sample_moment", &CMC::CalculateSampleMoment, "degree"_a);
   py_class.def("calculate_correlation", &CMC::CalculateCorrelation, "index_1"_a, "index_2"_a);
   py_class.def("calculate_correlation_list", &CMC::CalculateCorrelationList, "origin"_a, "index_list"_a);
   
   m.def("make_classical_monte_carlo", [](const ModelType &model,
                                          const solver::CMCUpdater cmc_updater) {
      return make_classical_monte_carlo(model, cmc_updater);
   }, "model"_a, "cmc_updater"_a = solver::CMCUpdater::METROPOLIS);
   
}




} // namespace wrapper
} // namespace compnal


#endif /* COMPNAL_PYBIND11_SOLVER_HPP_ */
