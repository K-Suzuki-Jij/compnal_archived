//
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
//  monte_carlo_updater.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_MONTE_CARLO_UPDATER_HPP_
#define COMPNAL_SOLVER_UPDATER_MONTE_CARLO_UPDATER_HPP_

#include "ssf_polynomial_ising_infinite_range.hpp"
#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType, class ModelType>
void ExecuteMetropolis(std::vector<utility::SpinType> *sample,
                       const ModelType &model,
                       const std::int32_t num_sweeps,
                       const RealType beta,
                       const std::uint64_t seed) {
   
   //using RealType = typename ModelType::ValueType;
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   const std::int32_t system_size = model.GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<RealType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
   
   // Prepare system size index list (0, 1, 2, ..., N - 1)
   std::vector<std::int32_t> system_size_index_list(system_size);
   std::iota(system_size_index_list.begin(), system_size_index_list.end(), 0);
   
   // Set energy difference
   std::vector<RealType> energy_difference(system_size);
   SetEnergyDifference<RealType>(&energy_difference, &system_size_index_list, *sample, model);
   
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (energy_difference[index] <= 0.0 || std::exp(-beta*energy_difference[index]) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, &system_size_index_list, index, model);
         }
      }
   }
}

} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_MODEL_MONTE_CARLO_UPDATER_HPP_ */
