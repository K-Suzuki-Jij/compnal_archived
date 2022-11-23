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
//  ssf_updater.hpp
//  compnal
//
//  Created by kohei on 2022/11/23.
//  
//

#ifndef COMPNAL_SOLVER_CMC_UTILITY_SSF_UPDATER_HPP_
#define COMPNAL_SOLVER_CMC_UTILITY_SSF_UPDATER_HPP_

#include <random>

namespace compnal {
namespace solver {
namespace cmc_utility {

template<class SystemType>
void SSFUpdater(SystemType *system,
                const std::int32_t num_sweeps,
                const typename SystemType::ValueType beta,
                const std::uint64_t seed,
                const std::function<bool(typename SystemType::ValueType, typename SystemType::ValueType, typename SystemType::ValueType)> ssf_transition_probability) {
   
   const std::int32_t system_size = system->GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<typename SystemType::ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
   
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (ssf_transition_probability(system->GetEnergyDifference(index), beta, dist_real(random_number_engine))) {
            system->Flip(index);
         }
      }
   }
   
}

} // namespace cmc_utility
} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_CMC_UTILITY_SSF_UPDATER_HPP_ */
