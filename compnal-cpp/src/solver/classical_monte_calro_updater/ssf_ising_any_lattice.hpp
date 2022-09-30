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
//  ssf_ising_any_lattice.hpp
//  compnal
//
//  Created by kohei on 2022/09/30.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_ISING_ANY_LATTICE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_ISING_ANY_LATTICE_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::Ising<lattice::AnyLattice, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::Ising<lattice::AnyLattice, RealType>::OPType> &sample,
                         const model::Ising<lattice::AnyLattice, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   const auto &adjacency_list = model.GetAdjacencyList();
   const auto system_size = model.GetSystemSize();
   
   for (std::int32_t i = 0; i < system_size; ++i) {
      const auto spin = sample[i];
      for (const auto &it: adjacency_list[i]) {
         (*energy_difference)[it.first] += -2*it.second*sample[it.first]*spin;
      }
   }
}


template<typename RealType>
void UpdateConfiguration(std::vector<typename model::Ising<lattice::AnyLattice, RealType>::OPType> *sample,
                         std::vector<RealType> *energy_difference,
                         const std::int32_t index,
                         const std::vector<std::vector<std::pair<std::int32_t, RealType>>> &adjacency_list
                         ) {
   const auto spin = (*sample)[index];
   for (const auto &it: adjacency_list[index]) {
      (*energy_difference)[it.first] += 4*it.second*(*sample)[it.first]*spin;
   }
   
   (*sample)[index] *= -1;
   (*energy_difference)[index] *= -1;
}

} // namespace updater
} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_UPDATER_SSF_ISING_ANY_LATTICE_HPP_ */
