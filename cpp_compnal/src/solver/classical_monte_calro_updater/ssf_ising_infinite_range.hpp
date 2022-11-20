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
//  ssf_ising_infinite_range.hpp
//  compnal
//
//  Created by kohei on 2022/09/29.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_ISING_INFINITE_RANGE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_ISING_INFINITE_RANGE_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::Ising<lattice::InfiniteRange, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::Ising<lattice::InfiniteRange, RealType>::OPType> &sample,
                         const model::Ising<lattice::InfiniteRange, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   const std::int32_t system_size = model.GetSystemSize();
   const typename model::Ising<lattice::Chain, RealType>::LinearType linear = model.GetLinear();
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   for (std::int32_t i = 0; i < system_size; ++i) {
      (*energy_difference)[i] += -2*linear*sample[i];
      for (std::int32_t j = i + 1; j < system_size; ++j) {
         (*energy_difference)[i] += -2*quadratic*sample[i]*sample[j];
         (*energy_difference)[j] += -2*quadratic*sample[i]*sample[j];
      }
   }
   
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::Ising<lattice::InfiniteRange, RealType>::OPType> *sample,
                         std::vector<typename model::Ising<lattice::InfiniteRange, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::Ising<lattice::InfiniteRange, RealType> &model) {
   
   const std::int32_t system_size = model.GetSystemSize();
   const auto spin = (*sample)[index];
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();
   for (std::int32_t i = 0; i < index; ++i) {
      (*energy_difference)[i] += 4*quadratic*(*sample)[i]*spin;
   }
   for (std::int32_t i = index + 1; i < system_size; ++i) {
      (*energy_difference)[i] += 4*quadratic*(*sample)[i]*spin;
   }
   (*energy_difference)[index] *= -1;
   (*sample)[index] *= -1;

}

} // namespace updater
} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_UPDATER_SSF_ISING_INFINITE_RANGE_HPP_ */
