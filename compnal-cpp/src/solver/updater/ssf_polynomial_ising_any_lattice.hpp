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
//  ssf_polynomial_ising_any_lattice.hpp
//  compnal
//
//  Created by kohei on 2022/08/16.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_ANY_LATTICE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_ANY_LATTICE_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> &sample,
                         const std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> &sign_list,
                         const model::PolynomialIsing<lattice::AnyLattice, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType;
   const std::vector<std::vector<std::int32_t>> &key_list = model.GetKeyList();
   const std::vector<ValueType> &value_list = model.GetValueList();
   
   for (std::size_t i = 0; i < key_list.size(); ++i) {
      for (const auto &index: key_list[i]) {
         (*energy_difference)[index] += -2*value_list[i]*sign_list[i];
      }
   }
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> *sample,
                         std::vector<RealType> *energy_difference,
                         std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> *sign_list,
                         const std::int32_t index,
                         const std::vector<std::vector<std::size_t>> &adjacency_list,
                         const std::vector<std::vector<std::int32_t>> &key_list,
                         const std::vector<RealType> &value_list) {
   
   for (const auto &interaction_index: adjacency_list[index]) {
      const RealType value = 4*value_list[interaction_index]*(*sign_list)[interaction_index];
      (*sign_list)[interaction_index] *= -1;
      for (const auto &update_index: key_list[interaction_index]) {
         if (update_index == index) continue;
         (*energy_difference)[update_index] += value;
      }
   }
   
   (*energy_difference)[index] *= -1;
   (*sample)[index] *= -1;
}

} // namespace updater
} // namespace solver
} // namespace compnal



#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_ANY_LATTICE_HPP_ */
