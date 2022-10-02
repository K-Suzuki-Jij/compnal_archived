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
//  ssf_ising_chain.hpp
//  compnal
//
//  Created by kohei on 2022/09/25.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_ISING_CHAIN_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_ISING_CHAIN_HPP_

#include "../../lattice/chain.hpp"
#include "../../model/classical/ising.hpp"
#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::Ising<lattice::Chain, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::Ising<lattice::Chain, RealType>::OPType> &sample,
                         const model::Ising<lattice::Chain, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   const std::int32_t system_size = model.GetSystemSize();
   const typename model::Ising<lattice::Chain, RealType>::LinearType linear = model.GetLinear();
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();

   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      for (std::int32_t index = 0; index < system_size - 1; ++index) {
         (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index + 1] - 2*linear*sample[index];
         (*energy_difference)[index + 1] += -2*quadratic*sample[index]*sample[index + 1];
      }
      (*energy_difference)[system_size - 1] += -2*quadratic*sample[system_size - 1]*sample[0] - 2*linear*sample[system_size - 1];
      (*energy_difference)[0] += -2*quadratic*sample[system_size - 1]*sample[0];
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      for (std::int32_t index = 0; index < system_size - 1; ++index) {
         (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index + 1] - 2*linear*sample[index];
         (*energy_difference)[index + 1] += -2*quadratic*sample[index]*sample[index + 1];
      }
      (*energy_difference)[system_size - 1] += -2*linear*sample[system_size - 1];
   }
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }
   
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::Ising<lattice::Chain, RealType>::OPType> *sample,
                         std::vector<typename model::Ising<lattice::Chain, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::Ising<lattice::Chain, RealType> &model) {
   
   const std::int32_t system_size = model.GetSystemSize();
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      if (0 < index && index < system_size - 1) {
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
         (*energy_difference)[index] *= -1;
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
      }
      else if (index == 0) {
         (*energy_difference)[0] *= -1;
         (*energy_difference)[1] += 4*quadratic*(*sample)[0]*(*sample)[1];
         (*energy_difference)[system_size - 1] += 4*quadratic*(*sample)[0]*(*sample)[system_size - 1];
      }
      else {
         (*energy_difference)[0] *= 4*quadratic*(*sample)[system_size - 1]*(*sample)[0];
         (*energy_difference)[system_size - 2] += 4*quadratic*(*sample)[system_size - 2]*(*sample)[system_size - 1];
         (*energy_difference)[system_size - 1] *= -1;
      }
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      if (0 < index && index < system_size - 1) {
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
         (*energy_difference)[index] *= -1;
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
      }
      else if (index == 0) {
         (*energy_difference)[0] *= -1;
         (*energy_difference)[1] += 4*quadratic*(*sample)[0]*(*sample)[1];
      }
      else {
         (*energy_difference)[system_size - 2] += 4*quadratic*(*sample)[system_size - 2]*(*sample)[system_size - 1];
         (*energy_difference)[system_size - 1] *= -1;
      }
   }
   else {
      throw std::runtime_error("Unsupported BoundaryCondition");
   }

   (*sample)[index] *= -1;
   
}


} // namespace updater
} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_UPDATER_SSF_ISING_CHAIN_HPP_ */
