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
//  ssf_polynomial_ising_chain.hpp
//  compnal
//
//  Created by kohei on 2022/08/15.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CHAIN_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CHAIN_HPP_

#include "../../lattice/chain.hpp"
#include "../../model/classical/polynomial_ising.hpp"
#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<std::pair<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType, RealType>> *sample_delta,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
}

template<typename RealType>
void UpdateConfiguration(std::vector<std::pair<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType, RealType>> *sample_delta,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
}

template<typename RealType>
void SetEnergyDifference(std::vector<RealType> *energy_difference,
                         const std::vector<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType> &sample,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::Chain, RealType>::OPType;
   const std::int32_t system_size = model.GetSystemSize();
   const std::vector<RealType> &interaction = model.GetInteraction();
   
   if (model.lattice.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<RealType>::epsilon()) {
            continue;
         }
         const RealType target_ineraction = interaction[degree];
         for (std::int32_t index = 0; index < system_size; ++index) {
            RealType val = 0;
            for (std::int32_t i = 0; i < degree; ++i) {
               OPType sign = 1;
               for (std::int32_t j = 0; j < degree; ++j) {
                  std::int32_t connected_index = index - degree + 1 + i + j;
                  if (connected_index < 0) {
                     connected_index += system_size;
                  }
                  else if (connected_index >= system_size) {
                     connected_index -= system_size;
                  }
                  sign *= sample[connected_index];
               }
               val += sign*target_ineraction;
            }
            (*energy_difference)[index] = -2.0*val;
         }
      }
   }
   else if (model.lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<RealType>::epsilon()) {
            continue;
         }
         const RealType target_ineraction = interaction[degree];
         for (std::int32_t index = 0; index < system_size; ++index) {
            RealType val = 0;
            for (std::int32_t i = 0; i < degree; ++i) {
               if (index - degree + 1 + i < 0 || index + i >= system_size) {
                  continue;
               }
               OPType sign = 1;
               for (std::int32_t j = 0; j < degree; ++j) {
                  std::int32_t connected_index = index - degree + 1 + i + j;
                  sign *= sample[connected_index];
               }
               val += sign*target_ineraction;
            }
            (*energy_difference)[index] = -2.0*val;
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }
   
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType> *sample,
                         std::vector<RealType> *energy_difference,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
   using OPType = typename model::PolynomialIsing<lattice::Chain, RealType>::OPType;
   const std::vector<RealType> &interaction = model.GetInteraction();
   const std::int32_t system_size = model.GetSystemSize();
   
   if (model.lattice.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<RealType>::epsilon()) {
            continue;
         }
         const RealType target_ineraction = interaction[degree];
         for (std::int32_t i = 0; i < degree; ++i) {
            OPType sign = 1;
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = index - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += system_size;
               }
               else if (connected_index >= system_size) {
                  connected_index -= system_size;
               }
               sign *= (*sample)[connected_index];
            }
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = index - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += system_size;
               }
               else if (connected_index >= system_size) {
                  connected_index -= system_size;
               }
               if (connected_index != index) {
                  (*energy_difference)[connected_index] += 4*target_ineraction*sign;
               }
            }
         }
      }
   }
   else if (model.lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<RealType>::epsilon()) {
            continue;
         }
         const RealType target_ineraction = interaction[degree];
         for (std::int32_t i = 0; i < degree; ++i) {
            if (index - degree + 1 + i < 0 || index + i >= system_size) {
               continue;
            }
            OPType sign = 1;
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = index - degree + 1 + i + j;
               sign *= (*sample)[connected_index];
            }
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = index - degree + 1 + i + j;
               if (connected_index != index) {
                  (*energy_difference)[connected_index] += 4*target_ineraction*sign;
               }
            }
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }
   (*energy_difference)[index] *= -1;
   (*sample)[index] *= -1;
   
}

} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CHAIN_HPP_ */
