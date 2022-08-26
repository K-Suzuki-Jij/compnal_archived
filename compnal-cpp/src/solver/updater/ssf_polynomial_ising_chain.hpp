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
void SetEnergyDifference(std::vector<std::pair<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType, RealType>> *sample_energy_difference_pair,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample_energy_difference_pair->size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
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
                  sign *= (*sample_energy_difference_pair)[connected_index].first;
               }
               val += sign*target_ineraction;
            }
            (*sample_energy_difference_pair)[index].second = -2.0*val;
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
                  sign *= (*sample_energy_difference_pair)[connected_index].first;
               }
               val += sign*target_ineraction;
            }
            (*sample_energy_difference_pair)[index].second = -2.0*val;
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }

   
}

template<typename RealType>
void UpdateConfiguration(std::vector<std::pair<typename model::PolynomialIsing<lattice::Chain, RealType>::OPType, RealType>> *sample_energy_difference_pair,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::Chain, RealType> &model) {
   
   using OPType = typename model::PolynomialIsing<lattice::Chain, RealType>::OPType;
   const std::vector<RealType> &interaction = model.GetInteraction();
   const std::int32_t degree = model.GetDegree();
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
               sign *= (*sample_energy_difference_pair)[connected_index].first;
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
                  (*sample_energy_difference_pair)[connected_index].second += 4*target_ineraction*sign;
               }
            }
         }
      }
   }
   else if (model.lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      if (degree == 3) {
         const RealType target_ineraction_deg2 = interaction[2];
         const RealType target_ineraction_deg3 = interaction[3];
         if (2 <= index && index <= system_size - 3) {
            const OPType sign3_1 =
            (*sample_energy_difference_pair)[index - 2].first*
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            const OPType sign3_2 =
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first*
            (*sample_energy_difference_pair)[index + 1].first;
            const OPType sign3_3 =
            (*sample_energy_difference_pair)[index - 0].first*
            (*sample_energy_difference_pair)[index + 1].first*
            (*sample_energy_difference_pair)[index + 2].first;
            
            const OPType sign2_1 =
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            const OPType sign2_2 =
            (*sample_energy_difference_pair)[index - 0].first*
            (*sample_energy_difference_pair)[index + 1].first;
            
            (*sample_energy_difference_pair)[index - 2].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg3*(sign3_1 + sign3_2);
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg2*sign2_1;
            (*sample_energy_difference_pair)[index + 1].second += 4*target_ineraction_deg3*(sign3_2 + sign3_3);
            (*sample_energy_difference_pair)[index + 1].second += 4*target_ineraction_deg2*sign2_2;
            (*sample_energy_difference_pair)[index + 2].second += 4*target_ineraction_deg3*sign3_3;
         }
         else if (index == 0) {
            const OPType sign3_1 =
            (*sample_energy_difference_pair)[0].first*
            (*sample_energy_difference_pair)[1].first*
            (*sample_energy_difference_pair)[2].first;
            const OPType sign2_1 =
            (*sample_energy_difference_pair)[0].first*
            (*sample_energy_difference_pair)[1].first;
            
            (*sample_energy_difference_pair)[1].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[1].second += 4*target_ineraction_deg2*sign2_1;
            (*sample_energy_difference_pair)[2].second += 4*target_ineraction_deg3*sign3_1;
         }
         else if (index == 1) {
            const OPType sign3_1 =
            (*sample_energy_difference_pair)[0].first*
            (*sample_energy_difference_pair)[1].first*
            (*sample_energy_difference_pair)[2].first;
            const OPType sign3_2 =
            (*sample_energy_difference_pair)[1].first*
            (*sample_energy_difference_pair)[2].first*
            (*sample_energy_difference_pair)[3].first;
            const OPType sign2_1 =
            (*sample_energy_difference_pair)[0].first*
            (*sample_energy_difference_pair)[1].first;
            const OPType sign2_2 =
            (*sample_energy_difference_pair)[1].first*
            (*sample_energy_difference_pair)[2].first;
            
            (*sample_energy_difference_pair)[0].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[0].second += 4*target_ineraction_deg2*sign2_1;
            (*sample_energy_difference_pair)[2].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[2].second += 4*target_ineraction_deg3*sign3_2;
            (*sample_energy_difference_pair)[2].second += 4*target_ineraction_deg2*sign2_2;
            (*sample_energy_difference_pair)[3].second += 4*target_ineraction_deg3*sign3_2;
         }
         else if (index == system_size - 2) {
            const OPType sign3_1 =
            (*sample_energy_difference_pair)[index - 2].first*
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            const OPType sign3_2 =
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first*
            (*sample_energy_difference_pair)[index + 1].first;
            const OPType sign2_1 =
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            const OPType sign2_2 =
            (*sample_energy_difference_pair)[index + 0].first*
            (*sample_energy_difference_pair)[index + 1].first;
            
            (*sample_energy_difference_pair)[index - 2].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg3*sign3_2;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg2*sign2_1;
            (*sample_energy_difference_pair)[index + 1].second += 4*target_ineraction_deg3*sign3_2;
            (*sample_energy_difference_pair)[index + 1].second += 4*target_ineraction_deg2*sign2_2;
         }
         else if (index == system_size - 1) {
            const OPType sign3_1 =
            (*sample_energy_difference_pair)[index - 2].first*
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            const OPType sign2_1 =
            (*sample_energy_difference_pair)[index - 1].first*
            (*sample_energy_difference_pair)[index - 0].first;
            
            (*sample_energy_difference_pair)[index - 2].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg3*sign3_1;
            (*sample_energy_difference_pair)[index - 1].second += 4*target_ineraction_deg2*sign2_1;
         }
         else {
            throw std::runtime_error("Invalid index.");
         }
      }
      else {
         for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
            if (std::abs(interaction[degree]) <= std::numeric_limits<RealType>::epsilon()) {
               continue;
            }
            const RealType target_ineraction = interaction[degree];
            
            for (std::int32_t i = std::max(index - degree + 1, 0); i <= index; ++i) {
               if (i > system_size - degree) {
                  break;
               }
               OPType sign = 1;
               for (std::int32_t j = i; j < i + degree; ++j) {
                  sign *= (*sample_energy_difference_pair)[j].first;
               }
               for (std::int32_t j = i; j < index; ++j) {
                  (*sample_energy_difference_pair)[j].second += 4*target_ineraction*sign;
               }
               for (std::int32_t j = index + 1; j < i + degree; ++j) {
                  (*sample_energy_difference_pair)[j].second += 4*target_ineraction*sign;
               }
            }
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }
   (*sample_energy_difference_pair)[index].first *= -1;
   (*sample_energy_difference_pair)[index].second *= -1;
   
}


} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CHAIN_HPP_ */
