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
//  ssf_polynomial_ising_square.hpp
//  compnal
//
//  Created by kohei on 2022/08/16.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_SQUARE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_SQUARE_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::PolynomialIsing<lattice::Square, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::PolynomialIsing<lattice::Square, RealType>::OPType> &sample,
                         const model::PolynomialIsing<lattice::Square, RealType> &model) {
   
   if (static_cast<std::int32_t>(energy_difference->size()) != model.GetSystemSize()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::Chain, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::Chain, RealType>::ValueType;
   const std::vector<ValueType> &interaction = model.GetInteraction();
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<ValueType>::epsilon()) {
            continue;
         }
         const ValueType target_ineraction = interaction[degree];
         
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
               ValueType val = 0;
               for (std::int32_t i = 0; i < degree; ++i) {
                  OPType sign = 1;
                  for (std::int32_t j = 0; j < degree; ++j) {
                     std::int32_t connected_index = coo_x - degree + 1 + i + j;
                     if (connected_index < 0) {
                        connected_index += x_size;
                     }
                     else if (connected_index >= x_size) {
                        connected_index -= x_size;
                     }
                     sign *= sample[coo_y*x_size + connected_index];
                  }
                  val += sign*target_ineraction;
               }
               (*energy_difference)[coo_y*x_size + coo_x] = -2.0*val;
            }
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
               ValueType val = 0;
               for (std::int32_t i = 0; i < degree; ++i) {
                  OPType sign = 1;
                  for (std::int32_t j = 0; j < degree; ++j) {
                     std::int32_t connected_index = coo_y - degree + 1 + i + j;
                     if (connected_index < 0) {
                        connected_index += y_size;
                     }
                     else if (connected_index >= y_size) {
                        connected_index -= y_size;
                     }
                     sign *= sample[connected_index*x_size + coo_x];
                  }
                  val += sign*target_ineraction;
               }
               (*energy_difference)[coo_y*x_size + coo_x] = -2.0*val;
            }
         }
      }
   }
   
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<ValueType>::epsilon()) {
            continue;
         }
         const ValueType target_ineraction = interaction[degree];
         
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
               ValueType val = 0;
               for (std::int32_t i = 0; i < degree; ++i) {
                  if (coo_x - degree + 1 + i < 0 || coo_x + i >= x_size) {
                     continue;
                  }
                  OPType sign = 1;
                  for (std::int32_t j = 0; j < degree; ++j) {
                     std::int32_t connected_index = coo_x - degree + 1 + i + j;
                     sign *= (sample)[coo_y*x_size + connected_index];
                  }
                  val += sign*target_ineraction;
               }
               (*energy_difference)[coo_y*x_size + coo_x] = -2.0*val;
            }
         }

         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
               ValueType val = 0;
               for (std::int32_t i = 0; i < degree; ++i) {
                  if (coo_y - degree + 1 + i < 0 || coo_y + i >= y_size) {
                     continue;
                  }
                  OPType sign = 1;
                  for (std::int32_t j = 0; j < degree; ++j) {
                     std::int32_t connected_index = coo_y - degree + 1 + i + j;
                     sign *= (sample)[connected_index*x_size + coo_x];
                  }
                  val += sign*target_ineraction;
               }
               (*energy_difference)[coo_y*x_size + coo_x] = -2.0*val;
            }
         }
      }
   }
   
   
   else {
      throw std::runtime_error("Unsupported BinaryCondition");
   }
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::PolynomialIsing<lattice::Square, RealType>::OPType> *sample,
                         std::vector<typename model::PolynomialIsing<lattice::Square, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::Square, RealType> &model) {
   
   using OPType = typename model::PolynomialIsing<lattice::Chain, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::Chain, RealType>::ValueType;
   const std::vector<ValueType> &interaction = model.GetInteraction();
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   const std::int32_t coo_x = index%x_size;
   const std::int32_t coo_y = index/x_size;
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<ValueType>::epsilon()) {
            continue;
         }
         const ValueType target_ineraction = interaction[degree];
         
         // x-direction
         for (std::int32_t i = 0; i < degree; ++i) {
            OPType sign = 1;
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = coo_x - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += x_size;
               }
               else if (connected_index >= x_size) {
                  connected_index -= x_size;
               }
               sign *= (*sample)[coo_y*x_size + connected_index];
            }
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = coo_x - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += x_size;
               }
               else if (connected_index >= x_size) {
                  connected_index -= x_size;
               }
               if (connected_index != coo_x){
                  (*energy_difference)[coo_y*x_size + connected_index] += 4*target_ineraction*sign;
               }
            }
         }
         
         // y-direction
         for (std::int32_t i = 0; i < degree; ++i) {
            OPType sign = 1;
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = coo_y - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += y_size;
               }
               else if (connected_index >= y_size) {
                  connected_index -= y_size;
               }
               sign *= (*sample)[connected_index*x_size + coo_x];
            }
            for (std::int32_t j = 0; j < degree; ++j) {
               std::int32_t connected_index = coo_y - degree + 1 + i + j;
               if (connected_index < 0) {
                  connected_index += y_size;
               }
               else if (connected_index >= y_size) {
                  connected_index -= y_size;
               }
               if (connected_index != coo_y){
                  (*energy_difference)[connected_index*x_size + coo_x] += 4*target_ineraction*sign;
               }
            }
         }
      }
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      for (std::int32_t degree = 1; degree < interaction.size(); ++degree) {
         if (std::abs(interaction[degree]) <= std::numeric_limits<ValueType>::epsilon()) {
            continue;
         }
         const ValueType target_ineraction = interaction[degree];
         
         // x-direction
         for (std::int32_t i = std::max(coo_x - degree + 1, 0); i <= coo_x; ++i) {
            if (i > x_size - degree) {
               break;
            }
            OPType sign = 1;
            for (std::int32_t j = i; j < i + degree; ++j) {
               sign *= (*sample)[coo_y*x_size + j];
            }
            for (std::int32_t j = i; j < coo_x; ++j) {
               (*energy_difference)[coo_y*x_size + j] += 4*target_ineraction*sign;
            }
            for (std::int32_t j = coo_x + 1; j < i + degree; ++j) {
               (*energy_difference)[coo_y*x_size + j] += 4*target_ineraction*sign;
            }
         }
         
         // y-direction
         for (std::int32_t i = std::max(coo_y - degree + 1, 0); i <= coo_y; ++i) {
            if (i > y_size - degree) {
               break;
            }
            OPType sign = 1;
            for (std::int32_t j = i; j < i + degree; ++j) {
               sign *= (*sample)[j*x_size + coo_x];
            }
            for (std::int32_t j = i; j < coo_y; ++j) {
               (*energy_difference)[j*x_size + coo_x] += 4*target_ineraction*sign;
            }
            for (std::int32_t j = coo_y + 1; j < i + degree; ++j) {
               (*energy_difference)[j*x_size + coo_x] += 4*target_ineraction*sign;
            }
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BoundaryCondition");
   }
   (*energy_difference)[index] *= -1;
   (*sample)[index] *= -1;
   
}

} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_SQUARE_HPP_ */
