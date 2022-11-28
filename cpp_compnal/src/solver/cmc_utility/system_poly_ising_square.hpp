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
//  system_poly_ising_square.hpp
//  compnal
//
//  Created by kohei on 2022/11/23.
//  
//

#ifndef COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_POLY_ISING_SQUARE_HPP_
#define COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_POLY_ISING_SQUARE_HPP_

#include "../../lattice/all.hpp"
#include "../../model/all.hpp"
#include "base_system.hpp"

namespace compnal {
namespace solver {
namespace cmc_utility {

template<typename RealType>
class CMCSystem<model::PolynomialIsing<lattice::Square, RealType>> {
   
   using ModelType = model::PolynomialIsing<lattice::Square, RealType>;

public:
   using ValueType = typename ModelType::ValueType;
   
   CMCSystem(const ModelType &model, const uint64_t seed):
   system_size_(model.GetSystemSize()),
   x_size_(model.GetLattice().GetXSize()),
   y_size_(model.GetLattice().GetYSize()),
   bc_(model.GetBoundaryCondition()),
   interaction_(model.GetInteraction()) {
      sample_ = GenerateRandomSpin(seed);
      energy_difference_ = GenerateEnergyDifference(sample_);
   }
   
   void Flip(const std::int32_t index) {
      const std::int32_t coo_x = index%x_size_;
      const std::int32_t coo_y = index/x_size_;
      if (bc_ == lattice::BoundaryCondition::PBC) {
         for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
            if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
               continue;
            }
            const typename ModelType::ValueType target_ineraction = interaction_[degree];
            
            for (std::int32_t i = 0; i < degree; ++i) {
               typename ModelType::OPType sign_x = 1;
               typename ModelType::OPType sign_y = 1;
               for (std::int32_t j = 0; j < degree; ++j) {
                  // x-direction
                  std::int32_t connected_index_x = coo_x - degree + 1 + i + j;
                  if (connected_index_x < 0) {
                     connected_index_x += x_size_;
                  }
                  else if (connected_index_x >= x_size_) {
                     connected_index_x -= x_size_;
                  }
                  sign_x *= sample_[coo_y*x_size_ + connected_index_x];
                  
                  // y-direction
                  std::int32_t connected_index_y = coo_y - degree + 1 + i + j;
                  if (connected_index_y < 0) {
                     connected_index_y += y_size_;
                  }
                  else if (connected_index_y >= y_size_) {
                     connected_index_y -= y_size_;
                  }
                  sign_y *= sample_[connected_index_y*x_size_ + coo_x];
               }
               for (std::int32_t j = 0; j < degree; ++j) {
                  // x-direction
                  std::int32_t connected_index_x = coo_x - degree + 1 + i + j;
                  if (connected_index_x < 0) {
                     connected_index_x += x_size_;
                  }
                  else if (connected_index_x >= x_size_) {
                     connected_index_x -= x_size_;
                  }
                  if (connected_index_x != coo_x){
                     energy_difference_[coo_y*x_size_ + connected_index_x] += 4*target_ineraction*sign_x;
                  }
                  
                  // y-direction
                  std::int32_t connected_index_y = coo_y - degree + 1 + i + j;
                  if (connected_index_y < 0) {
                     connected_index_y += y_size_;
                  }
                  else if (connected_index_y >= y_size_) {
                     connected_index_y -= y_size_;
                  }
                  if (connected_index_y != coo_y){
                     energy_difference_[connected_index_y*x_size_ + coo_x] += 4*target_ineraction*sign_y;
                  }
               }
            }
         }
      }
      else if (bc_ == lattice::BoundaryCondition::OBC) {
         for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
            if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
               continue;
            }
            const typename ModelType::ValueType target_ineraction = interaction_[degree];
            
            // x-direction
            for (std::int32_t i = std::max(coo_x - degree + 1, 0); i <= coo_x; ++i) {
               if (i > x_size_ - degree) {
                  break;
               }
               typename ModelType::OPType sign = 1;
               for (std::int32_t j = i; j < i + degree; ++j) {
                  sign *= sample_[coo_y*x_size_ + j];
               }
               for (std::int32_t j = i; j < i + degree; ++j) {
                  if (j == coo_x) {continue;}
                  energy_difference_[coo_y*x_size_ + j] += 4*target_ineraction*sign;
               }
            }
            
            // y-direction
            for (std::int32_t i = std::max(coo_y - degree + 1, 0); i <= coo_y; ++i) {
               if (i > y_size_ - degree) {
                  break;
               }
               typename ModelType::OPType sign = 1;
               for (std::int32_t j = i; j < i + degree; ++j) {
                  sign *= sample_[j*x_size_ + coo_x];
               }
               for (std::int32_t j = i; j < i + degree; ++j) {
                  if (j == coo_y) {continue;}
                  energy_difference_[j*x_size_ + coo_x] += 4*target_ineraction*sign;
               }
            }
         }
      }
      else {
         throw std::runtime_error("Unsupported BoundaryCondition");
      }
      energy_difference_[index] *= -1;
      sample_[index] *= -1;
      /*
      if (interaction_.size() == 4) {
         //const typename ModelType::ValueType J_2 = interaction_[2];
         const typename ModelType::ValueType J_3 = interaction_[3];
         if (bc_ == lattice::BoundaryCondition::PBC) {
            // x-direction
            const std::int32_t y_ind = coo_y*x_size_;
            if (coo_x == 0) {
               const typename ModelType::OPType sign_1 = sample_[y_ind + x_size_ - 2]*sample_[y_ind + x_size_ - 1]*sample_[y_ind];
               energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_1;
               energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_1;
               //energy_difference_[y_ind              ] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[y_ind + x_size_ - 1]*sample_[y_ind]*sample_[y_ind + 1];
               energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_2;
               //energy_difference_[y_ind              ] += 4*J_3*sign_2;
               energy_difference_[y_ind + 1          ] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[y_ind]*sample_[y_ind + 1]*sample_[y_ind + 2];
               //energy_difference_[y_ind              ] += 4*J_3*sign_3;
               energy_difference_[y_ind + 1          ] += 4*J_3*sign_3;
               energy_difference_[y_ind + 2          ] += 4*J_3*sign_3;
            }
            else if (coo_x == 1) {
               const typename ModelType::OPType sign_1 = sample_[y_ind + x_size_ - 1]*sample_[y_ind]*sample_[y_ind + 1];
               energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_1;
               energy_difference_[y_ind              ] += 4*J_3*sign_1;
               //energy_difference_[y_ind + 1          ] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[y_ind]*sample_[y_ind + 1]*sample_[y_ind + 2];
               energy_difference_[y_ind              ] += 4*J_3*sign_2;
               //energy_difference_[y_ind + 1          ] += 4*J_3*sign_2;
               energy_difference_[y_ind + 2          ] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[y_ind + 1]*sample_[y_ind + 2]*sample_[y_ind + 3];
               //energy_difference_[y_ind + 1          ] += 4*J_3*sign_3;
               energy_difference_[y_ind + 2          ] += 4*J_3*sign_3;
               energy_difference_[y_ind + 3          ] += 4*J_3*sign_3;
            }
            else if (coo_x == x_size_ - 1) {
               const typename ModelType::OPType sign_1 = sample_[y_ind + x_size_ - 3]*sample_[y_ind + x_size_ - 2]*sample_[y_ind + x_size_ - 1];
               energy_difference_[y_ind + x_size_ - 3] += 4*J_3*sign_1;
               energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_1;
               //energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[y_ind + x_size_ - 2]*sample_[y_ind + x_size_ - 1]*sample_[y_ind];
               energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_2;
               //energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_2;
               energy_difference_[y_ind              ] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[y_ind + x_size_ - 1]*sample_[y_ind]*sample_[y_ind + 1];
               //energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_3;
               energy_difference_[y_ind              ] += 4*J_3*sign_3;
               energy_difference_[y_ind + 1          ] += 4*J_3*sign_3;
            }
            else if (coo_x == x_size_ - 2) {
               const typename ModelType::OPType sign_1 = sample_[y_ind + x_size_ - 4]*sample_[y_ind + x_size_ - 3]*sample_[y_ind + x_size_ - 2];
               energy_difference_[y_ind + x_size_ - 4] += 4*J_3*sign_1;
               energy_difference_[y_ind + x_size_ - 3] += 4*J_3*sign_1;
               //energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[y_ind + x_size_ - 3]*sample_[y_ind + x_size_ - 2]*sample_[y_ind + x_size_ - 1];
               energy_difference_[y_ind + x_size_ - 3] += 4*J_3*sign_2;
               //energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_2;
               energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[y_ind + x_size_ - 2]*sample_[y_ind + x_size_ - 1]*sample_[y_ind];
               //energy_difference_[y_ind + x_size_ - 2] += 4*J_3*sign_3;
               energy_difference_[y_ind + x_size_ - 1] += 4*J_3*sign_3;
               energy_difference_[y_ind              ] += 4*J_3*sign_3;
            }
            else {
               const typename ModelType::OPType sign_1 = sample_[y_ind + coo_x - 2]*sample_[y_ind + coo_x - 1]*sample_[y_ind + coo_x];
               energy_difference_[y_ind + coo_x - 2] += 4*J_3*sign_1;
               energy_difference_[y_ind + coo_x - 1] += 4*J_3*sign_1;
               //energy_difference_[y_ind + coo_x    ] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[y_ind + coo_x - 1]*sample_[y_ind + coo_x]*sample_[y_ind + coo_x + 1];
               energy_difference_[y_ind + coo_x - 1] += 4*J_3*sign_2;
               //energy_difference_[y_ind + coo_x    ] += 4*J_3*sign_2;
               energy_difference_[y_ind + coo_x + 1] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[y_ind + coo_x]*sample_[y_ind + coo_x + 1]*sample_[y_ind + coo_x + 2];
               //energy_difference_[y_ind + coo_x    ] += 4*J_3*sign_3;
               energy_difference_[y_ind + coo_x + 1] += 4*J_3*sign_3;
               energy_difference_[y_ind + coo_x + 2] += 4*J_3*sign_3;
            }
            // y-direction
            if (coo_y == 0) {
               const typename ModelType::OPType sign_1 = sample_[(coo_y - 2)*x_size_ + coo_x]*sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_x];
               energy_difference_[(coo_y - 2)*x_size_ + coo_x] += 4*J_3*sign_1;
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += 4*J_3*sign_1;
               //energy_difference_[                      coo_x] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_x]*sample_[x_size_ + coo_x];
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += 4*J_3*sign_2;
               //energy_difference_[                      coo_x] += 4*J_3*sign_2;
               energy_difference_[            x_size_ + coo_x] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[coo_x]*sample_[x_size_ + coo_x]*sample_[2*x_size_ + coo_x];
               //energy_difference_[            coo_x] += 4*J_3*sign_3;
               energy_difference_[  x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[2*x_size_ + coo_x] += 4*J_3*sign_3;
            }
            else if (coo_y == 1) {
               const typename ModelType::OPType sign_1 = sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_x]*sample_[x_size_ + coo_x];
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += 4*J_3*sign_1;
               energy_difference_[                      coo_x] += 4*J_3*sign_1;
               //energy_difference_[            x_size_ + coo_x] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[coo_x]*sample_[x_size_ + coo_x]*sample_[2*x_size_ + coo_x];
               energy_difference_[            coo_x] += 4*J_3*sign_2;
               //energy_difference_[  x_size_ + coo_x] += 4*J_3*sign_2;
               energy_difference_[2*x_size_ + coo_x] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[x_size_ + coo_x]*sample_[2*x_size_ + coo_x]*sample_[3*x_size_ + coo_x];
               //energy_difference_[  x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[2*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[3*x_size_ + coo_x] += 4*J_3*sign_3;
            }
            else if (coo_y == y_size_ - 1) {
               const typename ModelType::OPType sign_1 = sample_[(y_size_ - 3)*x_size_ + coo_x]*sample_[(y_size_ - 2)*x_size_ + coo_x]*sample_[(y_size_ - 1)*x_size_ + coo_x];
               energy_difference_[(y_size_ - 3)*x_size_ + coo_x] += 4*J_3*sign_1;
               energy_difference_[(y_size_ - 2)*x_size_ + coo_x] += 4*J_3*sign_1;
               //energy_difference_[(y_size_ - 1)*x_size_ + coo_x] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[(y_size_ - 2)*x_size_ + coo_x]*sample_[(y_size_ - 1)*x_size_ + coo_x]*sample_[coo_x];
               energy_difference_[(y_size_ - 2)*x_size_ + coo_x] += 4*J_3*sign_2;
               //energy_difference_[(y_size_ - 1)*x_size_ + coo_x] += 4*J_3*sign_2;
               energy_difference_[                        coo_x] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[(y_size_ - 1)*x_size_ + coo_x]*sample_[coo_x]*sample_[x_size_ + coo_x];
               //energy_difference_[(y_size_ - 1)*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[                        coo_x] += 4*J_3*sign_3;
               energy_difference_[              x_size_ + coo_x] += 4*J_3*sign_3;
            }
            else if (coo_y == y_size_ - 2) {
               const typename ModelType::OPType sign_1 = sample_[(y_size_ - 4)*x_size_ + coo_x]*sample_[(y_size_ - 3)*x_size_ + coo_x]*sample_[(y_size_ - 2)*x_size_ + coo_x];
               energy_difference_[(y_size_ - 4)*x_size_ + coo_x] += 4*J_3*sign_1;
               energy_difference_[(y_size_ - 3)*x_size_ + coo_x] += 4*J_3*sign_1;
               //energy_difference_[(y_size_ - 2)*x_size_ + coo_x] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[(y_size_ - 3)*x_size_ + coo_x]*sample_[(y_size_ - 2)*x_size_ + coo_x]*sample_[(y_size_ - 1)*x_size_ + coo_x];
               energy_difference_[(y_size_ - 3)*x_size_ + coo_x] += 4*J_3*sign_2;
               //energy_difference_[(y_size_ - 2)*x_size_ + coo_x] += 4*J_3*sign_2;
               energy_difference_[(y_size_ - 1)*x_size_ + coo_x] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[(y_size_ - 2)*x_size_ + coo_x]*sample_[(y_size_ - 1)*x_size_ + coo_x]*sample_[coo_x];
               //energy_difference_[(y_size_ - 2)*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[(y_size_ - 1)*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[                        coo_x] += 4*J_3*sign_3;
            }
            else {
               const typename ModelType::OPType sign_1 = sample_[(coo_y - 2)*x_size_ + coo_x]*sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_y*x_size_ + coo_x];
               energy_difference_[(coo_y - 2)*x_size_ + coo_x] += 4*J_3*sign_1;
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += 4*J_3*sign_1;
               //energy_difference_[      coo_y*x_size_ + coo_x] += 4*J_3*sign_1;
               const typename ModelType::OPType sign_2 = sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_y*x_size_ + coo_x]*sample_[(coo_y + 1)*x_size_ + coo_x];
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += 4*J_3*sign_2;
               //energy_difference_[      coo_y*x_size_ + coo_x] += 4*J_3*sign_2;
               energy_difference_[(coo_y + 1)*x_size_ + coo_x] += 4*J_3*sign_2;
               const typename ModelType::OPType sign_3 = sample_[coo_y*x_size_ + coo_x]*sample_[(coo_y + 1)*x_size_ + coo_x]*sample_[(coo_y + 2)*x_size_ + coo_x];
               //energy_difference_[      coo_y*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[(coo_y + 1)*x_size_ + coo_x] += 4*J_3*sign_3;
               energy_difference_[(coo_y + 2)*x_size_ + coo_x] += 4*J_3*sign_3;
            }
         }
         else if (bc_ == lattice::BoundaryCondition::OBC) {
            // x-direction
            const std::int32_t y_ind = coo_y*x_size_;
            if (2 <= coo_x) {
               typename ModelType::ValueType val = 4*J_3*sample_[y_ind + coo_x - 2]*sample_[y_ind + coo_x - 1]*sample_[y_ind + coo_x];
               energy_difference_[y_ind + coo_x - 2] += val;
               energy_difference_[y_ind + coo_x - 1] += val;
               //energy_difference_[y_ind + coo_x    ] += val;
            }
            if (1 <= coo_x && coo_x <= x_size_ - 2) {
               typename ModelType::ValueType val = 4*J_3*sample_[y_ind + coo_x - 1]*sample_[y_ind + coo_x]*sample_[y_ind + coo_x + 1];
               energy_difference_[y_ind + coo_x - 1] += val;
               //energy_difference_[y_ind + coo_x    ] += val;
               energy_difference_[y_ind + coo_x + 1] += val;
            }
            if (coo_x <= x_size_ - 3) {
               typename ModelType::ValueType val = 4*J_3*sample_[y_ind + coo_x]*sample_[y_ind + coo_x + 1]*sample_[y_ind + coo_x + 2];
               //energy_difference_[y_ind + coo_x    ] += val;
               energy_difference_[y_ind + coo_x + 1] += val;
               energy_difference_[y_ind + coo_x + 2] += val;
            }
            // y-direction
            if (2 <= coo_y) {
               typename ModelType::ValueType val = 4*J_3*sample_[(coo_y - 2)*x_size_ + coo_x]*sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_y*x_size_ + coo_x];
               energy_difference_[(coo_y - 2)*x_size_ + coo_x] += val;
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += val;
               //energy_difference_[            x_size_ + coo_x] += val;
            }
            if (1 <= coo_y && coo_y <= y_size_ - 2) {
               typename ModelType::ValueType val = 4*J_3*sample_[(coo_y - 1)*x_size_ + coo_x]*sample_[coo_y*x_size_ + coo_x]*sample_[(coo_y + 1)*x_size_ + coo_x];
               energy_difference_[(coo_y - 1)*x_size_ + coo_x] += val;
               //energy_difference_[            x_size_ + coo_x] += val;
               energy_difference_[(coo_y + 1)*x_size_ + coo_x] += val;
            }
            if (coo_y <= y_size_ - 3) {
               typename ModelType::ValueType val = 4*J_3*sample_[coo_y*x_size_ + coo_x]*sample_[(coo_y + 1)*x_size_ + coo_x]*sample_[(coo_y + 2)*x_size_ + coo_x];
               //energy_difference_[            x_size_ + coo_x] += val;
               energy_difference_[(coo_y + 1)*x_size_ + coo_x] += val;
               energy_difference_[(coo_y + 2)*x_size_ + coo_x] += val;
            }
         }
         else {
            throw std::runtime_error("Unsupported BoundaryCondition");
         }
         energy_difference_[index] *= -1;
         sample_[index] *= -1;
      }
      else {
         if (bc_ == lattice::BoundaryCondition::PBC) {
            for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
               if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
                  continue;
               }
               const typename ModelType::ValueType target_ineraction = interaction_[degree];
               
               for (std::int32_t i = 0; i < degree; ++i) {
                  typename ModelType::OPType sign_x = 1;
                  typename ModelType::OPType sign_y = 1;
                  for (std::int32_t j = 0; j < degree; ++j) {
                     // x-direction
                     std::int32_t connected_index_x = coo_x - degree + 1 + i + j;
                     if (connected_index_x < 0) {
                        connected_index_x += x_size_;
                     }
                     else if (connected_index_x >= x_size_) {
                        connected_index_x -= x_size_;
                     }
                     sign_x *= sample_[coo_y*x_size_ + connected_index_x];
                     
                     // y-direction
                     std::int32_t connected_index_y = coo_y - degree + 1 + i + j;
                     if (connected_index_y < 0) {
                        connected_index_y += y_size_;
                     }
                     else if (connected_index_y >= y_size_) {
                        connected_index_y -= y_size_;
                     }
                     sign_y *= sample_[connected_index_y*x_size_ + coo_x];
                  }
                  for (std::int32_t j = 0; j < degree; ++j) {
                     // x-direction
                     std::int32_t connected_index_x = coo_x - degree + 1 + i + j;
                     if (connected_index_x < 0) {
                        connected_index_x += x_size_;
                     }
                     else if (connected_index_x >= x_size_) {
                        connected_index_x -= x_size_;
                     }
                     if (connected_index_x != coo_x){
                        energy_difference_[coo_y*x_size_ + connected_index_x] += 4*target_ineraction*sign_x;
                     }
                     
                     // y-direction
                     std::int32_t connected_index_y = coo_y - degree + 1 + i + j;
                     if (connected_index_y < 0) {
                        connected_index_y += y_size_;
                     }
                     else if (connected_index_y >= y_size_) {
                        connected_index_y -= y_size_;
                     }
                     if (connected_index_y != coo_y){
                        energy_difference_[connected_index_y*x_size_ + coo_x] += 4*target_ineraction*sign_y;
                     }
                  }
               }
            }
         }
         else if (bc_ == lattice::BoundaryCondition::OBC) {
            for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
               if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
                  continue;
               }
               const typename ModelType::ValueType target_ineraction = interaction_[degree];
               
               // x-direction
               for (std::int32_t i = std::max(coo_x - degree + 1, 0); i <= coo_x; ++i) {
                  if (i > x_size_ - degree) {
                     break;
                  }
                  typename ModelType::OPType sign = 1;
                  for (std::int32_t j = i; j < i + degree; ++j) {
                     sign *= sample_[coo_y*x_size_ + j];
                  }
                  for (std::int32_t j = i; j < i + degree; ++j) {
                     if (j == coo_x) {continue;}
                     energy_difference_[coo_y*x_size_ + j] += 4*target_ineraction*sign;
                  }
               }
               
               // y-direction
               for (std::int32_t i = std::max(coo_y - degree + 1, 0); i <= coo_y; ++i) {
                  if (i > y_size_ - degree) {
                     break;
                  }
                  typename ModelType::OPType sign = 1;
                  for (std::int32_t j = i; j < i + degree; ++j) {
                     sign *= sample_[j*x_size_ + coo_x];
                  }
                  for (std::int32_t j = i; j < i + degree; ++j) {
                     if (j == coo_y) {continue;}
                     energy_difference_[j*x_size_ + coo_x] += 4*target_ineraction*sign;
                  }
               }
            }
         }
         else {
            throw std::runtime_error("Unsupported BoundaryCondition");
         }
         energy_difference_[index] *= -1;
         sample_[index] *= -1;
      }
       */
   }
   
   const std::vector<typename ModelType::OPType> &GetSample() const {
      return sample_;
   }
   
   typename ModelType::ValueType GetEnergyDifference(const std::int32_t index) const {
      return energy_difference_[index];
   }
   
   std::int32_t GetSystemSize() const {
      return system_size_;
   }
   
private:
   const std::int32_t system_size_;
   const std::int32_t x_size_;
   const std::int32_t y_size_;
   const lattice::BoundaryCondition bc_;
   const std::vector<typename ModelType::ValueType> interaction_;
   
   std::vector<typename ModelType::OPType> sample_;
   std::vector<typename ModelType::ValueType> energy_difference_;
   
   std::vector<typename ModelType::ValueType> GenerateEnergyDifference(const std::vector<typename ModelType::OPType> &sample) const {
      std::vector<typename ModelType::ValueType> energy_difference(system_size_);
      if (bc_ == lattice::BoundaryCondition::PBC) {
         for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
            if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
               continue;
            }
            const typename ModelType::ValueType target_ineraction = interaction_[degree];

            for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
               for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
                  typename ModelType::ValueType val = 0;
                  for (std::int32_t i = 0; i < degree; ++i) {
                     typename ModelType::OPType sign_x = 1;
                     typename ModelType::OPType sign_y = 1;
                     for (std::int32_t j = 0; j < degree; ++j) {
                        // x-direction
                        std::int32_t connected_index_x = coo_x - degree + 1 + i + j;
                        if (connected_index_x < 0) {
                           connected_index_x += x_size_;
                        }
                        else if (connected_index_x >= x_size_) {
                           connected_index_x -= x_size_;
                        }
                        sign_x *= sample[coo_y*x_size_ + connected_index_x];
                        
                        // y-direction
                        std::int32_t connected_index_y = coo_y - degree + 1 + i + j;
                        if (connected_index_y < 0) {
                           connected_index_y += y_size_;
                        }
                        else if (connected_index_y >= y_size_) {
                           connected_index_y -= y_size_;
                        }
                        sign_y *= sample[connected_index_y*x_size_ + coo_x];
                     }
                     val += sign_x*target_ineraction + sign_y*target_ineraction;
                  }
                  energy_difference[coo_y*x_size_ + coo_x] = -2.0*val;
               }
            }
         }
      }
      else if (bc_ == lattice::BoundaryCondition::OBC) {
         for (std::int32_t degree = 1; degree < interaction_.size(); ++degree) {
            if (std::abs(interaction_[degree]) <= std::numeric_limits<typename ModelType::ValueType>::epsilon()) {
               continue;
            }
            const typename ModelType::ValueType target_ineraction = interaction_[degree];

            for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
               for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
                  // x-direction
                  typename ModelType::ValueType val = 0;
                  for (std::int32_t i = 0; i < degree; ++i) {
                     if (coo_x - degree + 1 + i < 0 || coo_x + i >= x_size_) {
                        continue;
                     }
                     typename ModelType::OPType sign = 1;
                     for (std::int32_t j = 0; j < degree; ++j) {
                        std::int32_t connected_index = coo_x - degree + 1 + i + j;
                        sign *= (sample)[coo_y*x_size_ + connected_index];
                     }
                     val += sign*target_ineraction;
                  }
                  
                  // y-direction
                  for (std::int32_t i = 0; i < degree; ++i) {
                     if (coo_y - degree + 1 + i < 0 || coo_y + i >= y_size_) {
                        continue;
                     }
                     typename ModelType::OPType sign = 1;
                     for (std::int32_t j = 0; j < degree; ++j) {
                        std::int32_t connected_index = coo_y - degree + 1 + i + j;
                        sign *= (sample)[connected_index*x_size_ + coo_x];
                     }
                     val += sign*target_ineraction;
                  }
                  energy_difference[coo_y*x_size_ + coo_x] = -2.0*val;
               }
            }
         }
      }
      else {
         throw std::runtime_error("Unsupported BinaryCondition");
      }
      return energy_difference;
   }
   
   std::vector<typename ModelType::OPType> GenerateRandomSpin(const std::uint64_t seed) const {
      std::vector<typename ModelType::OPType> sample(system_size_);
      std::uniform_int_distribution<utility::SpinType> dist(0, 1);
      utility::RandType random_number_engine(seed);
      for (std::size_t i = 0; i < sample.size(); i++) {
         sample[i] = 2*dist(random_number_engine) - 1;
      }
      return sample;
   }
   
};

template<typename RealType>
CMCSystem(const model::PolynomialIsing<lattice::Square, RealType>, const uint64_t) -> CMCSystem<model::PolynomialIsing<lattice::Square, RealType>>;


} // namespace cmc_utility
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_POLY_ISING_SQUARE_HPP_ */
