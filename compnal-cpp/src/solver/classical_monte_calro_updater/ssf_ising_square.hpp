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
//  ssf_ising_square.hpp
//  compnal
//
//  Created by kohei on 2022/09/25.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_ISING_SQUARE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_ISING_SQUARE_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::Ising<lattice::Square, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::Ising<lattice::Square, RealType>::OPType> &sample,
                         const model::Ising<lattice::Square, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   const typename model::Ising<lattice::Square, RealType>::LinearType linear = model.GetLinear();
   const typename model::Ising<lattice::Square, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      // x-direction
      for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
         for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
            (*energy_difference)[coo_y*x_size + coo_x] += -2*quadratic*sample[coo_y*x_size + coo_x]*sample[coo_y*x_size + coo_x + 1] - 2*linear*sample[coo_y*x_size + coo_x];
            (*energy_difference)[coo_y*x_size + coo_x + 1] += -2*quadratic*sample[coo_y*x_size + coo_x]*sample[coo_y*x_size + coo_x + 1];
         }
         (*energy_difference)[coo_y*x_size + x_size - 1] += -2*linear*sample[coo_y*x_size + x_size - 1];
      }
      // y-direction
      for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
         for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
            (*energy_difference)[coo_y*x_size + coo_x] += -2*quadratic*sample[coo_y*x_size + coo_x]*sample[(coo_y + 1)*x_size + coo_x] - 2*linear*sample[coo_y*x_size + coo_x];
            (*energy_difference)[(coo_y + 1)*x_size + coo_x] += -2*quadratic*sample[coo_y*x_size + coo_x]*sample[(coo_y + 1)*x_size + coo_x];
         }
         (*energy_difference)[(y_size - 1)*x_size + coo_x] += -2*linear*sample[(y_size - 1)*x_size + coo_x];
      }
   }
   else {
      throw std::runtime_error("Unsupported BoundaryCondition");
   }
}


template<typename RealType>
void UpdateConfiguration(std::vector<typename model::Ising<lattice::Square, RealType>::OPType> *sample,
                         std::vector<typename model::Ising<lattice::Square, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::Ising<lattice::Square, RealType> &model) {
   
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   const std::int32_t coo_x = index%x_size;
   const std::int32_t coo_y = index/x_size;
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      // x-direction
      if (0 < coo_x && coo_x < x_size - 1) {
         (*energy_difference)[coo_y*x_size + coo_x - 1] += 4*quadratic*(*sample)[coo_y*x_size + coo_x - 1]*(*sample)[coo_y*x_size + coo_x];
         (*energy_difference)[coo_y*x_size + coo_x + 1] += 4*quadratic*(*sample)[coo_y*x_size + coo_x]*(*sample)[coo_y*x_size + coo_x + 1];
      }
      else if (coo_x == 0) {
         (*energy_difference)[coo_y*x_size + coo_x + 1] += 4*quadratic*(*sample)[coo_y*x_size + coo_x]*(*sample)[coo_y*x_size + coo_x + 1];
      }
      else {
         (*energy_difference)[coo_y*x_size + coo_x - 1] += 4*quadratic*(*sample)[coo_y*x_size + coo_x - 1]*(*sample)[coo_y*x_size + coo_x];
      }
      
      // y-direction
      if (0 < coo_y && coo_y < y_size - 1) {
         (*energy_difference)[(coo_y - 1)*x_size + coo_x] += 4*quadratic*(*sample)[(coo_y - 1)*x_size + coo_x]*(*sample)[coo_y*x_size + coo_x];
         (*energy_difference)[(coo_y + 1)*x_size + coo_x] += 4*quadratic*(*sample)[coo_y*x_size + coo_x]*(*sample)[(coo_y + 1)*x_size + coo_x];
      }
      else if (coo_y == 0) {
         (*energy_difference)[(coo_y + 1)*x_size + coo_x] += 4*quadratic*(*sample)[coo_y*x_size + coo_x]*(*sample)[(coo_y + 1)*x_size + coo_x];
      }
      else {
         (*energy_difference)[(coo_y - 1)*x_size + coo_x] += 4*quadratic*(*sample)[(coo_y - 1)*x_size + coo_x]*(*sample)[coo_y*x_size + coo_x];
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


#endif /* COMPNAL_SOLVER_UPDATER_SSF_ISING_SQUARE_HPP_ */
