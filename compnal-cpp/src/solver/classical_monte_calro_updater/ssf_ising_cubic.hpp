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
//  ssf_ising_cubic.hpp
//  compnal
//
//  Created by kohei on 2022/09/28.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_ISING_CUBIC_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_ISING_CUBIC_HPP_

#include "../../lattice/cubic.hpp"
#include "../../model/classical/ising.hpp"
#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::Ising<lattice::Cubic, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::Ising<lattice::Cubic, RealType>::OPType> &sample,
                         const model::Ising<lattice::Cubic, RealType> &model) {
   
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   if (energy_difference->size() != sample.size()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   const std::int32_t z_size = model.GetLattice().GetZSize();
   const typename model::Ising<lattice::Square, RealType>::LinearType linear = model.GetLinear();
   const typename model::Ising<lattice::Square, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      // x-direction
      for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
         for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
            for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index + 1] - 2*linear*sample[index];
               (*energy_difference)[index + 1] += -2*quadratic*sample[index]*sample[index + 1];
            }
            (*energy_difference)[coo_z*x_size*y_size + coo_y*x_size + x_size - 1] += -2*quadratic*sample[coo_z*x_size*y_size + coo_y*x_size + x_size - 1]*sample[coo_z*x_size*y_size + coo_y*x_size + 0] -2*linear*sample[coo_z*x_size*y_size + coo_y*x_size + x_size - 1];
            (*energy_difference)[coo_z*x_size*y_size + coo_y*x_size + 0] += -2*quadratic*sample[coo_z*x_size*y_size + coo_y*x_size + x_size - 1]*sample[coo_z*x_size*y_size + coo_y*x_size + 0];
         }
      }
      
      // y-direction
      for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
         for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
            for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               const std::int32_t index_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index_p1] - 2*linear*sample[index];
               (*energy_difference)[index_p1] += -2*quadratic*sample[index]*sample[index_p1];
            }
            (*energy_difference)[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x] += -2*quadratic*sample[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x]*sample[coo_z*x_size*y_size + 0*x_size + coo_x] -2*linear*sample[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x];
            (*energy_difference)[coo_z*x_size*y_size + 0*x_size + coo_x] += -2*quadratic*sample[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x]*sample[coo_z*x_size*y_size + 0*x_size + coo_x];
         }
      }
      
      // z-direction
      for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size - 1; ++coo_z) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               const std::int32_t index_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index_p1] - 2*linear*sample[index];
               (*energy_difference)[index_p1] += -2*quadratic*sample[index]*sample[index_p1];
            }
            (*energy_difference)[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x] += -2*quadratic*sample[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x]*sample[0*x_size*y_size + coo_y*x_size + coo_x] -2*linear*sample[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x];
            (*energy_difference)[0*x_size*y_size + coo_y*x_size + coo_x] += -2*quadratic*sample[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x]*sample[0*x_size*y_size + coo_y*x_size + coo_x];
         }
      }
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      // x-direction
      for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
         for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
            for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index + 1] - 2*linear*sample[index];
               (*energy_difference)[index + 1] += -2*quadratic*sample[index]*sample[index + 1];
            }
            (*energy_difference)[coo_z*x_size*y_size + coo_y*x_size + x_size - 1] += -2*linear*sample[coo_z*x_size*y_size + coo_y*x_size + x_size - 1];
         }
      }
      
      // y-direction
      for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
         for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
            for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               const std::int32_t index_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index_p1] - 2*linear*sample[index];
               (*energy_difference)[index_p1] += -2*quadratic*sample[index]*sample[index_p1];
            }
            (*energy_difference)[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x] += -2*linear*sample[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x];
         }
      }
      
      // z-direction
      for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size - 1; ++coo_z) {
               const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
               const std::int32_t index_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
               (*energy_difference)[index] += -2*quadratic*sample[index]*sample[index_p1] - 2*linear*sample[index];
               (*energy_difference)[index_p1] += -2*quadratic*sample[index]*sample[index_p1];
            }
            (*energy_difference)[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x] += -2*linear*sample[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x];
         }
      }
   }
   else {
      throw std::runtime_error("Unsupported BoundaryCondition");
   }

}


template<typename RealType>
void UpdateConfiguration(std::vector<typename model::Ising<lattice::Cubic, RealType>::OPType> *sample,
                         std::vector<typename model::Ising<lattice::Cubic, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::Ising<lattice::Cubic, RealType> &model) {
   
   const std::int32_t x_size = model.GetLattice().GetXSize();
   const std::int32_t y_size = model.GetLattice().GetYSize();
   const std::int32_t z_size = model.GetLattice().GetZSize();
   const std::int32_t coo_z = index/(x_size*y_size);
   const std::int32_t coo_y = (index - coo_z*x_size*y_size)/x_size;
   const std::int32_t coo_x = index - coo_z*x_size*y_size - coo_y*x_size;
   const typename model::Ising<lattice::Chain, RealType>::QuadraticType quadratic = model.GetQuadratic();
   
   if (model.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
      // x-direction
      if (0 < coo_x && coo_x < x_size - 1) {
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
      }
      else if (coo_x == 0) {
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
         (*energy_difference)[coo_z*x_size*y_size + coo_y*x_size + x_size - 1] += 4*quadratic*(*sample)[index]*(*sample)[coo_z*x_size*y_size + coo_y*x_size + x_size - 1];
      }
      else {
         (*energy_difference)[coo_z*x_size*y_size + coo_y*x_size + 0] += 4*quadratic*(*sample)[coo_z*x_size*y_size + coo_y*x_size + x_size - 1]*(*sample)[coo_z*x_size*y_size + coo_y*x_size + 0];
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
      }
      
      // y-direction
      const std::int32_t index_y_m1 = coo_z*x_size*y_size + (coo_y - 1)*x_size + coo_x;
      const std::int32_t index_y_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
      if (0 < coo_y && coo_y < y_size - 1) {
         (*energy_difference)[index_y_m1] += 4*quadratic*(*sample)[index_y_m1]*(*sample)[index];
         (*energy_difference)[index_y_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_y_p1];
      }
      else if (coo_y == 0) {
         (*energy_difference)[index_y_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_y_p1];
         (*energy_difference)[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x] += 4*quadratic*(*sample)[coo_z*x_size*y_size + 0*x_size + coo_x]*(*sample)[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x];
      }
      else {
         (*energy_difference)[index_y_m1] += 4*quadratic*(*sample)[index_y_m1]*(*sample)[index];
         (*energy_difference)[coo_z*x_size*y_size + 0*x_size + coo_x] += 4*quadratic*(*sample)[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x]*(*sample)[coo_z*x_size*y_size + 0*x_size + coo_x];
      }
      
      // z-direction
      const std::int32_t index_z_m1 = (coo_z - 1)*x_size*y_size + coo_y*x_size + coo_x;
      const std::int32_t index_z_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
      if (0 < coo_z && coo_z < z_size - 1) {
         (*energy_difference)[index_z_m1] += 4*quadratic*(*sample)[index_z_m1]*(*sample)[index];
         (*energy_difference)[index_z_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_z_p1];
      }
      else if (coo_z == 0) {
         (*energy_difference)[index_z_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_z_p1];
         (*energy_difference)[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x] += 4*quadratic*(*sample)[0*x_size*y_size + coo_y*x_size + coo_x]*(*sample)[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x];
      }
      else {
         (*energy_difference)[index_z_m1] += 4*quadratic*(*sample)[index_z_m1]*(*sample)[index];
         (*energy_difference)[0*x_size*y_size + coo_y*x_size + coo_x] += 4*quadratic*(*sample)[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x]*(*sample)[0*x_size*y_size + coo_y*x_size + coo_x];
      }
   }
   else if (model.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
      // x-direction
      if (0 < coo_x && coo_x < x_size - 1) {
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
      }
      else if (coo_x == 0) {
         (*energy_difference)[index + 1] += 4*quadratic*(*sample)[index]*(*sample)[index + 1];
      }
      else {
         (*energy_difference)[index - 1] += 4*quadratic*(*sample)[index - 1]*(*sample)[index];
      }
      
      // y-direction
      const std::int32_t index_y_m1 = coo_z*x_size*y_size + (coo_y - 1)*x_size + coo_x;
      const std::int32_t index_y_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
      if (0 < coo_y && coo_y < y_size - 1) {
         (*energy_difference)[index_y_m1] += 4*quadratic*(*sample)[index_y_m1]*(*sample)[index];
         (*energy_difference)[index_y_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_y_p1];
      }
      else if (coo_y == 0) {
         (*energy_difference)[index_y_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_y_p1];
      }
      else {
         (*energy_difference)[index_y_m1] += 4*quadratic*(*sample)[index_y_m1]*(*sample)[index];
      }
      
      // z-direction
      const std::int32_t index_z_m1 = (coo_z - 1)*x_size*y_size + coo_y*x_size + coo_x;
      const std::int32_t index_z_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
      if (0 < coo_z && coo_z < z_size - 1) {
         (*energy_difference)[index_z_m1] += 4*quadratic*(*sample)[index_z_m1]*(*sample)[index];
         (*energy_difference)[index_z_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_z_p1];
      }
      else if (coo_z == 0) {
         (*energy_difference)[index_z_p1] += 4*quadratic*(*sample)[index]*(*sample)[index_z_p1];
      }
      else {
         (*energy_difference)[index_z_m1] += 4*quadratic*(*sample)[index_z_m1]*(*sample)[index];
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


#endif /* COMPNAL_SOLVER_UPDATER_SSF_ISING_CUBIC_HPP_ */
