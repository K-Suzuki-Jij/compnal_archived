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
//  system_ising_cubic.hpp
//  compnal
//
//  Created by kohei on 2022/11/23.
//  
//

#ifndef COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_ISING_CUBIC_HPP_
#define COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_ISING_CUBIC_HPP_

#include "../../lattice/all.hpp"
#include "../../model/all.hpp"
#include "base_system.hpp"

namespace compnal {
namespace solver {
namespace cmc_utility {

template<typename RealType>
class CMCSystem<model::Ising<lattice::Cubic, RealType>> {
   
   using ModelType = model::Ising<lattice::Cubic, RealType>;

public:
   using ValueType = typename ModelType::ValueType;
   
   CMCSystem(const ModelType &model, const uint64_t seed):
   system_size_(model.GetSystemSize()),
   x_size_(model.GetLattice().GetXSize()),
   y_size_(model.GetLattice().GetYSize()),
   z_size_(model.GetLattice().GetZSize()),
   bc_(model.GetBoundaryCondition()),
   quadratic_(model.GetQuadratic()),
   linear_(model.GetLinear()) {
      sample_ = GenerateRandomSpin(seed);
      energy_difference_ = GenerateEnergyDifference(sample_);
   }
   
   void Flip(const std::int32_t index) {
      const std::int32_t coo_z = index/(x_size_*y_size_);
      const std::int32_t coo_y = (index - coo_z*x_size_*y_size_)/x_size_;
      const std::int32_t coo_x = index - coo_z*x_size_*y_size_ - coo_y*x_size_;
      if (bc_ == lattice::BoundaryCondition::PBC) {
         // x-direction
         if (0 < coo_x && coo_x < x_size_ - 1) {
            energy_difference_[index - 1] += 4*quadratic_*sample_[index - 1]*sample_[index];
            energy_difference_[index + 1] += 4*quadratic_*sample_[index]*sample_[index + 1];
         }
         else if (coo_x == 0) {
            energy_difference_[index + 1] += 4*quadratic_*sample_[index]*sample_[index + 1];
            energy_difference_[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1] += 4*quadratic_*sample_[index]*sample_[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1];
         }
         else {
            energy_difference_[coo_z*x_size_*y_size_ + coo_y*x_size_ + 0] += 4*quadratic_*sample_[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1]*sample_[coo_z*x_size_*y_size_ + coo_y*x_size_ + 0];
            energy_difference_[index - 1] += 4*quadratic_*sample_[index - 1]*sample_[index];
         }
         
         // y-direction
         const std::int32_t index_y_m1 = coo_z*x_size_*y_size_ + (coo_y - 1)*x_size_ + coo_x;
         const std::int32_t index_y_p1 = coo_z*x_size_*y_size_ + (coo_y + 1)*x_size_ + coo_x;
         if (0 < coo_y && coo_y < y_size_ - 1) {
            energy_difference_[index_y_m1] += 4*quadratic_*sample_[index_y_m1]*sample_[index];
            energy_difference_[index_y_p1] += 4*quadratic_*sample_[index]*sample_[index_y_p1];
         }
         else if (coo_y == 0) {
            energy_difference_[index_y_p1] += 4*quadratic_*sample_[index]*sample_[index_y_p1];
            energy_difference_[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x] += 4*quadratic_*sample_[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x]*sample_[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x];
         }
         else {
            energy_difference_[index_y_m1] += 4*quadratic_*sample_[index_y_m1]*sample_[index];
            energy_difference_[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x] += 4*quadratic_*sample_[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x]*sample_[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x];
         }
         
         // z-direction
         const std::int32_t index_z_m1 = (coo_z - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
         const std::int32_t index_z_p1 = (coo_z + 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
         if (0 < coo_z && coo_z < z_size_ - 1) {
            energy_difference_[index_z_m1] += 4*quadratic_*sample_[index_z_m1]*sample_[index];
            energy_difference_[index_z_p1] += 4*quadratic_*sample_[index]*sample_[index_z_p1];
         }
         else if (coo_z == 0) {
            energy_difference_[index_z_p1] += 4*quadratic_*sample_[index]*sample_[index_z_p1];
            energy_difference_[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x] += 4*quadratic_*sample_[0*x_size_*y_size_ + coo_y*x_size_ + coo_x]*sample_[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x];
         }
         else {
            energy_difference_[index_z_m1] += 4*quadratic_*sample_[index_z_m1]*sample_[index];
            energy_difference_[0*x_size_*y_size_ + coo_y*x_size_ + coo_x] += 4*quadratic_*sample_[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x]*sample_[0*x_size_*y_size_ + coo_y*x_size_ + coo_x];
         }
      }
      else if (bc_ == lattice::BoundaryCondition::OBC) {
         // x-direction
         if (0 < coo_x && coo_x < x_size_ - 1) {
            energy_difference_[index - 1] += 4*quadratic_*sample_[index - 1]*sample_[index];
            energy_difference_[index + 1] += 4*quadratic_*sample_[index]*sample_[index + 1];
         }
         else if (coo_x == 0) {
            energy_difference_[index + 1] += 4*quadratic_*sample_[index]*sample_[index + 1];
         }
         else {
            energy_difference_[index - 1] += 4*quadratic_*sample_[index - 1]*sample_[index];
         }
         
         // y-direction
         const std::int32_t index_y_m1 = coo_z*x_size_*y_size_ + (coo_y - 1)*x_size_ + coo_x;
         const std::int32_t index_y_p1 = coo_z*x_size_*y_size_ + (coo_y + 1)*x_size_ + coo_x;
         if (0 < coo_y && coo_y < y_size_ - 1) {
            energy_difference_[index_y_m1] += 4*quadratic_*sample_[index_y_m1]*sample_[index];
            energy_difference_[index_y_p1] += 4*quadratic_*sample_[index]*sample_[index_y_p1];
         }
         else if (coo_y == 0) {
            energy_difference_[index_y_p1] += 4*quadratic_*sample_[index]*sample_[index_y_p1];
         }
         else {
            energy_difference_[index_y_m1] += 4*quadratic_*sample_[index_y_m1]*sample_[index];
         }
         
         // z-direction
         const std::int32_t index_z_m1 = (coo_z - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
         const std::int32_t index_z_p1 = (coo_z + 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
         if (0 < coo_z && coo_z < z_size_ - 1) {
            energy_difference_[index_z_m1] += 4*quadratic_*sample_[index_z_m1]*sample_[index];
            energy_difference_[index_z_p1] += 4*quadratic_*sample_[index]*sample_[index_z_p1];
         }
         else if (coo_z == 0) {
            energy_difference_[index_z_p1] += 4*quadratic_*sample_[index]*sample_[index_z_p1];
         }
         else {
            energy_difference_[index_z_m1] += 4*quadratic_*sample_[index_z_m1]*sample_[index];
         }
      }
      else {
         throw std::runtime_error("Unsupported BoundaryCondition");
      }
      energy_difference_[index] *= -1;
      sample_[index] *= -1;
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
   const std::int32_t z_size_;
   const lattice::BoundaryCondition bc_;
   const typename ModelType::QuadraticType quadratic_;
   const typename ModelType::LinearType linear_;
   
   std::vector<typename ModelType::OPType> sample_;
   std::vector<typename ModelType::ValueType> energy_difference_;
   
   std::vector<typename ModelType::ValueType> GenerateEnergyDifference(const std::vector<typename ModelType::OPType> &sample) const {
      std::vector<typename ModelType::ValueType> energy_difference(system_size_);
      if (bc_ == lattice::BoundaryCondition::PBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size_; ++coo_z) {
               for (std::int32_t coo_x = 0; coo_x < x_size_ - 1; ++coo_x) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index + 1] - 2*linear_*sample[index];
                  energy_difference[index + 1] += -2*quadratic_*sample[index]*sample[index + 1];
               }
               energy_difference[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1] +=
               -2*quadratic_*sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1]*
               sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + 0] -2*linear_*
               sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1];
               
               energy_difference[coo_z*x_size_*y_size_ + coo_y*x_size_ + 0] +=
               -2*quadratic_*sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1]*
               sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + 0];
            }
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
            for (std::int32_t coo_z = 0; coo_z < z_size_; ++coo_z) {
               for (std::int32_t coo_y = 0; coo_y < y_size_ - 1; ++coo_y) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  const std::int32_t index_p1 = coo_z*x_size_*y_size_ + (coo_y + 1)*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index_p1] - 2*linear_*sample[index];
                  energy_difference[index_p1] += -2*quadratic_*sample[index]*sample[index_p1];
               }
               energy_difference[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x] +=
               -2*quadratic_*sample[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x]*
               sample[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x] -
               2*linear_*sample[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x];
               
               energy_difference[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x] +=
               -2*quadratic_*sample[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x]*
               sample[coo_z*x_size_*y_size_ + 0*x_size_ + coo_x];
            }
         }
         
         // z-direction
         for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
               for (std::int32_t coo_z = 0; coo_z < z_size_ - 1; ++coo_z) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  const std::int32_t index_p1 = (coo_z + 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index_p1] - 2*linear_*sample[index];
                  energy_difference[index_p1] += -2*quadratic_*sample[index]*sample[index_p1];
               }
               energy_difference[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x] += -2*quadratic_*sample[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x]*sample[0*x_size_*y_size_ + coo_y*x_size_ + coo_x] -2*linear_*sample[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x];
               energy_difference[0*x_size_*y_size_ + coo_y*x_size_ + coo_x] += -2*quadratic_*sample[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x]*sample[0*x_size_*y_size_ + coo_y*x_size_ + coo_x];
            }
         }
      }
      else if (bc_ == lattice::BoundaryCondition::OBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size_; ++coo_z) {
               for (std::int32_t coo_x = 0; coo_x < x_size_ - 1; ++coo_x) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index + 1] - 2*linear_*sample[index];
                  energy_difference[index + 1] += -2*quadratic_*sample[index]*sample[index + 1];
               }
               energy_difference[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1] += -2*linear_*sample[coo_z*x_size_*y_size_ + coo_y*x_size_ + x_size_ - 1];
            }
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
            for (std::int32_t coo_z = 0; coo_z < z_size_; ++coo_z) {
               for (std::int32_t coo_y = 0; coo_y < y_size_ - 1; ++coo_y) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  const std::int32_t index_p1 = coo_z*x_size_*y_size_ + (coo_y + 1)*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index_p1] - 2*linear_*sample[index];
                  energy_difference[index_p1] += -2*quadratic_*sample[index]*sample[index_p1];
               }
               energy_difference[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x] += -2*linear_*sample[coo_z*x_size_*y_size_ + (y_size_ - 1)*x_size_ + coo_x];
            }
         }
         
         // z-direction
         for (std::int32_t coo_x = 0; coo_x < x_size_; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size_; ++coo_y) {
               for (std::int32_t coo_z = 0; coo_z < z_size_ - 1; ++coo_z) {
                  const std::int32_t index = coo_z*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  const std::int32_t index_p1 = (coo_z + 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x;
                  energy_difference[index] += -2*quadratic_*sample[index]*sample[index_p1] - 2*linear_*sample[index];
                  energy_difference[index_p1] += -2*quadratic_*sample[index]*sample[index_p1];
               }
               energy_difference[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x] += -2*linear_*sample[(z_size_ - 1)*x_size_*y_size_ + coo_y*x_size_ + coo_x];
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
CMCSystem(const model::Ising<lattice::Cubic, RealType>, const uint64_t) -> CMCSystem<model::Ising<lattice::Cubic, RealType>>;


} // namespace cmc_utility
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_CMC_UTILITY_SYSTEM_ISING_CUBIC_HPP_ */
