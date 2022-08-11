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
//  Created by Kohei Suzuki on 2022/06/10.
//

#ifndef COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_
#define COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_

#include "lattice.hpp"
#include <vector>

namespace compnal {
namespace model {

//! @brief Polynomial Ising class.
template<typename RealType>
class PolynomialIsing {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   
   using ValueType = RealType;
   using SpinType = std::int8_t;
   
   PolynomialIsing(const int system_size,
                   const std::vector<RealType> &interaction,
                   const Lattice lattice) {
      SetSystemSize(system_size);
      SetInteraction(interaction);
      SetLattice(lattice);
      SetBoundaryCondtion(GetDefaultBoundaryCondition(lattice));
   }
   
   PolynomialIsing(const int system_size,
                   const std::vector<RealType> &interaction,
                   const Lattice lattice,
                   const BoundaryCondition boundary_condition) {
      SetSystemSize(system_size);
      SetInteraction(interaction);
      SetLattice(lattice);
      SetBoundaryCondtion(boundary_condition);
   }
   
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         throw std::runtime_error("system_size mut be larger than 0");
      }
      system_size_ = system_size;
   }
   
   void SetInteraction(const std::vector<RealType> &interaction) {
      std::vector<RealType> new_interaction = interaction;
      std::int64_t size = static_cast<std::int64_t>(new_interaction.size());
      for (std::int64_t i = size - 1; i >= 0; --i) {
         if (std::abs(new_interaction[i]) > std::numeric_limits<RealType>::epsilon()) {
            new_interaction.resize(i + 1);
            break;
         }
      }
      interaction_ = new_interaction;
      polynomial_degree_ = static_cast<int>(new_interaction.size());
   }
   
   void SetLattice(const Lattice lattice) {
      lattice_ = lattice;
   }
   
   void SetBoundaryCondtion(const BoundaryCondition boundary_condition) {
      boundary_condition_ = boundary_condition;
   }
   
   int GetSystemSize() const {
      return system_size_;
   }
   
   const std::vector<RealType> &GetInteraction() const {
      return interaction_;
   }
   
   Lattice GetLattice() const {
      return lattice_;
   }
   
   RealType CalculateMagnetization(const std::vector<SpinType> &sample) const {
      if (static_cast<int>(sample.size()) != system_size_) {
         throw std::runtime_error("Sample size is not equal to system_size");
      }
      RealType magnetization = 0;
      for (const auto &it: sample) {
         magnetization += it;
      }
      return magnetization/system_size_;
   }
   
   RealType CalculateEnergy(const std::vector<SpinType> &sample) const {
      
      if (static_cast<int>(sample.size()) != system_size_) {
         throw std::runtime_error("Sample size is not equal to system_size");
      }
      
      RealType energy = 0;
      
      if (lattice_ == Lattice::CHAIN) {
         throw std::runtime_error("CHAIN is under construction.");
      }
      else if (lattice_ == Lattice::SQUARE) {
         throw std::runtime_error("SQUARE is under construction.");
      }
      else if (lattice_ == Lattice::TRIANGLE) {
         throw std::runtime_error("TRIANGLE is under construction.");
      }
      else if (lattice_ == Lattice::HONEYCOMB) {
         throw std::runtime_error("HONEYCOMB is under construction.");
      }
      else if (lattice_ == Lattice::CUBIC) {
         throw std::runtime_error("CUBIC is under construction.");
      }
      else if (lattice_ == Lattice::INFINIT_RANGE) {
         for (std::size_t i = 0; i < interaction_.size(); ++i) {
            const int polynomial_degree = static_cast<int>(i) + 1;
            const RealType target_ineraction = interaction_[i];
            if (std::abs(target_ineraction) > std::numeric_limits<RealType>::epsilon()) {
               std::vector<int> indices(polynomial_degree);
               int start_index = 0;
               int size = 0;
               while (true) {
                  for (int i = start_index; i < system_size_; ++i) {
                     indices[size++] = i;
                     if (size == polynomial_degree) {
                        SpinType sign = 1;
                        for (int j = 0; j < polynomial_degree; ++j) {
                           sign *= sample[indices[j]];
                        }
                        energy += target_ineraction*sign;
                        break;
                     }
                  }
                  --size;
                  if (size < 0) {
                     break;
                  }
                  start_index = indices[size] + 1;
               }
            }
         }
         return energy;
      }
      else if (lattice_ == Lattice::ANY_TYPE) {
         throw std::runtime_error("ANY_TYPE is under construction.");
      }
      else {
         throw std::runtime_error("Unknown lattice type detected.");
      }
   }

private:
   int system_size_ = 0;
   int polynomial_degree_ = 0;
   std::vector<RealType> interaction_ = {};
   Lattice lattice_ = Lattice::NONE;
   BoundaryCondition boundary_condition_ = BoundaryCondition::NONE;
   
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_ */
