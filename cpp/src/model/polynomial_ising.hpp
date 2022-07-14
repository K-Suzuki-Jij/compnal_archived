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


template<typename RealType>
class PolynomialIsing {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   
   using ValueType = RealType;
   using SpinType = char;
   
   PolynomialIsing(const int system_size, const std::vector<RealType> &interaction, const Lattice lattice) {
      SetSystemSize(system_size);
      SetInteraction(interaction);
      lattice_ = lattice;
   }
   
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         throw std::runtime_error("system_size mut be larger than 0");
      }
      system_size_ = system_size;
   }
   
   void SetInteraction(const std::vector<RealType> &interaction) {
      interaction_ = interaction;
      polynomial_degree_ = static_cast<int>(interaction.size());
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
   

private:
   int system_size_ = 0;
   int polynomial_degree_;
   std::vector<RealType> interaction_ = {};
   Lattice lattice_ = Lattice::NONE;
   
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_ */
