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

#include "../../lattice/all.hpp"
#include "../base_polynomial_general_model.hpp"
#include <vector>
#include <unordered_map>

namespace compnal {
namespace model {

//! @brief Polynomial Ising class.
template<typename RealType, class LatticeType>
class PolynomialIsing {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   
   using ValueType = RealType;
   using SpinType = std::int8_t;
   
   const LatticeType lattice;
   
   PolynomialIsing(const LatticeType &input_lattice,
                   const std::unordered_map<int, RealType> &interaction
                   ): lattice(input_lattice) {
      for (const auto &it: interaction) {
         SetInteraction(it.first, it.second);
      }
   }
      
   void SetInteraction(const int degree, const RealType value) {
      if (std::abs(value) > std::numeric_limits<RealType>::epsilon()) {
         interaction_[degree] = value;
      }
   }
   
   const std::unordered_map<int, RealType> &GetInteraction() const {
      return interaction_;
   }
      
   RealType CalculateEnergy(const std::vector<SpinType> &sample) const {
      return CalculateEnergy(lattice, sample);
   }

private:
   std::unordered_map<int, RealType> interaction_;
   
   RealType CalculateEnergy(const lattice::InfiniteRange &infinite_range_lattice,
                            const std::vector<SpinType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
};

template<typename RealType>
class PolynomialIsing<RealType, lattice::AnyLattice> {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
  
public:
   using IndexType = typename BasePolynomialGeneralModel<RealType>::IndexType;
   using IndexHash = typename BasePolynomialGeneralModel<RealType>::IndexHash;
   using VectorHash = typename BasePolynomialGeneralModel<RealType>::VectorHash;
   using InteractionType = typename BasePolynomialGeneralModel<RealType>::InteractionType;
   
   const lattice::AnyLattice lattice;
   
   PolynomialIsing(const lattice::AnyLattice &input_lattice): lattice(input_lattice) {}
   
   void AddInteraction(const std::vector<IndexType> &index_list, const RealType value) {
      interaction.AddInteraction(index_list, value);
   }
   
   void GenerateIndexList() const {
      return interaction.GenerateIndexList();
   }
   
   const InteractionType &GetInteraction() const {
      return interaction.GetInteraction();
   }
   
   const std::unordered_set<IndexType, IndexHash> &GetIndexSet() const {
      return interaction.GetIndexSet();
   }
   
   const std::size_t GetSystemSize() const {
      return interaction.GetSystemSize();
   }
   
private:
   BasePolynomialGeneralModel<RealType> interaction;
   
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_ */
