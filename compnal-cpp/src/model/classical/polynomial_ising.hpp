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
#include "../polynomial_general_model.hpp"
#include "../../utility/type.hpp"
#include <vector>
#include <unordered_map>

namespace compnal {
namespace model {

//! @brief Polynomial Ising class.
template<class LatticeType, typename RealType>
class PolynomialIsing {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   
   using ValueType = RealType;
   using OPType = utility::SpinType;
   
   const LatticeType lattice;
   
   PolynomialIsing(const LatticeType &input_lattice,
                   const std::unordered_map<std::int32_t, RealType> &interaction
                   ): lattice(input_lattice) {
      for (const auto &it: interaction) {
         SetInteraction(it.first, it.second);
      }
   }
      
   void SetInteraction(const std::int32_t degree, const RealType value) {
      if (degree < 0) {
         throw std::runtime_error("degree must be larger than or equal to 0");
      }
      if (std::abs(value) > std::numeric_limits<RealType>::epsilon()) {
         interaction_.resize(degree + 1);
         interaction_[degree] = value;
      }
   }
   
   const std::vector<RealType> &GetInteraction() const {
      return interaction_;
   }
      
   RealType CalculateEnergy(const std::vector<OPType> &sample) const {
      return CalculateEnergy(lattice, sample);
   }
   
   std::int32_t GetSystemSize() const {
      return lattice.GetSystemSize();
   }
   
   static RealType CalculateAverage(const std::vector<std::vector<OPType>> &samples) {
      return CalculateMoment(samples, 1);
   }
   
   std::int32_t GetDegree() const {
      return static_cast<std::int32_t>(interaction_.size()) - 1;
   }
   
   static RealType CalculateMoment(const std::vector<std::vector<OPType>> &samples, const std::int32_t degree) {
      RealType val = RealType{0.0};
      if (degree == 1) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (std::int32_t i = 0; i < static_cast<std::int32_t>(samples.size()); ++i) {
            const auto &sample = samples[i];
            RealType temp = RealType{0.0};
            for (std::int32_t j = 0; j < static_cast<std::int32_t>(sample.size()); ++j) {
               temp += sample[j];
            }
            val += temp/sample.size();
         }
      }
      else if (degree == 2) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (std::int32_t i = 0; i < static_cast<std::int32_t>(samples.size()); ++i) {
            const auto &sample = samples[i];
            const std::int32_t sample_size = static_cast<std::int32_t>(sample.size());
            RealType temp = RealType{0.0};
            for (std::int32_t j_1 = 0; j_1 < sample_size; ++j_1) {
               for (std::int32_t j_2 = 0; j_2 < sample_size; ++j_2) {
                  temp += sample[j_1]*sample[j_2];
               }
            }
            val += temp/(sample_size*sample_size);
         }
      }
      else if (degree == 4) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (std::int32_t i = 0; i < static_cast<std::int32_t>(samples.size()); ++i) {
            const auto &sample = samples[i];
            const std::int32_t sample_size = static_cast<std::int32_t>(sample.size());
            RealType temp = RealType{0.0};
            for (std::int32_t j_1 = 0; j_1 < sample_size; ++j_1) {
               for (std::int32_t j_2 = 0; j_2 < sample_size; ++j_2) {
                  for (std::int32_t j_3 = 0; j_3 < sample_size; ++j_3) {
                     for (std::int32_t j_4 = 0; j_4 < sample_size; ++j_4) {
                        temp += sample[j_1]*sample[j_2]*sample[j_3]*sample[j_4];
                     }
                  }
               }
            }
            val += temp/(sample_size*sample_size*sample_size*sample_size);
         }
      }
      else {
         throw std::runtime_error("Under construction");
      }
      return val/samples.size();
   }

private:
   std::vector<RealType> interaction_;
   
   RealType CalculateEnergy(const lattice::InfiniteRange &infinite_range_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
};

template<class LatticeType, typename Key, typename RealType>
PolynomialIsing(LatticeType, std::unordered_map<Key, RealType>)
-> PolynomialIsing<LatticeType, RealType>;

template<typename RealType>
class PolynomialIsing<lattice::AnyLattice, RealType> {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
  
public:
   using ValueType = RealType;
   using OPType = utility::SpinType;
   using IndexType = typename PolynomialGeneralModel<RealType>::IndexType;
   using IndexHash = typename PolynomialGeneralModel<RealType>::IndexHash;
   using VectorHash = typename PolynomialGeneralModel<RealType>::VectorHash;
   using InteractionType = typename PolynomialGeneralModel<RealType>::InteractionType;
   
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
   PolynomialGeneralModel<RealType> interaction;
   
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_ISING_HPP_ */
