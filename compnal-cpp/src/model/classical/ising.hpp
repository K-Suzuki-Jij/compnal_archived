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
//  ising.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_MODEL_ISING_HPP_
#define COMPNAL_MODEL_ISING_HPP_

#include "../../lattice/all.hpp"
#include "../../utility/type.hpp"
#include "../quadratic_general_model.hpp"
#include <vector>

namespace compnal {
namespace model {

template<class LatticeType, typename RealType>
class Ising {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   
   using ValueType = RealType;
   using IndexType = std::int32_t;
   using OPType = utility::SpinType;
   using LinearType = RealType;
   using QuadraticType = RealType;
   
   Ising(const LatticeType &lattice,
         const LinearType interaction_deg_1,
         const QuadraticType interaction_deg_2):
   lattice_(lattice), linear_(interaction_deg_1), quadratic_(interaction_deg_2) {}
   
   std::pair<LinearType, QuadraticType> GetInteraction() const {
      return {linear_, quadratic_};
   }
   
   LinearType GetLinear() const {
      return linear_;
   }
   
   QuadraticType GetQuadratic() const {
      return quadratic_;
   }
   
   std::int32_t GetSystemSize() const {
      return lattice_.GetSystemSize();
   }
   
   lattice::BoundaryCondition GetBoundaryCondition() const {
      return lattice_.GetBoundaryCondition();
   }
      
   std::int32_t GetDegree() const {
      if (std::abs(quadratic_) > std::numeric_limits<RealType>::epsilon()) {
         return 2;
      }
      else if (std::abs(linear_) > std::numeric_limits<RealType>::epsilon()) {
         return 1;
      }
      else {
         return 0;
      }
   }
   
   const LatticeType &GetLattice() const {
      return lattice_;
   }
   
   RealType CalculateEnergy(const std::vector<OPType> &sample) const {
      return CalculateEnergy(lattice_, sample);
   }
   
   RealType CalculateMoment(const std::vector<std::vector<OPType>> &samples,
                            const std::int32_t degree,
                            const std::int32_t num_threads = utility::DEFAULT_NUM_THREADS) const {
      if (degree <= 0) {
         throw std::runtime_error("degree must be lager than 0.");
      }
      
      RealType val = 0;
#pragma omp parallel for schedule(guided) reduction(+: val) num_threads(num_threads)
      for (std::int32_t i = 0; i < static_cast<std::int32_t>(samples.size()); ++i) {
         RealType avg = CalculateMagnetization(samples[i]);
         RealType prod = 1;
         for (std::int32_t j = 0; j < degree; ++j) {
            prod = prod*avg;
         }
         val += prod;
      }
      return val/samples.size();
   }
   
   RealType CalculateOnsiteAverage(const std::vector<std::vector<OPType>> &samples,
                                         const IndexType index) const {
      if (index < 0) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::size_t i = 0; i < samples.size(); ++i) {
         val += samples[i][index];
      }
      return val/samples.size();
   }
   
   RealType CalculateCorrelation(const std::vector<std::vector<OPType>> &samples,
                                 const IndexType ind1,
                                 const IndexType ind2) const {
      const std::int32_t size = static_cast<std::int32_t>(samples.size());
      if (ind1 < 0 || ind2 < 0) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::int32_t i = 0; i < size; ++i) {
         val += samples[i][ind1]*samples[i][ind2];
      }
      return val/size;
   }
   
private:
   LatticeType lattice_;
   LinearType linear_ = 0;
   QuadraticType quadratic_ = 0;

   
   RealType CalculateEnergy(const lattice::Chain &chain_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateEnergy(const lattice::Cubic &cubic_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateEnergy(const lattice::Square &square_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateEnergy(const lattice::InfiniteRange &infinite_range_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateMagnetization(const std::vector<OPType> &sample) const {
      RealType val = 0;
      for (std::size_t i = 0; i < sample.size(); ++i) {
         val += sample[i];
      }
      return val/sample.size();
   }
};

template<typename RealType>
class Ising<lattice::AnyLattice, RealType> {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   using ValueType = RealType;
   using OPType = utility::SpinType;
   using IndexType = typename QuadraticGeneralModel<RealType>::IndexType;
   using IndexHash = typename QuadraticGeneralModel<RealType>::IndexHash;
   using PairHash  = typename QuadraticGeneralModel<RealType>::PairHash;
   using LinearType = typename QuadraticGeneralModel<RealType>::LinearType;
   using QuadraticType = typename QuadraticGeneralModel<RealType>::QuadraticType;
   
   Ising(const lattice::AnyLattice &lattice,
         const LinearType &linear,
         const QuadraticType &quadratic):
   lattice_(lattice), interaction_(linear, quadratic) {}
      
   const std::vector<IndexType> &GetIndexList() const {
      return interaction_.GetIndexList();
   }
   
   RealType GetConstant() const {
      return interaction_.GetConstant();
   }
   
   const LinearType &GetLinear() const {
      return interaction_.GetLinear();
   }
   
   const std::vector<RealType> &GetQuadraticValue() const {
      return interaction_.GetQuadraticValue();
   }
   
   const std::vector<std::pair<std::int32_t, std::int32_t>> &GetQuadraticKey() const {
      return interaction_.GetQuadraticKey();
   }
   
   const std::unordered_map<IndexType, std::int32_t, IndexHash> &GetIndexMap() const {
      return interaction_.GetIndexMap();
   }
   
   std::int32_t GetSystemSize() const {
      return interaction_.GetSystemSize();
   }
   
   std::int32_t GetDegree() const {
      return interaction_.GetDegree();
   }
   
   lattice::BoundaryCondition GetBoundaryCondition() const {
      return lattice_.GetBoundaryCondition();
   }
   
   RealType CalculateEnergy(const std::vector<OPType> &sample) const {
      if (sample.size() != interaction_.GetSystemSize()) {
         throw std::runtime_error("The sample size is not equal to the system size");
      }
      const std::unordered_map<IndexType, std::int64_t, IndexHash> &index_map = interaction_.GetIndexMap();
      RealType val = interaction_.GetConstant();
      
      for (const auto &it: interaction_.GetLinearInteraction()) {
         val += it.second*sample[index_map.at(it.first)];
      }
      
      for (const auto &it: interaction_.GetQuadraticInteraction()) {
         val += it.second*sample[index_map.at(it.first.first)]*sample[index_map.at(it.first.second)];
      }
      
      return val;
   }
   
   RealType CalculateMoment(const std::vector<std::vector<OPType>> &samples,
                            const std::int32_t degree,
                            const std::int32_t num_threads = utility::DEFAULT_NUM_THREADS) const {
      if (degree <= 0) {
         throw std::runtime_error("degree must be lager than 0.");
      }
      
      RealType val = 0;
#pragma omp parallel for schedule(guided) reduction(+: val) num_threads(num_threads)
      for (std::int32_t i = 0; i < static_cast<std::int32_t>(samples.size()); ++i) {
         RealType avg = CalculateMagnetization(samples[i]);
         RealType prod = 1;
         for (std::int32_t j = 0; j < degree; ++j) {
            prod = prod*avg;
         }
         val += prod;
      }
      return val/samples.size();
   }
   
   RealType CalculateCorrelation(const std::vector<std::vector<OPType>> &samples,
                                 const IndexType ind1,
                                 const IndexType ind2) const {
      const std::unordered_map<IndexType, std::int64_t, IndexHash> &index_map = interaction_.GetIndexMap();
      if (index_map.count(ind1) == 0 || index_map.count(ind2) == 0) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::size_t i = 0; i < samples.size(); ++i) {
         val += samples[i][index_map.at(ind1)]*samples[i][index_map.at(ind2)];
      }
      return val/samples.size();
   }
   
private:
   QuadraticGeneralModel<RealType> interaction_;
   lattice::AnyLattice lattice_;
   
   RealType CalculateMagnetization(const std::vector<OPType> &sample) const {
      RealType val = 0;
      for (std::size_t i = 0; i < sample.size(); ++i) {
         val += sample[i];
      }
      return val/sample.size();
   }
   
};

template<class LatticeType, typename RealType>
auto make_ising(const LatticeType &lattice,
                const typename Ising<LatticeType, RealType>::LinearType &linear,
                const typename Ising<LatticeType, RealType>::QuadraticType &quadratic) {
   return Ising<LatticeType, RealType>{lattice, linear, quadratic};
}

} // namespace model
} // namespace compnal

#endif /* COMPNAL_MODEL_ISING_HPP_ */
