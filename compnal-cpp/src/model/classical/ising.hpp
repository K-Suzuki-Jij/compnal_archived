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
   
   Ising(const LatticeType &lattice,
         const RealType interaction_deg_1,
         const RealType interaction_deg_2):
   lattice_(lattice), interaction_deg_1_(interaction_deg_1), interaction_deg_2_(interaction_deg_2) {}
   
   std::tuple<RealType, RealType, RealType> GetInteraction() const {
      return {interaction_deg_0_, interaction_deg_1_, interaction_deg_2_};
   }
   
   std::int32_t GetSystemSize() const {
      return lattice_.GetSystemSize();
   }
   
   lattice::BoundaryCondition GetBoundaryCondition() const {
      return lattice_.GetBoundaryCondition();
   }
      
   std::int32_t GetDegree() const {
      return 2;
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
   
   RealType CalculateOnsiteSampleAverage(const std::vector<std::vector<OPType>> &samples,
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
   RealType interaction_deg_0_ = 0;
   RealType interaction_deg_1_ = 0;
   RealType interaction_deg_2_ = 0;

   
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
   
   RealType CalculateEnergy(const lattice::Honeycomb &honeycomb_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateEnergy(const lattice::Square &square_lattice,
                            const std::vector<OPType> &sample) const {
      RealType energy = 0;
      return energy;
   }
   
   RealType CalculateEnergy(const lattice::Triangle &triangle_lattice,
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
   
   Ising(const lattice::AnyLattice &lattice,
         const std::unordered_map<IndexType, RealType, IndexHash> &linear,
         const std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash> &quadratic):
   lattice_(lattice), interaction_(linear, quadratic) {}
   
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
                const RealType interaction_deg_1,
                const RealType interaction_deg_2) {
   return Ising<LatticeType, RealType>{lattice, interaction_deg_1, interaction_deg_2};
}

template<typename RealType>
auto make_ising(const lattice::AnyLattice &lattice,
                const std::unordered_map<
                typename QuadraticGeneralModel<RealType>::IndexType, RealType,
                typename QuadraticGeneralModel<RealType>::IndexHash
                > &linear,
                const std::unordered_map<
                std::pair<
                typename QuadraticGeneralModel<RealType>::IndexType,
                typename QuadraticGeneralModel<RealType>::IndexType
                >,
                RealType,
                typename QuadraticGeneralModel<RealType>::PairHash
                > &quadratic) {
   
   return Ising<lattice::AnyLattice, RealType>{lattice, linear, quadratic};
   
}

} // namespace model
} // namespace compnal

#endif /* COMPNAL_MODEL_ISING_HPP_ */
