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
//  ssf_polynomial_ising_cubic.hpp
//  compnal
//
//  Created by kohei on 2022/08/16.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CUBIC_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CUBIC_HPP_

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::PolynomialIsing<lattice::Cubic, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::PolynomialIsing<lattice::Cubic, RealType>::OPType> &sample,
                         const model::PolynomialIsing<lattice::Cubic, RealType> &model) {
   
}

template<typename RealType>
void UpdateConfiguration(std::vector<typename model::PolynomialIsing<lattice::Cubic, RealType>::OPType> *sample,
                         std::vector<typename model::PolynomialIsing<lattice::Cubic, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::Cubic, RealType> &model) {
   
}

} // namespace updater
} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_CUBIC_HPP_ */
