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
//  Created by Kohei Suzuki on 2022/07/22.
//

#ifndef COMPNAL_SOLVER_UPDATER_POLY_ISING_CHAIN_HPP_
#define COMPNAL_SOLVER_UPDATER_POLY_ISING_CHAIN_HPP_

namespace compnal {
namespace solver {

template<typename SpinType, typename RealType>
void ResetEnergyDifferenceChain(std::vector<RealType> *energy_difference,
                                std::vector<int> *system_size_index_list,
                                const std::vector<SpinType> &sample,
                                const model::PolynomialIsing<RealType> &model) {
   
   if (static_cast<int>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   const int system_size = model.GetSystemSize();
   const auto &interaction = model.GetInteraction();
   energy_difference->clear();
   energy_difference->resize(system_size);
   
   
}




} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_POLY_ISING_CHAIN_HPP_ */
