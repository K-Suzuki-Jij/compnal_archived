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
//  Created by Kohei Suzuki on 2022/07/18.
//

#ifndef COMPNAL_SOLVER_UPDATER_ALL_HPP_
#define COMPNAL_SOLVER_UPDATER_ALL_HPP_

#include "updater_poly_ising_infinite_range.hpp"

namespace compnal {
namespace solver {

template<typename SpinType, typename RealType>
void ResetEnergyDifference(std::vector<RealType> *energy_difference,
                           std::vector<int> *system_size_index_list,
                           const std::vector<SpinType> &sample,
                           const model::PolynomialIsing<RealType> &model) {
   
   if (static_cast<int>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   if (model.GetLattice() == model::Lattice::CHAIN) {
      throw std::runtime_error("CHAIN is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::SQUARE) {
      throw std::runtime_error("SQUARE is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::TRIANGLE) {
      throw std::runtime_error("TRIANGLE is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::HONEYCOMB) {
      throw std::runtime_error("HONEYCOMB is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::CUBIC) {
      throw std::runtime_error("CUBIC is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::INFINIT_RANGE) {
      ResetEnergyDifferenceInfiniteRange(energy_difference, system_size_index_list, sample, model);
   }
   else if (model.GetLattice() == model::Lattice::ANY_TYPE) {
      throw std::runtime_error("ANY_TYPE is under construction.");
   }
   else {
      throw std::runtime_error("Unknown lattice type detected.");
   }
   
}

template<typename SpinType, typename RealType>
void UpdateConfiguration(std::vector<SpinType> *sample,
                         std::vector<RealType> *energy_difference,
                         std::vector<int> *system_size_index_list,
                         const int index,
                         const model::PolynomialIsing<RealType> &model) {
   
   if (model.GetLattice() == model::Lattice::CHAIN) {
      throw std::runtime_error("CHAIN is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::SQUARE) {
      throw std::runtime_error("SQUARE is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::TRIANGLE) {
      throw std::runtime_error("TRIANGLE is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::HONEYCOMB) {
      throw std::runtime_error("HONEYCOMB is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::CUBIC) {
      throw std::runtime_error("CUBIC is under construction.");
   }
   else if (model.GetLattice() == model::Lattice::INFINIT_RANGE) {
      UpdaterInfiniteRange(sample, energy_difference, system_size_index_list, index, model);
   }
   else if (model.GetLattice() == model::Lattice::ANY_TYPE) {
      throw std::runtime_error("ANY_TYPE is under construction.");
   }
   else {
      throw std::runtime_error("Unknown lattice type detected.");
   }
   
}



} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_UPDATER_ALL_HPP_ */
