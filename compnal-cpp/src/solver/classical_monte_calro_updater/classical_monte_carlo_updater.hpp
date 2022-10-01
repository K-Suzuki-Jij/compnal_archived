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
//  classical_monte_carlo_updater.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_CLASSICAL_MONTE_CARLO_UPDATER_HPP_
#define COMPNAL_SOLVER_UPDATER_CLASSICAL_MONTE_CARLO_UPDATER_HPP_

#include "ssf_polynomial_ising_any_lattice.hpp"
#include "ssf_polynomial_ising_chain.hpp"
#include "ssf_polynomial_ising_cubic.hpp"
#include "ssf_polynomial_ising_square.hpp"
#include "ssf_polynomial_ising_infinite_range.hpp"

#include "ssf_ising_chain.hpp"
#include "ssf_ising_square.hpp"
#include "ssf_ising_cubic.hpp"
#include "ssf_ising_infinite_range.hpp"
#include "ssf_ising_any_lattice.hpp"

#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<class ModelType>
void ExecuteMetropolis(std::vector<typename ModelType::OPType> *sample,
                       const ModelType &model,
                       const std::int32_t num_sweeps,
                       const typename ModelType::ValueType beta,
                       const std::uint64_t seed) {
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }

   // Set energy difference
   std::vector<typename ModelType::ValueType> energy_difference(model.GetSystemSize());
   SetEnergyDifference<typename ModelType::ValueType>(&energy_difference, *sample, model);
   
   const std::int32_t system_size = model.GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<typename ModelType::ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
      
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (energy_difference[index] <= 0.0 || std::exp(-beta*energy_difference[index]) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, index, model);
         }
      }
   }
}

template<typename RealType>
void ExecuteMetropolis(std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> *sample,
                       const model::PolynomialIsing<lattice::AnyLattice, RealType> &model,
                       const std::int32_t num_sweeps,
                       const typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType beta,
                       const std::uint64_t seed) {
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType;
   const std::int32_t system_size = model.GetSystemSize();
   std::vector<std::int8_t> sign_list(model.GetKeyList().size());
   std::vector<ValueType> energy_difference(model.GetSystemSize());
   
   for (std::size_t i = 0; i < model.GetKeyList().size(); ++i) {
      std::int8_t sign = 1;
      for (const auto &index: model.GetKeyList()[i]) {
         sign *= (*sample)[index];
      }
      sign_list[i] = sign;
   }

   // Set energy difference
   SetEnergyDifference<ValueType>(&energy_difference, *sample, sign_list, model);
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
   
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (energy_difference[index] <= 0.0 || std::exp(-beta*energy_difference[index]) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, &sign_list, index, model.GetAdjacencyList(), model.GetKeyList(), model.GetValueList());
         }
      }
   }
   
}

template<typename RealType>
void ExecuteMetropolis(std::vector<typename model::Ising<lattice::AnyLattice, RealType>::OPType> *sample,
                       const model::Ising<lattice::AnyLattice, RealType> &model,
                       const std::int32_t num_sweeps,
                       const typename model::Ising<lattice::AnyLattice, RealType>::ValueType beta,
                       const std::uint64_t seed) {
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::Ising<lattice::AnyLattice, RealType>::OPType;
   using ValueType = typename model::Ising<lattice::AnyLattice, RealType>::ValueType;

   // Set energy difference
   std::vector<ValueType> energy_difference(model.GetSystemSize());
   SetEnergyDifference<ValueType>(&energy_difference, *sample, model);
   
   const std::int32_t system_size = model.GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
      
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (energy_difference[index] <= 0.0 || std::exp(-beta*energy_difference[index]) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, index, model.GetRowPtr(), model.GetColPtr(), model.GetValPtr());
         }
      }
   }
   
}

template<class ModelType>
void ExecuteHeatBath(std::vector<typename ModelType::OPType> *sample,
                     const ModelType &model,
                     const std::int32_t num_sweeps,
                     const typename ModelType::ValueType beta,
                     const std::uint64_t seed) {
      
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   // Set energy difference
   std::vector<typename ModelType::ValueType> energy_difference(model.GetSystemSize());
   SetEnergyDifference<typename ModelType::ValueType>(&energy_difference, *sample, model);

   const std::int32_t system_size = model.GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<typename ModelType::ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
   
   // Do Heat Bath update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (1/(1 + std::exp(beta*(energy_difference)[index])) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, index, model);
         }
      }
   }
}

template<typename RealType>
void ExecuteHeatBath(std::vector<typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType> *sample,
                     const model::PolynomialIsing<lattice::AnyLattice, RealType> &model,
                     const std::int32_t num_sweeps,
                     const typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType beta,
                     const std::uint64_t seed) {
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::AnyLattice, RealType>::ValueType;
   const std::int32_t system_size = model.GetSystemSize();
   std::vector<std::int8_t> sign_list(model.GetKeyList().size());
   std::vector<ValueType> energy_difference(model.GetSystemSize());
   
   for (std::size_t i = 0; i < model.GetKeyList().size(); ++i) {
      std::int8_t sign = 1;
      for (const auto &index: model.GetKeyList()[i]) {
         sign *= (*sample)[index];
      }
      sign_list[i] = sign;
   }

   // Set energy difference
   SetEnergyDifference<ValueType>(&energy_difference, *sample, sign_list, model);
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
   
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (1/(1 + std::exp(beta*(energy_difference)[index])) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, &sign_list, index, model.GetAdjacencyList(), model.GetKeyList(), model.GetValueList());
         }
      }
   }
   
}

template<typename RealType>
void ExecuteHeatBath(std::vector<typename model::Ising<lattice::AnyLattice, RealType>::OPType> *sample,
                     const model::Ising<lattice::AnyLattice, RealType> &model,
                     const std::int32_t num_sweeps,
                     const typename model::Ising<lattice::AnyLattice, RealType>::ValueType beta,
                     const std::uint64_t seed) {
   
   if (sample->size() != static_cast<std::size_t>(model.GetSystemSize())) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::Ising<lattice::AnyLattice, RealType>::OPType;
   using ValueType = typename model::Ising<lattice::AnyLattice, RealType>::ValueType;

   // Set energy difference
   std::vector<ValueType> energy_difference(model.GetSystemSize());
   SetEnergyDifference<ValueType>(&energy_difference, *sample, model);
   
   const std::int32_t system_size = model.GetSystemSize();
   
   // Set random number engine
   utility::RandType random_number_engine(seed);
   std::uniform_real_distribution<ValueType> dist_real(0, 1);
   std::uniform_int_distribution<std::int32_t> dist_system_size(0, system_size - 1);
      
   // Do Metropolis update
   for (std::int32_t sweep_count = 0; sweep_count < num_sweeps; sweep_count++) {
      for (std::int32_t i = 0; i < system_size; i++) {
         const std::int32_t index = dist_system_size(random_number_engine);
         if (1/(1 + std::exp(beta*(energy_difference)[index])) > dist_real(random_number_engine)) {
            UpdateConfiguration(sample, &energy_difference, index, model.GetRowPtr(), model.GetColPtr(), model.GetValPtr());
         }
      }
   }
   
   
}


} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_CLASSICAL_MONTE_CARLO_UPDATER_HPP_ */
