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
//  Created by Kohei Suzuki on 2022/06/11.
//

#ifndef COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_
#define COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_

#include "../utility/type.hpp"
#include "updater/monte_carlo_updater.hpp"
#include <vector>
#include <random>
#include <sstream>

namespace compnal {
namespace solver {

enum class CMCUpdater {
   
   METROPOLIS,
   HEAT_BATH
   
};

template<class ModelType>
class ClassicalMonteCarlo {
   
   using RealType = typename ModelType::ValueType;
   using OPType = typename ModelType::OPType;
   
public:
   
   const ModelType model;
   const CMCUpdater cmc_updater;
   
   ClassicalMonteCarlo(const ModelType &model,
                       const CMCUpdater cmc_updater=CMCUpdater::METROPOLIS):
   model(model), cmc_updater(cmc_updater) {
      seed_ = std::random_device()();
      samples_.resize(num_samples_);
      for (std::size_t i = 0; i < samples_.size(); ++i) {
         RandomizeConfiguration(&samples_[i], seed_);
      }
   }
   
   void SetNumSweeps(const std::int32_t num_sweeps) {
      if (num_sweeps < 0) {
         throw std::runtime_error("num_sweeps must be non-negative integer.");
      }
      num_sweeps_ = num_sweeps;
   }
   
   void SetNumSamples(const std::int32_t num_samples) {
      if (num_samples < 0) {
         throw std::runtime_error("num_samples must be non-negative integer.");
      }
      num_samples_ = num_samples;
   }
   
   void SetTemperature(const RealType temperature) {
      if (temperature <= 0) {
         throw std::runtime_error("Temperature must be positive number");
      }
      beta_ = 1/temperature;
   }
   
   void SetInverseTemperature(const RealType beta) {
      if (beta < 0) {
         throw std::runtime_error("Inverse temperature must be positive number");
      }
      beta_ = beta;
   }
   
   std::int32_t GetNumSweeps() const {
      return num_sweeps_;
   }
   
   std::int32_t GetNumSamples() const {
      return num_samples_;
   }
   
   const std::vector<std::vector<OPType>> &GetSamples() const {
      return samples_;
   }
   
   const std::vector<OPType> &GetSample(const std::int64_t i) const {
      return samples_.at(i);
   }
   
   RealType GetTemperature() const {
      return 1/beta_;
   }
   
   RealType GetInverseTemperature() const {
      return beta_;
   }
   
   std::uint64_t GetSeed() const {
      return seed_;
   }
   
   void Run() {
      Run(std::random_device()());
   }
   
   void Run(const std::uint64_t seed) {
      if (num_samples_ < 0) {
         throw std::runtime_error("num_samples must be non-negative integer.");
      }
      if (num_sweeps_ < 0) {
         throw std::runtime_error("num_sweeps must be non-negative integer.");
      }
      
      seed_ = seed;
      
      utility::RandType random_number_engine(seed_);
      std::vector<std::uint64_t> configuration_seed_list(num_samples_);
      std::vector<std::uint64_t> execute_seed_list(num_samples_);

      for (int32_t i = 0; i < num_samples_; ++i) {
         configuration_seed_list[i] = random_number_engine();
         execute_seed_list[i] = random_number_engine();
      }
      
      samples_.resize(num_samples_);
      for (std::size_t i = 0; i < samples_.size(); ++i) {
         samples_[i].resize(model.GetSystemSize());
      }
      samples_.shrink_to_fit();
      
      if (this->cmc_updater == CMCUpdater::METROPOLIS) {
#pragma omp parallel for schedule(guided)
         for (std::int32_t sample_count = 0; sample_count < num_samples_; sample_count++) {
            std::vector<std::pair<OPType, RealType>> sample_energy_difference_pair(model.GetSystemSize());
            RandomizeConfiguration(&sample_energy_difference_pair, configuration_seed_list[sample_count]);
            updater::ExecuteMetropolis(&sample_energy_difference_pair, model, num_sweeps_, beta_, execute_seed_list[sample_count]);
            for (std::size_t i = 0; i < sample_energy_difference_pair.size(); ++i) {
               samples_[sample_count][i] = sample_energy_difference_pair[i].first;
            }
         }
      }
      else if (this->cmc_updater == CMCUpdater::HEAT_BATH) {
#pragma omp parallel for schedule(guided)
         for (std::int32_t sample_count = 0; sample_count < num_samples_; sample_count++) {
            std::vector<std::pair<OPType, RealType>> sample_energy_difference_pair(model.GetSystemSize());
            RandomizeConfiguration(&sample_energy_difference_pair, configuration_seed_list[sample_count]);
            updater::ExecuteHeatBath(&sample_energy_difference_pair, model, num_sweeps_, beta_, execute_seed_list[sample_count]);
            for (std::size_t i = 0; i < sample_energy_difference_pair.size(); ++i) {
               samples_[sample_count][i] = sample_energy_difference_pair[i].first;
            }
         }
      }
      else {
         throw std::runtime_error("Unknown Updater");
      }

   }
   
   RealType CalculateSampleAverage() const {
      return model.CalculateAverage(samples_);
   }
   
   RealType CalculateSampleMoment(const std::int32_t degree) const {
      return model.CalculateMoment(samples_, degree);
   }
   
private:
   std::int32_t num_sweeps_  = 1000;
   std::int32_t num_samples_ = 1;
   RealType beta_   = 1;
   std::vector<std::vector<OPType>> samples_;
   std::uint64_t seed_ = 0;
   
   void RandomizeConfiguration(std::vector<utility::SpinType> *sample,
                               const std::uint64_t seed) const {
      std::uniform_int_distribution<utility::SpinType> dist(0, 1);
      sample->resize(model.GetSystemSize());
      utility::RandType random_number_engine(seed);
      for (std::size_t i = 0; i < sample->size(); i++) {
         (*sample)[i] = 2*dist(random_number_engine) - 1;
      }
   }
   
   void RandomizeConfiguration(std::vector<std::pair<OPType, RealType>> *sample_energy_difference_pair,
                               const std::uint64_t seed) const {
      std::uniform_int_distribution<utility::SpinType> dist(0, 1);
      utility::RandType random_number_engine(seed);
      for (auto &&it: (*sample_energy_difference_pair)) {
         it.first = 2*dist(random_number_engine) - 1;
      }
   }
   
};


template<class ModelType>
auto make_classical_monte_carlo(const ModelType &model,
                                const CMCUpdater cmc_updater=CMCUpdater::METROPOLIS) {
   return ClassicalMonteCarlo<ModelType>{model, cmc_updater};
}


} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
