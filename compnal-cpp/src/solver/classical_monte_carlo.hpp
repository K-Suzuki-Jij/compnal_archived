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

const CMCUpdater DEFAULT_CMC_UPDATER = CMCUpdater::METROPOLIS;

template<class ModelType>
class ClassicalMonteCarlo {
   
public:

   using RealType  = typename ModelType::ValueType;
   using IndexType = typename ModelType::IndexType;
   using OPType    = typename ModelType::OPType;
   
   ClassicalMonteCarlo(const ModelType  &model,
                       const CMCUpdater cmc_updater = DEFAULT_CMC_UPDATER):
   cmc_updater_(cmc_updater), seed_(std::random_device()()), model_(model) {
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
   
   void SetNumThreads(const std::int32_t num_threads) {
      if (num_threads < 0) {
         throw std::runtime_error("num_threads must be non-negative integer.");
      }
      num_threads_ = num_threads;
   }
   
   void SetInverseTemperature(const RealType beta) {
      if (beta < 0) {
         throw std::runtime_error("Inverse temperature must be positive number");
      }
      beta_ = beta;
   }
   
   void SetTemperature(const RealType temperature) {
      if (temperature <= 0) {
         throw std::runtime_error("Temperature must be positive number");
      }
      beta_ = 1/temperature;
   }
   
   void SetCMCUpdater(const CMCUpdater cmc_updater) {
      cmc_updater_ = cmc_updater;
   }
   
   std::int32_t GetNumSweeps() const {
      return num_sweeps_;
   }
   
   std::int32_t GetNumSamples() const {
      return num_samples_;
   }
   
   std::int32_t GetNumThreads() const {
      return num_threads_;
   }
   
   RealType GetInverseTemperature() const {
      return beta_;
   }
   
   RealType GetTemperature() const {
      return 1/beta_;
   }
   
   CMCUpdater GetCMCUpdater() const {
      return cmc_updater_;
   }

   std::uint64_t GetSeed() const {
      return seed_;
   }
   
   const std::vector<std::vector<OPType>> &GetSamples() const {
      return samples_;
   }
   
   const std::vector<OPType> &GetSample(const std::size_t i) const {
      return samples_.at(i);
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
      
      samples_.clear();
      samples_.shrink_to_fit();
      samples_.resize(num_samples_);
      
      if (cmc_updater_ == CMCUpdater::METROPOLIS) {
#pragma omp parallel for schedule(guided) num_threads(num_threads_)
         for (std::int32_t sample_count = 0; sample_count < num_samples_; sample_count++) {
            samples_[sample_count].resize(model_.GetSystemSize());
            RandomizeConfiguration(&samples_[sample_count], configuration_seed_list[sample_count]);
            updater::ExecuteMetropolis(&samples_[sample_count], model_, num_sweeps_, beta_, execute_seed_list[sample_count]);
         }
      }
      else if (cmc_updater_ == CMCUpdater::HEAT_BATH) {
#pragma omp parallel for schedule(guided) num_threads(num_threads_)
         for (std::int32_t sample_count = 0; sample_count < num_samples_; sample_count++) {
            samples_[sample_count].resize(model_.GetSystemSize());
            RandomizeConfiguration(&samples_[sample_count], configuration_seed_list[sample_count]);
            updater::ExecuteHeatBath(&samples_[sample_count], model_, num_sweeps_, beta_, execute_seed_list[sample_count]);
         }
      }
      else {
         throw std::runtime_error("Unknown Updater");
      }

   }
   
   RealType CalculateSampleAverage() const {
      return model_.CalculateMoment(samples_, 1, num_threads_);
   }
   
   RealType CalculateSampleMoment(const std::int32_t degree) const {
      return model_.CalculateMoment(samples_, degree, num_threads_);
   }
   
   RealType CalculateCorrelation(const IndexType ind1, const IndexType ind2) const {
      return model_.CalculateCorrelation(samples_, ind1, ind2);
   }
   
   std::vector<RealType> CalculateCorrelationList(const IndexType origin,
                                                  const std::vector<IndexType> &index_list) const {
      std::int32_t size = static_cast<std::int32_t>(index_list.size());
      std::vector<RealType> value_list(size);
#pragma omp parallel for schedule(guided) num_threads(num_threads_)
      for (std::int32_t i = 0; i < size; ++i) {
         value_list[i] = model_.CalculateCorrelation(samples_, origin, index_list[i]);
      }
      return value_list;
   }
   
private:
   std::int32_t num_sweeps_  = utility::DEFAULT_NUM_SWEEPS;
   std::int32_t num_samples_ = utility::DEFAULT_NUM_SAMPLES;
   std::int32_t num_threads_ = utility::DEFAULT_NUM_THREADS;
   RealType beta_            = utility::DEFAULT_BETA<RealType>;
   CMCUpdater cmc_updater_   = DEFAULT_CMC_UPDATER;
   std::uint64_t seed_       = std::random_device()();
   ModelType model_;
   std::vector<std::vector<OPType>> samples_;

   void RandomizeConfiguration(std::vector<utility::SpinType> *sample,
                               const std::uint64_t seed) const {
      std::uniform_int_distribution<utility::SpinType> dist(0, 1);
      utility::RandType random_number_engine(seed);
      for (std::size_t i = 0; i < sample->size(); i++) {
         (*sample)[i] = 2*dist(random_number_engine) - 1;
      }
   }
      
};


template<class ModelType>
auto make_classical_monte_carlo(const ModelType  &model,
                                const CMCUpdater cmc_updater = DEFAULT_CMC_UPDATER) {
   return ClassicalMonteCarlo<ModelType>{model, cmc_updater};
}


} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
