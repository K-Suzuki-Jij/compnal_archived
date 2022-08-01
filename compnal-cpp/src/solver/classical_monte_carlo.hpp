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

#include "../model/polynomial_ising.hpp"
#include "classical_monte_carlo_updater/updater_all.hpp"
#include <vector>
#include <random>
#include <sstream>

namespace compnal {
namespace solver {

enum class Updater {
   
   METROPOLIS,
   HEAT_BATH
   
};


template<class ModelType>
class ClassicalMonteCarlo {
   
   using RealType = typename ModelType::ValueType;
   using SpinType = typename ModelType::SpinType;
   
public:
   
   const ModelType model;
   const Updater updater;
   
   ClassicalMonteCarlo(const ModelType &model_input,
                       const Updater updater_input):
   model(model_input), updater(updater_input) {}
   
   void SetNumSweeps(const int num_sweeps) {
      if (num_sweeps < 0) {
         throw std::runtime_error("num_sweeps must be non-negative integer.");
      }
      num_sweeps_ = num_sweeps;
   }
   
   void SetNumSamples(const int num_samples) {
      if (num_samples < 0) {
         throw std::runtime_error("num_samples must be non-negative integer.");
      }
      num_samples_ = num_samples;
      SetRandomSeed();
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
   
   void SetSeed(const std::uint64_t seed) {
      if (seed < 0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "seed=" << seed << "is not allowed." << std::endl;
         ss << "seed must be non-negative integer.";
         throw std::runtime_error(ss.str());
      }
      seed_list_.resize(num_samples_);
      std::mt19937_64 mt(seed);
      for (std::size_t i = 0; i < seed_list_.size(); ++i) {
         seed_list_[i] = mt();
      }
   }
   
   void SetRandomSeed() {
      SetSeed(std::random_device()());
   }
   
   int GetNumSweeps() const {
      return num_sweeps_;
   }
   
   int GetNumSamples() const {
      return num_samples_;
   }
   
   const std::vector<std::vector<SpinType>> &GetSamples() const {
      return samples_;
   }
   
   const std::vector<RealType> &GetEnergies() const {
      return energies_;
   }
   
   const std::vector<RealType> &GetMagnetizations() const {
      return magnetizations_;
   }
   
   template<typename T>
   const std::vector<SpinType> &GetSample(const T index) const {
      return samples_.at(index);
   }
   
   void ClearResults() {
      samples_.clear();
      samples_.shrink_to_fit();
      energies_.clear();
      energies_.shrink_to_fit();
      magnetizations_.clear();
      magnetizations_.shrink_to_fit();
   }
   
   void Run() {
      SetRandomSeed();
      Sample(model);
      energies_ = CalculateEnergies();
      magnetizations_ = CalculateMagnetizations();
   }
   
   void Run(const std::uint64_t seed) {
      SetSeed(seed);
      Sample(model);
      energies_ = CalculateEnergies();
      magnetizations_ = CalculateMagnetizations();
   }
   
   RealType CalculateSampleAverage(const int degree = 1) const {
      RealType val = RealType{0.0};
      if (degree == 1) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (int i = 0; i < static_cast<int>(samples_.size()); ++i) {
            const auto &sample = samples_[i];
            RealType temp = RealType{0.0};
            for (int j = 0; j < static_cast<int>(sample.size()); ++j) {
               temp += sample[j];
            }
            val += temp/sample.size();
         }
      }
      else if (degree == 2) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (int i = 0; i < static_cast<int>(samples_.size()); ++i) {
            const auto &sample = samples_[i];
            const int sample_size = static_cast<int>(sample.size());
            RealType temp = RealType{0.0};
            for (int j_1 = 0; j_1 < sample_size; ++j_1) {
               for (int j_2 = 0; j_2 < sample_size; ++j_2) {
                  temp += sample[j_1]*sample[j_2];
               }
            }
            val += temp/(sample_size*sample_size);
         }
      }
      else if (degree == 4) {
#pragma omp parallel for schedule(guided) reduction(+: val)
         for (int i = 0; i < static_cast<int>(samples_.size()); ++i) {
            const auto &sample = samples_[i];
            const int sample_size = static_cast<int>(sample.size());
            RealType temp = RealType{0.0};
            for (int j_1 = 0; j_1 < sample_size; ++j_1) {
               for (int j_2 = 0; j_2 < sample_size; ++j_2) {
                  for (int j_3 = 0; j_3 < sample_size; ++j_3) {
                     for (int j_4 = 0; j_4 < sample_size; ++j_4) {
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
      return val/samples_.size();
   }
   
private:
   int num_sweeps_  = 1000;
   int num_samples_ = 1;
   RealType beta_   = 1;
   std::vector<RealType> magnetizations_;
   std::vector<RealType> energies_;
   std::vector<std::uint64_t> seed_list_;
   std::vector<std::vector<SpinType>> samples_;
      
   std::vector<RealType> CalculateMagnetizations() const {
      if (static_cast<int>(samples_.size()) != num_samples_) {
         throw std::runtime_error("The size of samples is not equal to the num_samples_");
      }
      std::vector<RealType> magnetizations(num_samples_);
#pragma omp parallel for schedule(guided)
      for (int i = 0; i < num_samples_; ++i) {
         magnetizations[i] = model.CalculateMagnetization(samples_[i]);
      }
      return magnetizations;
   }
   
   std::vector<RealType> CalculateEnergies() const {
      if (static_cast<int>(samples_.size()) != num_samples_) {
         throw std::runtime_error("The size of samples is not equal to the num_samples_");
      }
      std::vector<RealType> energies(num_samples_);
#pragma omp parallel for schedule(guided)
      for (int i = 0; i < num_samples_; ++i) {
         energies[i] = model.CalculateEnergy(samples_[i]);
      }
      return energies;
   }
   
   void Sample(const model::PolynomialIsing<RealType> &model_input) {
      if (num_samples_ < 0) {
         throw std::runtime_error("num_samples must be non-negative integer.");
      }
      if (num_sweeps_ < 0) {
         throw std::runtime_error("num_sweeps must be non-negative integer.");
      }
      ClearResults();
      samples_.resize(num_samples_);
      
#pragma omp parallel for schedule(guided)
      for (int sample_count = 0; sample_count < num_samples_; sample_count++) {
         const auto system_size = model_input.GetSystemSize();
         std::mt19937_64 mt(seed_list_[sample_count]);
         std::uniform_real_distribution<RealType> dist_real(0, 1);
         std::uniform_int_distribution<int> dist_system_size(0, system_size - 1);
         
         // Set initial spin configurations
         auto &spin_configuration = samples_[sample_count];
         RandomizeSpins(&spin_configuration , &mt, system_size);
         
         // Prepare system size index list (0, 1, 2, ..., N - 1)
         std::vector<int> system_size_index_list(system_size);
         std::iota(system_size_index_list.begin(), system_size_index_list.end(), 0);
         
         // Set energy difference
         std::vector<RealType> energy_difference(system_size);
         ResetEnergyDifference(&energy_difference, &system_size_index_list, spin_configuration, model);
         
         for (int sweep_count = 0; sweep_count < num_sweeps_; sweep_count++) {
            for (int i = 0; i < system_size; i++) {
               const int index = dist_system_size(mt);
               if (energy_difference[index] <= 0.0 || std::exp(-beta_*energy_difference[index]) > dist_real(mt)) {
                  UpdateConfiguration(&spin_configuration, &energy_difference, &system_size_index_list, index, model);
               }
            }
         }
      }
   }
   
   void RandomizeSpins(std::vector<SpinType> *sample, std::mt19937_64 *random_number_engine, const int system_size) const {
      std::uniform_int_distribution<SpinType> dist(0, 1);
      sample->resize(system_size);
      for (int i = 0; i < system_size; i++) {
         (*sample)[i] = 2*dist(*random_number_engine) - 1;
      }
   }
   
};


} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
