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
#include <vector>
#include <random>

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
   model(model_input), updater(updater_input) {
      SetRandomSeed();
   }
   
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
   
   void SetInverseTemperature(const RealType beta) {
      beta_ = beta;
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
   
   template<typename T>
   const std::vector<SpinType> &GetSample(const T index) const {
      return samples_.at(index);
   }
   
   void ClearSamples() {
      samples_.clear();
      samples_.shrink_to_fit();
   }
   
   void Run() {
      Sample(model);
   }
   
   
private:
   int num_sweeps_  = 1000;
   int num_samples_ = 1;
   RealType beta_   = 100;
   std::vector<int> seed_list_;
   std::vector<std::vector<SpinType>> samples_;
   
   void SetRandomSeed() {
      seed_list_.resize(num_samples_);
      std::random_device rnd;
      for (std::size_t i = 0; i < seed_list_.size(); ++i) {
         seed_list_[i] = rnd();
      }
   }
   
   void Sample(const model::PolynomialIsing<RealType> &model_input) {
      if (num_samples_ < 0) {
         throw std::runtime_error("num_samples must be non-negative integer.");
      }
      if (num_sweeps_ < 0) {
         throw std::runtime_error("num_sweeps must be non-negative integer.");
      }
      ClearSamples();
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
         
         // Set energy difference
         std::vector<RealType> energy_difference(system_size);
         ResetEnergyDifference(&energy_difference, spin_configuration, model);
         
         // Prepare system size index list (0, 1, 2, ..., N - 1)
         std::vector<int> system_size_index_list(system_size);
         std::iota(system_size_index_list.begin(), system_size_index_list.end(), 0);
         
         for (int sweep_count = 0; sweep_count < num_sweeps_; sweep_count++) {
            for (int i = 0; i < system_size; i++) {
               const int index = dist_system_size(mt);
               if (energy_difference[index] <= 0.0 || std::exp(-beta_*energy_difference[index]) > dist_real(mt)) {
                  UpdateSpinAndEnergyDifference(&spin_configuration, &energy_difference, &system_size_index_list, index, model);
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
   
   void ResetEnergyDifference(std::vector<RealType> *energy_difference,
                              const std::vector<SpinType> &sample,
                              const model::PolynomialIsing<RealType> &model_input) const {
      
      if (static_cast<int>(sample.size()) != model_input.GetSystemSize()) {
         throw std::runtime_error("The sample size is not equal to the system size.");
      }
      
      const int system_size = model_input.GetSystemSize();
      const auto &interaction = model_input.GetInteraction();
      energy_difference->clear();
      energy_difference->resize(system_size);
      
      if (model_input.GetLattice() == model::Lattice::CHAIN) {
         throw std::runtime_error("CHAIN is under construction.");
      }
      else if (model_input.GetLattice() == model::Lattice::SQUARE) {
         throw std::runtime_error("SQUARE is under construction.");
      }
      else if (model_input.GetLattice() == model::Lattice::TRIANGLE) {
         throw std::runtime_error("TRIANGLE is under construction.");
      }
      else if (model_input.GetLattice() == model::Lattice::HONEYCOMB) {
         throw std::runtime_error("HONEYCOMB is under construction.");
      }
      else if (model_input.GetLattice() == model::Lattice::CUBIC) {
         throw std::runtime_error("CUBIC is under construction.");
      }
      else if (model_input.GetLattice() == model::Lattice::INFINIT_RANGE) {
         std::vector<int> system_size_index_list(system_size);
         std::iota(system_size_index_list.begin(), system_size_index_list.end(), 0);
         for (int index = 0; index < system_size; ++index) {
            RealType val = 0.0;
            const SpinType target_spin = sample[index];
            
            // Erase index in seed.
            std::swap(system_size_index_list[index], system_size_index_list.back());
            
            if (interaction.size() >= 1) {
               val += interaction[0]*sample[index];
            }
            for (int p = 1; p < static_cast<int>(interaction.size()); ++p) {
               if (std::abs(interaction[p]) > std::numeric_limits<RealType>::epsilon()) {
                  const RealType target_ineraction = interaction[p];
                  
                  std::vector<int> indices(p);
                  int start_index = 0;
                  int size = 0;
                  
                  while (true) {
                     for (int i = start_index; i < system_size - 1; ++i) {
                        indices[size++] = i;
                        if (size == p) {
                           for (int j = 0; j < p; ++j) {
                              val += target_ineraction*sample[system_size_index_list[indices[j]]]*target_spin;
                           }
                           break;
                        }
                     }
                     --size;
                     if (size < 0) {
                        break;
                     }
                     start_index = indices[size] + 1;
                  }
                  
               }
            }
            (*energy_difference)[index] = -2.0*val;
            
            // Restore index to seed.
            std::swap(system_size_index_list[index], system_size_index_list.back());
         }
      }
      else if (model_input.GetLattice() == model::Lattice::ANY_TYPE) {
         throw std::runtime_error("ANY_TYPE is under construction.");
      }
      else {
         throw std::runtime_error("Unknown lattice type detected.");
      }
   }
   
   void UpdateSpinAndEnergyDifference(std::vector<SpinType> *sample,
                                      std::vector<RealType> *energy_difference,
                                      std::vector<int> *system_size_index_list,
                                      const int index,
                                      const model::PolynomialIsing<RealType> &model_input) const {
      
      const int num_roop = model_input.GetSystemSize() - 1;
      const auto &interaction = model_input.GetInteraction();
      const SpinType target_spin = (*sample)[index];
      
      if (interaction.size() == 1) {
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
      }
      else if (interaction.size() == 2) {
         const RealType target_ineraction_deg2 = interaction[1];
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         system_size_index_list->pop_back();
         
         for (int i2 = 0; i2 < num_roop; ++i2) {
            const int index_i2 = (*system_size_index_list)[i2];
            (*energy_difference)[index_i2] += 4*target_ineraction_deg2*target_spin*(*sample)[index_i2];
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         system_size_index_list->push_back(index);
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
      else if (interaction.size() == 3) {
         const RealType target_ineraction_deg2 = interaction[1];
         const RealType target_ineraction_deg3 = interaction[2];
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         system_size_index_list->pop_back();
         
         const RealType temp_val_deg2 = 4*target_ineraction_deg2*target_spin;
         for (int i2 = 0; i2 < num_roop; ++i2) {
            const int index_i2 = (*system_size_index_list)[i2];
            const SpinType spin_i2 = (*sample)[index_i2];
            const RealType temp_val_deg3 = 4*target_ineraction_deg3*target_spin*spin_i2;
            (*energy_difference)[index_i2] += temp_val_deg2*spin_i2;
            for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
               const int index_i3 = (*system_size_index_list)[i3];
               (*energy_difference)[index_i3] += temp_val_deg3*(*sample)[index_i3];
            }
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         system_size_index_list->push_back(index);
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
      else if (interaction.size() == 4) {
         const RealType target_ineraction_deg2 = interaction[1];
         const RealType target_ineraction_deg3 = interaction[2];
         const RealType target_ineraction_deg4 = interaction[3];
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         system_size_index_list->pop_back();
         
         const RealType temp_val_deg2 = 4*target_ineraction_deg2*target_spin;
         for (int i2 = 0; i2 < num_roop; ++i2) {
            const int index_i2 = (*system_size_index_list)[i2];
            const SpinType spin_i2 = (*sample)[index_i2];
            const RealType temp_val_deg3 = 4*target_ineraction_deg3*target_spin*spin_i2;
            (*energy_difference)[index_i2] += temp_val_deg2*spin_i2;
            for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
               const int index_i3 = (*system_size_index_list)[i3];
               const SpinType spin_i3 = (*sample)[index_i3];
               const RealType temp_val_deg4 = 4*target_ineraction_deg4*target_spin*spin_i2*spin_i3;
               (*energy_difference)[index_i3] += temp_val_deg3*spin_i3;
               for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
                  const int index_i4 = (*system_size_index_list)[i4];
                  (*energy_difference)[index_i4] += temp_val_deg4*(*sample)[index_i4];
               }
            }
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         system_size_index_list->push_back(index);
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
      else if (interaction.size() == 5) {
         const RealType target_ineraction_deg2 = interaction[1];
         const RealType target_ineraction_deg3 = interaction[2];
         const RealType target_ineraction_deg4 = interaction[3];
         const RealType target_ineraction_deg5 = interaction[4];
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         system_size_index_list->pop_back();
         
         const RealType temp_val_deg2 = 4*target_ineraction_deg2*target_spin;
         for (int i2 = 0; i2 < num_roop; ++i2) {
            const int index_i2 = (*system_size_index_list)[i2];
            const SpinType spin_i2 = (*sample)[index_i2];
            const RealType temp_val_deg3 = 4*target_ineraction_deg3*target_spin*spin_i2;
            (*energy_difference)[index_i2] += temp_val_deg2*spin_i2;
            for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
               const int index_i3 = (*system_size_index_list)[i3];
               const SpinType spin_i3 = (*sample)[index_i3];
               const RealType temp_val_deg4 = 4*target_ineraction_deg4*target_spin*spin_i2*spin_i3;
               (*energy_difference)[index_i3] += temp_val_deg3*spin_i3;
               for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
                  const int index_i4 = (*system_size_index_list)[i4];
                  const SpinType spin_i4 = (*sample)[index_i4];
                  const RealType temp_val_deg5 = 4*target_ineraction_deg5*target_spin*spin_i2*spin_i3*spin_i4;
                  (*energy_difference)[index_i4] += temp_val_deg4*spin_i4;
                  for (int i5 = i4 + 1; i5 < num_roop; ++i5) {
                     const int index_i5 = (*system_size_index_list)[i5];
                     (*energy_difference)[index_i5] += temp_val_deg5*(*sample)[index_i5];
                  }
               }
            }
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         system_size_index_list->push_back(index);
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
      else if (interaction.size() == 6) {
         const RealType target_ineraction_deg2 = interaction[1];
         const RealType target_ineraction_deg3 = interaction[2];
         const RealType target_ineraction_deg4 = interaction[3];
         const RealType target_ineraction_deg5 = interaction[4];
         const RealType target_ineraction_deg6 = interaction[5];
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         system_size_index_list->pop_back();
         
         const RealType temp_val_deg2 = 4*target_ineraction_deg2*target_spin;
         for (int i2 = 0; i2 < num_roop; ++i2) {
            const int index_i2 = (*system_size_index_list)[i2];
            const SpinType spin_i2 = (*sample)[index_i2];
            const RealType temp_val_deg3 = 4*target_ineraction_deg3*target_spin*spin_i2;
            (*energy_difference)[index_i2] += temp_val_deg2*spin_i2;
            for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
               const int index_i3 = (*system_size_index_list)[i3];
               const SpinType spin_i3 = (*sample)[index_i3];
               const RealType temp_val_deg4 = 4*target_ineraction_deg4*target_spin*spin_i2*spin_i3;
               (*energy_difference)[index_i3] += temp_val_deg3*spin_i3;
               for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
                  const int index_i4 = (*system_size_index_list)[i4];
                  const SpinType spin_i4 = (*sample)[index_i4];
                  const RealType temp_val_deg5 = 4*target_ineraction_deg5*target_spin*spin_i2*spin_i3*spin_i4;
                  (*energy_difference)[index_i4] += temp_val_deg4*spin_i4;
                  for (int i5 = i4 + 1; i5 < num_roop; ++i5) {
                     const int index_i5 = (*system_size_index_list)[i5];
                     const SpinType spin_i5 = (*sample)[index_i5];
                     const RealType temp_val_deg6 = 4*target_ineraction_deg6*target_spin*spin_i2*spin_i3*spin_i4*spin_i5;
                     (*energy_difference)[index_i5] += temp_val_deg5*spin_i5;
                     for (int i6 = i5 + 1; i6 < num_roop; ++i6) {
                        const int index_i6 = (*system_size_index_list)[i6];
                        (*energy_difference)[index_i6] += temp_val_deg6*(*sample)[index_i6];
                     }
                  }
               }
            }
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         system_size_index_list->push_back(index);
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
      else {
         // Erase index in system_index_list_.
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
         
         for (int p = 1; p < static_cast<int>(interaction.size()); ++p) {
            if (std::abs(interaction[p]) > std::numeric_limits<RealType>::epsilon()) {
               const RealType target_ineraction = interaction[p];
               
               std::vector<int> indices(p);
               int start_index = 0;
               int size = 0;
               
               while (true) {
                  for (int i = start_index; i < num_roop; ++i) {
                     indices[size++] = i;
                     if (size == p) {
                        SpinType sign = 1;
                        for (int j = 0; j < p; ++j) {
                           sign *= (*sample)[(*system_size_index_list)[indices[j]]];
                        }
                        const RealType val = 4*target_spin*sign*target_ineraction;
                        for (int j = 0; j < p; ++j) {
                           (*energy_difference)[(*system_size_index_list)[indices[j]]] += val;
                        }
                        break;
                     }
                  }
                  --size;
                  if (size < 0) {
                     break;
                  }
                  start_index = indices[size] + 1;
               }
            }
         }
         
         (*energy_difference)[index] *= -1;
         (*sample)[index] *= -1;
         std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
      }
   }
   
   
   
};


} // namespace solver
} // namespace compnal


#endif /* COMPNAL_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
