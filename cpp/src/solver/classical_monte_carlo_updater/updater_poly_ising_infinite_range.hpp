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
//  Created by Kohei Suzuki on 2022/07/17.
//

#ifndef COMPNAL_SOLVER_UPDATER_POLY_ISING_INFINITE_RANGE_HPP_
#define COMPNAL_SOLVER_UPDATER_POLY_ISING_INFINITE_RANGE_HPP_

#include "../../model/polynomial_ising.hpp"
#include <vector>

namespace compnal {
namespace solver {

template<typename SpinType, typename RealType>
void UpdaterInfiniteRange(std::vector<SpinType> *sample,
                          std::vector<RealType> *energy_difference,
                          std::vector<int> *system_size_index_list,
                          const int index,
                          const model::PolynomialIsing<RealType> &model_input) {
   
   const int num_roop = model_input.GetSystemSize() - 1;
   const auto &interaction = model_input.GetInteraction();
   const SpinType target_spin = (*sample)[index];
   
   std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
   system_size_index_list->pop_back();
   
   if (interaction.size() == 2) {
      const RealType target_ineraction_deg2 = interaction[1];
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         (*energy_difference)[index_i2] += 4*target_ineraction_deg2*target_spin*(*sample)[index_i2];
      }
   }
   else if (interaction.size() == 3) {
      const RealType target_ineraction_deg2 = interaction[1];
      const RealType target_ineraction_deg3 = interaction[2];
      const SpinType spin_prod1 = target_spin;
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         const SpinType spin_prod2 = spin_prod1*(*sample)[index_i2];
         const RealType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         RealType sum_val_i2 = 0;
         for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
            const int index_i3 = (*system_size_index_list)[i3];
            const RealType diff_val_i3 = 4*target_ineraction_deg3*spin_prod2*(*sample)[index_i3];
            (*energy_difference)[index_i3] += diff_val_i3;
            sum_val_i2 += diff_val_i3;
         }
         (*energy_difference)[index_i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (interaction.size() == 4) {
      const RealType target_ineraction_deg2 = interaction[1];
      const RealType target_ineraction_deg3 = interaction[2];
      const RealType target_ineraction_deg4 = interaction[3];
      const SpinType spin_prod1 = target_spin;
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         const SpinType spin_prod2 = spin_prod1*(*sample)[index_i2];
         const RealType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         RealType sum_val_i2 = 0;
         for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
            const int index_i3 = (*system_size_index_list)[i3];
            const SpinType spin_prod3 = spin_prod2*(*sample)[index_i3];
            const RealType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            RealType sum_val_i3 = 0;
            for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
               const int index_i4 = (*system_size_index_list)[i4];
               const RealType diff_val_i4 = 4*target_ineraction_deg4*spin_prod3*(*sample)[index_i4];
               (*energy_difference)[index_i4] += diff_val_i4;
               sum_val_i3 += diff_val_i4;
            }
            (*energy_difference)[index_i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[index_i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (interaction.size() == 5) {
      const RealType target_ineraction_deg2 = interaction[1];
      const RealType target_ineraction_deg3 = interaction[2];
      const RealType target_ineraction_deg4 = interaction[3];
      const RealType target_ineraction_deg5 = interaction[4];
      const SpinType spin_prod1 = target_spin;
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         const SpinType spin_prod2 = spin_prod1*(*sample)[index_i2];
         const RealType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         RealType sum_val_i2 = 0;
         for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
            const int index_i3 = (*system_size_index_list)[i3];
            const SpinType spin_prod3 = spin_prod2*(*sample)[index_i3];
            const RealType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            RealType sum_val_i3 = 0;
            for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
               const int index_i4 = (*system_size_index_list)[i4];
               const SpinType spin_prod4 = spin_prod3*(*sample)[index_i4];
               const RealType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               RealType sum_val_i4 = 0;
               for (int i5 = i4 + 1; i5 < num_roop; ++i5) {
                  const int index_i5 = (*system_size_index_list)[i5];
                  const RealType diff_val_i5 = 4*target_ineraction_deg5*spin_prod4*(*sample)[index_i5];
                  (*energy_difference)[index_i5] += diff_val_i5;
                  sum_val_i4 += diff_val_i5;
               }
               (*energy_difference)[index_i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[index_i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[index_i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (interaction.size() == 6) {
      const RealType target_ineraction_deg2 = interaction[1];
      const RealType target_ineraction_deg3 = interaction[2];
      const RealType target_ineraction_deg4 = interaction[3];
      const RealType target_ineraction_deg5 = interaction[4];
      const RealType target_ineraction_deg6 = interaction[5];
      const SpinType spin_prod1 = target_spin;
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         const SpinType spin_prod2 = spin_prod1*(*sample)[index_i2];
         const RealType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         RealType sum_val_i2 = 0;
         for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
            const int index_i3 = (*system_size_index_list)[i3];
            const SpinType spin_prod3 = spin_prod2*(*sample)[index_i3];
            const RealType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            RealType sum_val_i3 = 0;
            for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
               const int index_i4 = (*system_size_index_list)[i4];
               const SpinType spin_prod4 = spin_prod3*(*sample)[index_i4];
               const RealType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               RealType sum_val_i4 = 0;
               for (int i5 = i4 + 1; i5 < num_roop; ++i5) {
                  const int index_i5 = (*system_size_index_list)[i5];
                  const SpinType spin_prod5 = spin_prod4*(*sample)[index_i5];
                  const RealType diff_val_i5 = 4*target_ineraction_deg5*spin_prod5;
                  RealType sum_val_i5 = 0;
                  for (int i6 = i5 + 1; i6 < num_roop; ++i6) {
                     const int index_i6 = (*system_size_index_list)[i6];
                     const RealType diff_val_i6 = 4*target_ineraction_deg6*spin_prod5*(*sample)[index_i6];
                     (*energy_difference)[index_i6] += diff_val_i6;
                     sum_val_i5 += diff_val_i6;
                  }
                  (*energy_difference)[index_i5] += sum_val_i5 + diff_val_i5;
                  sum_val_i4 += sum_val_i5 + diff_val_i5;
               }
               (*energy_difference)[index_i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[index_i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[index_i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (interaction.size() == 7) {
      const RealType target_ineraction_deg2 = interaction[1];
      const RealType target_ineraction_deg3 = interaction[2];
      const RealType target_ineraction_deg4 = interaction[3];
      const RealType target_ineraction_deg5 = interaction[4];
      const RealType target_ineraction_deg6 = interaction[5];
      const RealType target_ineraction_deg7 = interaction[6];
      const SpinType spin_prod1 = target_spin;
      for (int i2 = 0; i2 < num_roop; ++i2) {
         const int index_i2 = (*system_size_index_list)[i2];
         const SpinType spin_prod2 = spin_prod1*(*sample)[index_i2];
         const RealType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         RealType sum_val_i2 = 0;
         for (int i3 = i2 + 1; i3 < num_roop; ++i3) {
            const int index_i3 = (*system_size_index_list)[i3];
            const SpinType spin_prod3 = spin_prod2*(*sample)[index_i3];
            const RealType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            RealType sum_val_i3 = 0;
            for (int i4 = i3 + 1; i4 < num_roop; ++i4) {
               const int index_i4 = (*system_size_index_list)[i4];
               const SpinType spin_prod4 = spin_prod3*(*sample)[index_i4];
               const RealType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               RealType sum_val_i4 = 0;
               for (int i5 = i4 + 1; i5 < num_roop; ++i5) {
                  const int index_i5 = (*system_size_index_list)[i5];
                  const SpinType spin_prod5 = spin_prod4*(*sample)[index_i5];
                  const RealType diff_val_i5 = 4*target_ineraction_deg5*spin_prod5;
                  RealType sum_val_i5 = 0;
                  for (int i6 = i5 + 1; i6 < num_roop; ++i6) {
                     const int index_i6 = (*system_size_index_list)[i6];
                     const SpinType spin_prod6 = spin_prod5*(*sample)[index_i6];
                     const RealType diff_val_i6 = 4*target_ineraction_deg6*spin_prod6;
                     RealType sum_val_i6 = 0;
                     for (int i7 = i6 + 1; i7 < num_roop; ++i7) {
                        const int index_i7 = (*system_size_index_list)[i7];
                        const RealType diff_val_i7 = 4*target_ineraction_deg7*spin_prod6*(*sample)[index_i7];
                        (*energy_difference)[index_i7] += diff_val_i7;
                        sum_val_i6 += diff_val_i7;
                     }
                     (*energy_difference)[index_i6] += sum_val_i6 + diff_val_i6;
                     sum_val_i5 += sum_val_i6 + diff_val_i6;
                  }
                  (*energy_difference)[index_i5] += sum_val_i5 + diff_val_i5;
                  sum_val_i4 += sum_val_i5 + diff_val_i5;
               }
               (*energy_difference)[index_i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[index_i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[index_i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else {
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
   }
   (*energy_difference)[index] *= -1;
   (*sample)[index] *= -1;
   system_size_index_list->push_back(index);
   std::swap((*system_size_index_list)[index], (*system_size_index_list).back());
}


}
}

#endif /* COMPNAL_SOLVER_UPDATER_POLY_ISING_INFINITE_RANGE_HPP_ */
