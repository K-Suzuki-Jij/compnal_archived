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
//  ssf_polynomial_ising_infinite_range.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_INFINITE_RANGE_HPP_
#define COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_INFINITE_RANGE_HPP_

#include "../../lattice/infinite_range.hpp"
#include "../../model/classical/polynomial_ising.hpp"
#include "../../utility/type.hpp"

namespace compnal {
namespace solver {
namespace updater {

template<typename RealType>
void SetEnergyDifference(std::vector<typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::ValueType> *energy_difference,
                         const std::vector<typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::OPType> &sample,
                         const model::PolynomialIsing<lattice::InfiniteRange, RealType> &model) {
   
   if (static_cast<std::int32_t>(energy_difference->size()) != model.GetSystemSize()) {
      throw std::runtime_error("The size of energy_difference is not equal to the system size.");
   }
   if (static_cast<std::int32_t>(sample.size()) != model.GetSystemSize()) {
      throw std::runtime_error("The sample size is not equal to the system size.");
   }
   
   using OPType = typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::ValueType;
   const std::int32_t system_size = model.GetSystemSize();
   const std::vector<ValueType> &interaction = model.GetInteraction();

   for (std::int32_t index = 0; index < system_size; ++index) {
      ValueType val = 0.0;
      const OPType target_spin = sample[index];
   
      if (interaction.size() >= 2) {
         val += interaction[1]*target_spin;
      }
      for (std::int32_t p = 2; p < static_cast<std::int32_t>(interaction.size()); ++p) {
         if (std::abs(interaction[p]) > std::numeric_limits<ValueType>::epsilon()) {
            const ValueType target_ineraction = interaction[p];
            std::vector<std::int32_t> indices(p - 1);
            std::int32_t start_index = 0;
            std::int32_t size = 0;
   
            while (true) {
               for (std::int32_t i = start_index; i < system_size - 1; ++i) {
                  indices[size++] = i;
                  if (size == p - 1) {
                     OPType sign = 1;
                     for (std::int32_t j = 0; j < p - 1; ++j) {
                        if (indices[j] >= index) {
                           sign *= sample[indices[j] + 1];
                        }
                        else {
                           sign *= sample[indices[j]];
                        }
                     }
                     val += target_ineraction*sign*target_spin;
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
   }
   
}


template<typename RealType>
void UpdateConfiguration(std::vector<typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::OPType> *sample,
                         std::vector<typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::ValueType> *energy_difference,
                         const std::int32_t index,
                         const model::PolynomialIsing<lattice::InfiniteRange, RealType> &model) {
   
   using OPType = typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::OPType;
   using ValueType = typename model::PolynomialIsing<lattice::InfiniteRange, RealType>::ValueType;
   const std::vector<ValueType> &interaction = model.GetInteraction();
   const std::int32_t degree = model.GetDegree();
   const OPType target_spin = (*sample)[index];
   const std::int32_t system_size = model.GetSystemSize();
   
   if (degree == 2) {
      const ValueType target_ineraction_deg2 = interaction[2];
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         (*energy_difference)[i2] += 4*target_ineraction_deg2*target_spin*(*sample)[i2];
      }
   }
   else if (degree == 3) {
      const ValueType target_ineraction_deg2 = interaction[2];
      const ValueType target_ineraction_deg3 = interaction[3];
      const OPType spin_prod1 = target_spin;
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         const OPType spin_prod2 = spin_prod1*(*sample)[i2];
         const ValueType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         ValueType sum_val_i2 = 0;
         for (std::int32_t i3 = i2 + 1; i3 < index; ++i3) {
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod2*(*sample)[i3];
            (*energy_difference)[i3] += diff_val_i3;
            sum_val_i2 += diff_val_i3;
         }
         for (std::int32_t i3 = std::max(index + 1, i2 + 1); i3 < system_size; ++i3) {
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod2*(*sample)[i3];
            (*energy_difference)[i3] += diff_val_i3;
            sum_val_i2 += diff_val_i3;
         }
         (*energy_difference)[i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (degree == 4) {
      const ValueType target_ineraction_deg2 = interaction[2];
      const ValueType target_ineraction_deg3 = interaction[3];
      const ValueType target_ineraction_deg4 = interaction[4];
      const OPType spin_prod1 = target_spin;
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         const OPType spin_prod2 = spin_prod1*(*sample)[i2];
         const ValueType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         ValueType sum_val_i2 = 0;
         for (std::int32_t i3 = i2 + 1; i3 < system_size; ++i3) {
            if (i3 == index) {continue;}
            const OPType spin_prod3 = spin_prod2*(*sample)[i3];
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            ValueType sum_val_i3 = 0;
            for (std::int32_t i4 = i3 + 1; i4 < index; ++i4) {
               const ValueType diff_val_i4 = 4*target_ineraction_deg4*spin_prod3*(*sample)[i4];
               (*energy_difference)[i4] += diff_val_i4;
               sum_val_i3 += diff_val_i4;
            }
            for (std::int32_t i4 = std::max(index + 1, i3 + 1); i4 < system_size; ++i4) {
               const ValueType diff_val_i4 = 4*target_ineraction_deg4*spin_prod3*(*sample)[i4];
               (*energy_difference)[i4] += diff_val_i4;
               sum_val_i3 += diff_val_i4;
            }
            (*energy_difference)[i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (degree == 5) {
      const ValueType target_ineraction_deg2 = interaction[2];
      const ValueType target_ineraction_deg3 = interaction[3];
      const ValueType target_ineraction_deg4 = interaction[4];
      const ValueType target_ineraction_deg5 = interaction[5];
      const OPType spin_prod1 = target_spin;
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         const OPType spin_prod2 = spin_prod1*(*sample)[i2];
         const ValueType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         ValueType sum_val_i2 = 0;
         for (std::int32_t i3 = i2 + 1; i3 < system_size; ++i3) {
            if (i3 == index) {continue;}
            const OPType spin_prod3 = spin_prod2*(*sample)[i3];
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            ValueType sum_val_i3 = 0;
            for (std::int32_t i4 = i3 + 1; i4 < system_size; ++i4) {
               if (i4 == index) {continue;}
               const OPType spin_prod4 = spin_prod3*(*sample)[i4];
               const ValueType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               ValueType sum_val_i4 = 0;
               for (std::int32_t i5 = i4 + 1; i5 < index; ++i5) {
                  const ValueType diff_val_i5 = 4*target_ineraction_deg5*spin_prod4*(*sample)[i5];
                  (*energy_difference)[i5] += diff_val_i5;
                  sum_val_i4 += diff_val_i5;
               }
               for (std::int32_t i5 = std::max(index + 1, i4 + 1); i5 < system_size; ++i5) {
                  const ValueType diff_val_i5 = 4*target_ineraction_deg5*spin_prod4*(*sample)[i5];
                  (*energy_difference)[i5] += diff_val_i5;
                  sum_val_i4 += diff_val_i5;
               }
               (*energy_difference)[i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (degree == 6) {
      const ValueType target_ineraction_deg2 = interaction[2];
      const ValueType target_ineraction_deg3 = interaction[3];
      const ValueType target_ineraction_deg4 = interaction[4];
      const ValueType target_ineraction_deg5 = interaction[5];
      const ValueType target_ineraction_deg6 = interaction[6];
      const OPType spin_prod1 = target_spin;
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         const OPType spin_prod2 = spin_prod1*(*sample)[i2];
         const ValueType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         ValueType sum_val_i2 = 0;
         for (std::int32_t i3 = i2 + 1; i3 < system_size; ++i3) {
            if (i3 == index) {continue;}
            const OPType spin_prod3 = spin_prod2*(*sample)[i3];
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            ValueType sum_val_i3 = 0;
            for (std::int32_t i4 = i3 + 1; i4 < system_size; ++i4) {
               if (i4 == index) {continue;}
               const OPType spin_prod4 = spin_prod3*(*sample)[i4];
               const ValueType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               ValueType sum_val_i4 = 0;
               for (std::int32_t i5 = i4 + 1; i5 < system_size; ++i5) {
                  if (i5 == index) {continue;}
                  const OPType spin_prod5 = spin_prod4*(*sample)[i5];
                  const ValueType diff_val_i5 = 4*target_ineraction_deg5*spin_prod5;
                  ValueType sum_val_i5 = 0;
                  for (std::int32_t i6 = i5 + 1; i6 < index; ++i6) {
                     const ValueType diff_val_i6 = 4*target_ineraction_deg6*spin_prod5*(*sample)[i6];
                     (*energy_difference)[i6] += diff_val_i6;
                     sum_val_i5 += diff_val_i6;
                  }
                  for (std::int32_t i6 = std::max(index + 1, i5 + 1); i6 < system_size; ++i6) {
                     const ValueType diff_val_i6 = 4*target_ineraction_deg6*spin_prod5*(*sample)[i6];
                     (*energy_difference)[i6] += diff_val_i6;
                     sum_val_i5 += diff_val_i6;
                  }
                  (*energy_difference)[i5] += sum_val_i5 + diff_val_i5;
                  sum_val_i4 += sum_val_i5 + diff_val_i5;
               }
               (*energy_difference)[i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else if (degree == 7) {
      const ValueType target_ineraction_deg2 = interaction[2];
      const ValueType target_ineraction_deg3 = interaction[3];
      const ValueType target_ineraction_deg4 = interaction[4];
      const ValueType target_ineraction_deg5 = interaction[5];
      const ValueType target_ineraction_deg6 = interaction[6];
      const ValueType target_ineraction_deg7 = interaction[7];
      const OPType spin_prod1 = target_spin;
      for (std::int32_t i2 = 0; i2 < system_size; ++i2) {
         if (i2 == index) {continue;}
         const OPType spin_prod2 = spin_prod1*(*sample)[i2];
         const ValueType diff_val_i2 = 4*target_ineraction_deg2*spin_prod2;
         ValueType sum_val_i2 = 0;
         for (std::int32_t i3 = i2 + 1; i3 < system_size; ++i3) {
            if (i3 == index) {continue;}
            const OPType spin_prod3 = spin_prod2*(*sample)[i3];
            const ValueType diff_val_i3 = 4*target_ineraction_deg3*spin_prod3;
            ValueType sum_val_i3 = 0;
            for (std::int32_t i4 = i3 + 1; i4 < system_size; ++i4) {
               if (i4 == index) {continue;}
               const OPType spin_prod4 = spin_prod3*(*sample)[i4];
               const ValueType diff_val_i4 = 4*target_ineraction_deg4*spin_prod4;
               ValueType sum_val_i4 = 0;
               for (std::int32_t i5 = i4 + 1; i5 < system_size; ++i5) {
                  if (i5 == index) {continue;}
                  const OPType spin_prod5 = spin_prod4*(*sample)[i5];
                  const ValueType diff_val_i5 = 4*target_ineraction_deg5*spin_prod5;
                  ValueType sum_val_i5 = 0;
                  for (std::int32_t i6 = i5 + 1; i6 < system_size; ++i6) {
                     if (i6 == index) {continue;}
                     const OPType spin_prod6 = spin_prod5*(*sample)[i6];
                     const ValueType diff_val_i6 = 4*target_ineraction_deg6*spin_prod6;
                     ValueType sum_val_i6 = 0;
                     for (std::int32_t i7 = i6 + 1; i7 < index; ++i7) {
                        const ValueType diff_val_i7 = 4*target_ineraction_deg7*spin_prod6*(*sample)[i7];
                        (*energy_difference)[i7] += diff_val_i7;
                        sum_val_i6 += diff_val_i7;
                     }
                     for (std::int32_t i7 = std::max(index + 1, i6 + 1); i7 < system_size; ++i7) {
                        const ValueType diff_val_i7 = 4*target_ineraction_deg7*spin_prod6*(*sample)[i7];
                        (*energy_difference)[i7] += diff_val_i7;
                        sum_val_i6 += diff_val_i7;
                     }
                     (*energy_difference)[i6] += sum_val_i6 + diff_val_i6;
                     sum_val_i5 += sum_val_i6 + diff_val_i6;
                  }
                  (*energy_difference)[i5] += sum_val_i5 + diff_val_i5;
                  sum_val_i4 += sum_val_i5 + diff_val_i5;
               }
               (*energy_difference)[i4] += sum_val_i4 + diff_val_i4;
               sum_val_i3 += sum_val_i4 + diff_val_i4;
            }
            (*energy_difference)[i3] += sum_val_i3 + diff_val_i3;
            sum_val_i2 += sum_val_i3 + diff_val_i3;
         }
         (*energy_difference)[i2] += sum_val_i2 + diff_val_i2;
      }
   }
   else {
      for (std::int32_t p = 2; p < static_cast<std::int32_t>(interaction.size()); ++p) {
         if (std::abs(interaction[p]) > std::numeric_limits<ValueType>::epsilon()) {
            const ValueType target_ineraction = interaction[p];
            std::vector<std::int32_t> indices(p);
            std::int32_t start_index = 0;
            std::int32_t size = 0;
            
            while (true) {
               for (std::int32_t i = start_index; i < system_size - 1; ++i) {
                  indices[size++] = i;
                  if (size == p) {
                     OPType sign = 1;
                     for (std::int32_t j = 0; j < p; ++j) {
                        if (indices[j] >= index) {
                           sign *= (*sample)[indices[j] + 1];
                        }
                        else {
                           sign *= (*sample)[indices[j]];
                        }
                     }
                     const ValueType val = 4*target_spin*sign*target_ineraction;
                     for (std::int32_t j = 0; j < p; ++j) {
                        if (indices[j] >= index) {
                           (*energy_difference)[indices[j] + 1] += val;
                        }
                        else {
                           (*energy_difference)[indices[j]] += val;
                        }
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
   
}


} // namespace updater
} // namespace solver
} // namespace compnal

#endif /* COMPNAL_SOLVER_UPDATER_SSF_POLYNOMIAL_ISING_INFINITE_RANGE_HPP_ */
