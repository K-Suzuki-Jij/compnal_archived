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
//  Created by Kohei Suzuki on 2021/12/09.
//

#ifndef COMPNAL_BLAS_ORTHONORMALIZE_HPP_
#define COMPNAL_BLAS_ORTHONORMALIZE_HPP_

#include "../type/braket_vector.hpp"
#include <sstream>

namespace compnal {
namespace blas {

template<typename RealType>
void Orthonormalize(std::vector<type::BraketVector<RealType>> *vectors) {
   for (std::size_t i = 0; i < vectors->size(); ++i) {
      const std::int64_t dim = static_cast<std::int64_t>(vectors[i]->size());
      for (std::size_t j = 0; j < i; ++j) {
         const RealType inner_product = -1.0*CalculateInnerProduct(*vectors[i], (*vectors)[j]);
#pragma omp parallel for
         for (std::int64_t k = 0; k < dim; ++k) {
            vectors[i]->val[k] += inner_product*vectors[j]->val[k];
         }
      }
   }
}

template<typename RealType>
void Orthonormalize(type::BraketVector<RealType> *target_vector, const std::vector<type::BraketVector<RealType>> &vectors) {
   const std::int64_t dim = static_cast<std::int64_t>(target_vector->val.size());
   for (std::size_t i = 0; i < vectors.size(); ++i) {
      const RealType inner_product = -1.0*CalculateInnerProduct(*target_vector, vectors[i]);
#pragma omp parallel for
      for (std::int64_t j = 0; j < dim; ++j) {
         target_vector->val[j] += inner_product*vectors[i].val[j];
      }
   }
}

}  // namespace sparse_matrix
}  // namespace compnal


#endif /* COMPNAL_BLAS_ORTHONORMALIZE_HPP_ */
