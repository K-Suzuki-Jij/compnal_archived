//
//  orthonormalization.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/09.
//

#ifndef COMPNAL_SPARSE_MATRIX_ORTHONORMALIZE_HPP_
#define COMPNAL_SPARSE_MATRIX_ORTHONORMALIZE_HPP_

#include "braket_vector.hpp"
#include <sstream>

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
void Orthonormalize(std::vector<BraketVector<RealType>> *vectors) {
   for (std::size_t i = 0; i < vectors->size(); ++i) {
      const int64_t dim = static_cast<std::int64_t>(vectors[i]->size());
      for (std::size_t j = 0; j < i; ++j) {
         const RealType inner_product = -1.0*CalculateInnerProduct(*vectors[i], (*vectors)[j]);
#pragma omp parallel for
         for (std::int64_t k = 0; k < dim; ++k) {
            vectors[i]->val[k] += inner_product*vectors[j]->val[k];
         }
      }
      vectors[i]->Normalize();
   }
}

template<typename RealType>
void Orthonormalize(BraketVector<RealType> *target_vector, const std::vector<BraketVector<RealType>> &vectors) {
   const int64_t dim = static_cast<std::int64_t>(target_vector->val.size());
   for (std::size_t i = 0; i < vectors.size(); ++i) {
      const RealType inner_product = -1.0*CalculateInnerProduct(*target_vector, vectors[i]);
#pragma omp parallel for
      for (std::int64_t j = 0; j < dim; ++j) {
         target_vector->val[j] += inner_product*vectors[i].val[j];
      }
   }
   target_vector->Normalize();
}

}  // namespace sparse_matrix
}  // namespace compnal


#endif /* COMPNAL_SPARSE_MATRIX_ORTHONORMALIZE_HPP_ */
