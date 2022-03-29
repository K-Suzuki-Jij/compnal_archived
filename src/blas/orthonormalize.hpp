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

#include "braket_vector.hpp"
#include <sstream>
#include <limits>

namespace compnal {
namespace blas {

//! @brief Orthonormalize BraketVector list by Gram–Schmidt orthonormalization.
//! @tparam ElementType The value type of BraketVector.
//! @param vectors The pointer of BraketVector list.
template<typename ElementType>
void Orthonormalize(std::vector<BraketVector<ElementType>> *vectors, const bool flag_normalize = true) {
   std::int64_t vec_size = static_cast<std::int64_t>(vectors->size());
   for (std::int64_t i = 0; i < vec_size; ++i) {
      const std::int64_t dim = static_cast<std::int64_t>((*vectors)[i].Size());
      if (dim < vec_size) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Cannot orthonormalize all the vectors" << std::endl;
         throw std::runtime_error(ss.str());
      }
      for (std::int64_t j = 0; j < i; ++j) {
         if ((*vectors)[i].value_list.size() != (*vectors)[j].value_list.size()) {
            std::stringstream ss;
            ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
            ss << "Dimenstion of BraketVector list doen not match each other" << std::endl;
            ss << "Failed to orthonormalize." << std::endl;
            throw std::runtime_error(ss.str());
         }
         const ElementType inner_product = (*vectors)[i]*(*vectors)[j];
#pragma omp parallel for
         for (std::int64_t k = 0; k < dim; ++k) {
            (*vectors)[i].value_list[k] -= inner_product*(*vectors)[j].value_list[k];
         }
      }
      //Normalize
      if (flag_normalize) {
         const auto norm = (*vectors)[i].CalculateL2Norm();
         if (std::abs(norm) < 4*std::numeric_limits<ElementType>::epsilon()) {
            std::stringstream ss;
            ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
            ss << "Detected vectors pointing to the same direction." << std::endl;
            ss << "Failed to orthonormalize." << std::endl;
            throw std::runtime_error(ss.str());
         }
         (*vectors)[i].MultiplyByScalar(ElementType{1.0}/norm);
      }
   }
}

//! @brief Orthonormalize a BraketVector to BraketVector list by Gram–Schmidt orthonormalization.
//! @tparam ElementType The value type of BraketVector.
//! @param target_vector The pointer of the BraketVector to be orthonormalized.
//! @param vectors The BraketVector list. Note that this vectors must be orthonormalized.
template<typename ElementType>
void Orthonormalize(BraketVector<ElementType> *target_vector, const std::vector<BraketVector<ElementType>> &vectors, const bool flag_normalize = true) {
   const std::int64_t dim = static_cast<std::int64_t>(target_vector->value_list.size());
   for (std::size_t i = 0; i < vectors.size(); ++i) {
      if ((*target_vector).value_list.size() != vectors[i].value_list.size()) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Dimenstion of BraketVector list doen not match each other" << std::endl;
         ss << "Failed to orthonormalize." << std::endl;
         throw std::runtime_error(ss.str());
      }
      const ElementType inner_product = ElementType{1.0}*(*target_vector)*(vectors[i]);
#pragma omp parallel for
      for (std::int64_t j = 0; j < dim; ++j) {
         target_vector->value_list[j] -= inner_product*vectors[i].value_list[j];
      }
   }
   //Normalize
   if (flag_normalize) {
      const auto norm = (*target_vector).CalculateL2Norm();
      if (std::abs(norm) < 4*std::numeric_limits<ElementType>::epsilon()) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Detected vectors pointing to the same direction." << std::endl;
         ss << "Failed to orthonormalize." << std::endl;
         throw std::runtime_error(ss.str());
      }
      (*target_vector).MultiplyByScalar(ElementType{1.0}/norm);
   }
}

}  // namespace sparse_matrix
}  // namespace compnal


#endif /* COMPNAL_BLAS_ORTHONORMALIZE_HPP_ */
