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
//  Created by Kohei Suzuki on 2021/11/05.
//

#ifndef COMPNAL_BLAS_MATRIX_VECTOR_OPERATION_HPP_
#define COMPNAL_BLAS_MATRIX_VECTOR_OPERATION_HPP_

#include "braket_vector.hpp"
#include "compressed_row_storage.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace blas {

//! @brief Calculate matrix vector product.
//! \f$ \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! @tparam ReturnType The value type of output vector.
//! @tparam T1 The value type of the coefficient \f$ c\f$.
//! @tparam T2 The value type of the matirx \f$ \hat{M}\f$.
//! @tparam T3 The value type of the vector \f$ \boldsymbol{v}\f$.
//! @param vector_out The result of matrix vector product \f$
//! \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! @param coeff The coefficient \f$ c\f$.
//! @param matrix_in The matrix \f$ \hat{M} \f$.
//! @param vector_in The vector \f$ \boldsymbol{v} \f$.
template <typename ReturnType, typename T1, typename T2, typename T3>
void CalculateMatrixVectorProduct(BraketVector<ReturnType> *vector_out, const T1 coeff, const CRS<T2> &matrix_in,
                                  const BraketVector<T3> &vector_in) {
   if (matrix_in.col_dim != vector_in.Size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim << std::endl;
      ss << "The dimension of the input vector is " << vector_in.Size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->value_list.resize(matrix_in.row_dim);
   using T2T3 = decltype(std::declval<T2>() * std::declval<T3>());

#pragma omp parallel for
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      T2T3 temp = 0;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1]; ++j) {
         temp += matrix_in.val[j] * vector_in.value_list[matrix_in.col[j]];
      }
      vector_out->value_list[i] = temp * coeff;
   }
}

//! @brief Calculate matrix vector product.
//! \f$ \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! Note that the matrix \f$ \hat{M}\f$ must be symmetric and their elements are
//! stored only in lower triangle elements. In addition, the diagonal elements
//! must also stored even if they are zero.
//! @tparam ReturnType The value type of output vector.
//! @tparam T1 The value type of the coefficient \f$ c\f$.
//! @tparam T2 The value type of the matirx \f$ \hat{M}\f$.
//! @tparam T3 The value type of the vector \f$ \boldsymbol{v}\f$.
//! @param vector_out The pointer of result for matrix vector product \f$
//! \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! @param coeff The coefficient \f$ c\f$.
//! @param matrix_in The matrix \f$ \hat{M} \f$.
//! @param vector_in The vector \f$ \boldsymbol{v} \f$.
//! @param vectors_work The pointer of a temporary two-dimensional array used
//! for calculations. Their elements must be zero. Note that this working array
//! is used only when openmp is active.
template <typename ReturnType, typename T1, typename T2, typename T3>
void CalculateSymmetricMatrixVectorProduct(
    BraketVector<ReturnType> *vector_out, const T1 coeff, const CRS<T2> &matrix_in, const BraketVector<T3> &vector_in,
    std::vector<std::vector<decltype(std::declval<T2>() * std::declval<T3>())>> *vectors_work = nullptr) {
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "The input matrix is not symmetric" << std::endl;
      throw std::runtime_error(ss.str());
   }

   if (matrix_in.col_dim != vector_in.Size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim << std::endl;
      ss << "The dimension of the input vector is " << vector_in.Size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }

   vector_out->value_list.resize(matrix_in.row_dim);
   using T2T3 = decltype(std::declval<T2>() * std::declval<T3>());

#ifdef _OPENMP
   const int num_threads = omp_get_max_threads();
   if (vectors_work == nullptr || static_cast<int>(vectors_work->size()) != num_threads) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "Working vector (vectors_work) must be arrays of the number of "
            "parallel threads";
      throw std::runtime_error(ss.str());
   }

#pragma omp parallel
   {
      const int thread_num = omp_get_thread_num();
#pragma omp for schedule(guided)
      for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
         const T3 temp_vec_in = vector_in.value_list[i];
         T2T3 temp_val = matrix_in.val[matrix_in.row[i + 1] - 1] * temp_vec_in;
         for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
            temp_val += matrix_in.val[j] * vector_in.value_list[matrix_in.col[j]];
            (*vectors_work)[thread_num][matrix_in.col[j]] += matrix_in.val[j] * temp_vec_in;
         }
         (*vectors_work)[thread_num][i] += temp_val;
      }
   }

#pragma omp parallel for
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      T2T3 temp_val = 0.0;
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         temp_val += (*vectors_work)[thread_num][i];
         (*vectors_work)[thread_num][i] = 0.0;
      }
      vector_out->value_list[i] = temp_val * coeff;
   }

#else
   vector_out->Fill(0.0);
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      const T3 temp_vec_in = vector_in.val[i];
      T2T3 temp_val = matrix_in.val[matrix_in.row[i + 1] - 1] * temp_vec_in;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j] * vector_in.val[matrix_in.col[j]];
         vector_out->val[matrix_in.col[j]] += matrix_in.val[j] * temp_vec_in;
      }
      vector_out->val[i] += temp_val * coeff;
   }

#endif
}

//! @brief Calculate matrix vector product.
//! \f$ \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! @tparam T1 The value type of the coefficient \f$ c\f$.
//! @tparam T2 The value type of the matirx \f$ \hat{M}\f$.
//! @tparam T3 The value type of the vector \f$ \boldsymbol{v}\f$.
//! @param coeff The coefficient \f$ c\f$.
//! @param matrix_in The matrix \f$ \hat{M} \f$.
//! @param vector_in The vector \f$ \boldsymbol{v} \f$.
//! @return The result of matrix vector product \f$ \boldsymbol{v}_{\rm out} =
//! c\hat{M}\cdot\boldsymbol{v}\f$.
template <typename T1, typename T2, typename T3>
auto CalculateMatrixVectorProduct(const T1 coeff, const CRS<T2> &matrix_in, const BraketVector<T3> &vector_in)
    -> BraketVector<decltype(std::declval<T1>() * std::declval<T2>() * std::declval<T3>())> {
   BraketVector<decltype(std::declval<T1>() * std::declval<T2>() * std::declval<T3>())> vector_out;
   CalculateMatrixVectorProduct(&vector_out, coeff, matrix_in, vector_in);
   return vector_out;
}

//! @brief Calculate matrix vector product.
//! \f$ \boldsymbol{v}_{\rm out} = c\hat{M}\cdot\boldsymbol{v}\f$.
//! Note that the matrix \f$ \hat{M}\f$ must be symmetric and their elements are
//! stored only in lower triangle elements. In addition, the diagonal elements
//! must also stored even if they are zero.
//! @tparam T1 The value type of the coefficient \f$ c\f$.
//! @tparam T2 The value type of the matirx \f$ \hat{M}\f$.
//! @tparam T3 The value type of the vector \f$ \boldsymbol{v}\f$.
//! @param coeff The coefficient \f$ c\f$.
//! @param matrix_in The matrix \f$ \hat{M} \f$.
//! @param vector_in The vector \f$ \boldsymbol{v} \f$.
//! @param vectors_work The pointer of a temporary two-dimensional array used
//! for calculations. Their elements must be zero. Note that this working array
//! is used only when openmp is active.
//! @return The matrix vector product \f$ \boldsymbol{v}_{\rm out} =
//! c\hat{M}\cdot\boldsymbol{v}\f$.
template <typename T1, typename T2, typename T3>
auto CalculateSymmetricMatrixVectorProduct(
    const T1 coeff, const CRS<T2> &matrix_in, const BraketVector<T3> &vector_in,
    std::vector<std::vector<decltype(std::declval<T2>() * std::declval<T3>())>> *vectors_work = nullptr)
    -> BraketVector<decltype(std::declval<T1>() * std::declval<T2>() * std::declval<T3>())> {
   BraketVector<decltype(std::declval<T1>() * std::declval<T2>() * std::declval<T3>())> vector_out;
   CalculateSymmetricMatrixVectorProduct(&vector_out, coeff, matrix_in, vector_in, vectors_work);
   return vector_out;
}

}  // namespace blas
}  // namespace compnal

#endif /* COMPNAL_BLAS_MATRIX_VECTOR_OPERATION_HPP_ */
