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

#include "../type/compressed_row_storage.hpp"
#include "braket_vector.hpp"

namespace compnal {
namespace blas {

template<typename RealType>
void CalculateMatrixVectorProduct(type::BraketVector<RealType> *vector_out,
                                  const RealType coeef,
                                  const type::CRS<RealType> &matrix_in,
                                  const type::BraketVector<RealType> &vector_in) {
   
   if (matrix_in.col_dim != static_cast<std::int64_t>(vector_in.val.size())) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim  << std::endl;
      ss << "The dimension of the input vector is " << vector_in.val.size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(matrix_in.row_dim);
#pragma omp parallel for
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      RealType temp = 0.0;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i+1]; ++j) {
         temp += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
      }
      vector_out->val[i] = temp*coeef;
   }
}

template<typename RealType>
void CalculateSymmetricMatrixVectorProduct(type::BraketVector<RealType> *vector_out,
                                           const RealType coeef,
                                           const type::CRS<RealType> &matrix_in,
                                           const type::BraketVector<RealType> &vector_in,
                                           std::vector<std::vector<RealType>> *vectors_work) {
   
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not symmetric" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (matrix_in.col_dim != static_cast<std::int64_t>(vector_in.val.size())) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim  << std::endl;
      ss << "The dimension of the input vector is " << vector_in.val.size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(matrix_in.row_dim);
   
#ifdef _OPENMP
   const int num_threads = omp_get_max_threads();
   if (static_cast<int>(vectors_work->size()) != num_threads) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Working vector (vectors_work) must be arrays of the number of parallel threads";
      throw std::runtime_error(ss.str());
   }
   
#pragma omp parallel for schedule (guided)
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      const int      thread_num  = omp_get_thread_num();
      const RealType temp_vec_in = vector_in.val[i];
      RealType       temp_val    = matrix_in.val[matrix_in.row[i + 1] - 1]*temp_vec_in;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
         (*vectors_work)[thread_num][matrix_in.col[j]] += matrix_in.val[j]*temp_vec_in;
      }
      (*vectors_work)[thread_num][i] += temp_val;
   }
   
#pragma omp parallel for
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      RealType temp_val = 0.0;
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         temp_val += (*vectors_work)[thread_num][i];
         (*vectors_work)[thread_num][i] = 0.0;
      }
      vector_out->val[i] = temp_val*coeef;
   }
   
#else
   vector_out->Fill(0.0);
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      const RealType temp_vec_in = vector_in.val[i];
      RealType       temp_val    = matrix_in.val[matrix_in.row[i + 1] - 1]*temp_vec_in;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
         vector_out->val[matrix_in.col[j]] += matrix_in.val[j]*temp_vec_in;
      }
      vector_out->val[i] += temp_val*coeef;
   }
   
#endif
}

template<typename RealType>
void CalculateSymmetricMatrixVectorProduct(type::BraketVector<RealType> *vector_out,
                                           const RealType coeef,
                                           const type::CRS<RealType> &matrix_in,
                                           const type::BraketVector<RealType> &vector_in) {

   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not symmetric" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (matrix_in.col_dim != static_cast<std::int64_t>(vector_in.val.size())) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim  << std::endl;
      ss << "The dimension of the input vector is " << vector_in.val.size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(matrix_in.row_dim);
   
   vector_out->Fill(0.0);
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      const RealType temp_vec_in = vector_in.val[i];
      RealType       temp_val    = matrix_in.val[matrix_in.row[i + 1] - 1]*temp_vec_in;
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
         vector_out->val[matrix_in.col[j]] += matrix_in.val[j]*temp_vec_in;
      }
      vector_out->val[i] += temp_val*coeef;
   }

}


template<typename RealType>
type::BraketVector<RealType> CalculateMatrixVectorProduct(const RealType coeef,
                                                    const type::CRS<RealType> &matrix_in,
                                                    const type::BraketVector<RealType> &vector_in) {
   type::BraketVector<RealType> vector_out;
   CalculateMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in);
   return vector_out;
}

template<typename RealType>
type::BraketVector<RealType> CalculateSymmetricMatrixVectorProduct(const RealType coeef,
                                                             const type::CRS<RealType> &matrix_in,
                                                             const type::BraketVector<RealType> &vector_in,
                                                             std::vector<std::vector<RealType>> *vectors_work) {
   type::BraketVector<RealType> vector_out;
   CalculateSymmetricMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in, vectors_work);
   return vector_out;
}

template<typename RealType>
type::BraketVector<RealType> CalculateSymmetricMatrixVectorProduct(const RealType coeef,
                                                             const type::CRS<RealType> &matrix_in,
                                                             const type::BraketVector<RealType> &vector_in) {
   type::BraketVector<RealType> vector_out;
   CalculateSymmetricMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in);
   return vector_out;
}


} // namespace blas
} // namespace compnel

#endif /* COMPNAL_BLAS_MATRIX_VECTOR_OPERATION_HPP_ */
