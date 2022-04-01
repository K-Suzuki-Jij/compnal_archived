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
//  Created by Kohei Suzuki on 2021/11/12.
//

#ifndef COMPNAL_BLAS_CONJUGATE_GRADIENT_HPP_
#define COMPNAL_BLAS_CONJUGATE_GRADIENT_HPP_

#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"

#include <vector>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace blas {

template<typename RealType>
void ConjugateGradient(BraketVector<RealType> *vec_out,
                       const CRS<RealType>    &matrix_in,
                       const BraketVector<RealType> &vec_in,
                       const std::vector<BraketVector<RealType>> &subspace_vectors = {},
                       const bool flag_display_info = true,
                       const ParametersCG<RealType> &params = ParametersCG<RealType>()
                       ) {
   
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in.row_dim << ", col=" << matrix_in.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (static_cast<std::int64_t>(vec_in.value_list.size()) != matrix_in.row_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Matrix vector product (Ax=b) cannot be defined." << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   const auto start = std::chrono::system_clock::now();
   const std::int64_t dim = matrix_in.row_dim;
   BraketVector<RealType> rrr(dim);
   BraketVector<RealType> ppp(dim);
   BraketVector<RealType> yyy(dim);
   std::vector<std::vector<RealType>> vectors_work;
   
   if (params.flag_symmetric_crs) {
#ifdef _OPENMP
      vectors_work = std::vector<std::vector<RealType>>(omp_get_max_threads(), std::vector<RealType>(dim));
#endif
   }
   
   if (params.flag_use_initial_vec) {
      if (static_cast<std::int64_t>(vec_out->value_list.size()) != dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The dimension of the initial vector is not equal to that of the input matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   else {
      std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
      std::mt19937 random_number_engine;
      random_number_engine.seed(std::random_device()());
      vec_out->value_list.resize(dim);
      for (std::int64_t i = 0; i < dim; ++i) {
         vec_out->value_list[i] = uniform_rand(random_number_engine);
      }
   }
   Orthonormalize(vec_out, subspace_vectors);
   
   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&rrr, 1.0, matrix_in, *vec_out, &vectors_work);
   }
   else {
      CalculateMatrixVectorProduct(&rrr, 1.0, matrix_in, *vec_out);
   }
   
#pragma omp parallel for
   for (std::int64_t i = 0; i < dim; ++i) {
      rrr.value_list[i] = vec_in.value_list[i] - rrr.value_list[i];
      ppp.value_list[i] = rrr.value_list[i];
   }
   
   Orthonormalize(&rrr, subspace_vectors, false);
   Orthonormalize(&ppp, subspace_vectors, false);
   
   for (int step = 0; step < params.max_step; ++step) {
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&yyy, 1.0, matrix_in, ppp, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&yyy, 1.0, matrix_in, ppp);
      }
      
      Orthonormalize(&yyy, subspace_vectors, false);
      
      const RealType inner_prod = CalculateVectorVectorProduct(rrr, rrr);
      const RealType alpha      = inner_prod/CalculateVectorVectorProduct(ppp, yyy);
      
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         vec_out->value_list[i] += alpha*ppp.value_list[i];
         rrr.value_list[i]      -= alpha*yyy.value_list[i];
      }
      
      Orthonormalize(vec_out, subspace_vectors, false);
      Orthonormalize(&rrr   , subspace_vectors, false);
      
      const RealType residual_error = CalculateVectorVectorProduct(rrr, rrr);
      
      if (flag_display_info) {
         std::cout << "\rCG_Step[" << step << "]=" << std::scientific << std::setprecision(1);
         std::cout << residual_error << std::string(5, ' ') << std::flush;
      }
      
      if (residual_error < params.acc) {
         if (flag_display_info) {
            const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
            const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
            std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of conjugate_gradient:" << time_sec << "[sec]";
            std::cout << " (" << residual_error << ")" << std::flush;
            std::cout << std::endl;
         }
         return;
      }
      
      const RealType beta = residual_error/inner_prod;
      
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         ppp.value_list[i] = rrr.value_list[i] + beta*ppp.value_list[i];
      }
      Orthonormalize(&ppp, subspace_vectors, false);
   }
   
   std::stringstream ss;
   ss << "Error in " << __func__ << std::endl;
   ss << "Does not converge" << std::endl;
   throw std::runtime_error(ss.str());
   
}

} // namespace blas
} // namespace compnel

#endif /* COMPNAL_BLAS_CONJUGATE_GRADIENT_HPP_ */
