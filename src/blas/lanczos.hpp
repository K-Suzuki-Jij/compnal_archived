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
//  Created by Kohei Suzuki on 2022/03/05.
//

#ifndef COMPNAL_BLAS_LANCZOS_HPP_
#define COMPNAL_BLAS_LANCZOS_HPP_

#include "parameters.hpp"
#include "lapack.hpp"

namespace compnal {
namespace blas {

template <typename RealType>
void EigendecompositionLanczos(RealType                *target_value_out,
                               BraketVector<RealType>  *target_vector_out,
                               const CRS<RealType>     &matrix_in,
                               const std::vector<BraketVector<RealType>> &subspace_vectors = {},
                               const ParametersLanczos<RealType> &params = ParametersLanczos<RealType>()
                               ) {
   
   const auto start = std::chrono::system_clock::now();
   
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in.row_dim << ", col=" << matrix_in.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (matrix_in.row_dim == 0) {
      *target_value_out = 0.0;
      target_vector_out->Free();
      if (params.flag_display_info) {
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of Eigendecomposition by Lanczos:" << time_sec << "[sec]" << std::endl;
      }
      return;
   }
   
   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(target_value_out, &temp_gs_vector_out, matrix_in, static_cast<int>(subspace_vectors.size()));
      target_vector_out->Assign(temp_gs_vector_out);
      if (params.flag_display_info) {
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of Eigendecomposition by Lanczos:" << time_sec << "[sec]" << std::endl;
      }
      return;
   }
   
   int converge_step_number = 0;
   const std::int64_t dim = matrix_in.row_dim;
   RealType residual_error_final = 0.0;
   BraketVector<RealType> vector_0(dim);
   BraketVector<RealType> vector_1(dim);
   BraketVector<RealType> vector_2(dim);
   
   std::vector<std::vector<RealType>> rits_vector;
   std::vector<RealType> krylov_eigen_vector;
   std::vector<RealType> krylov_eigen_value(params.max_step + 1);
   std::vector<RealType> diagonal_value;
   std::vector<RealType> off_diagonal_value;
   std::vector<std::vector<RealType>> vectors_work;
   
   if (params.flag_symmetric_crs) {
#ifdef _OPENMP
      vectors_work = std::vector<std::vector<RealType>>(omp_get_max_threads(), std::vector<RealType>(dim));
#endif
   }
   
   std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
   const unsigned int seed = std::random_device()();
   std::mt19937 random_number_engine;
   
   if (params.flag_use_initial_vec) {
      vector_0.Assign(*target_vector_out);
   }
   else {
      random_number_engine.seed(seed);
      for (std::int64_t i = 0; i < dim; ++i) {
         vector_0.value_list[i] = uniform_rand(random_number_engine);
      }
   }
   
   Orthonormalize(&vector_0, subspace_vectors);
   
   if (params.flag_store_vec) {
      rits_vector.push_back(vector_0.value_list);
   }
   
   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
   }
   else {
      CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
   }
   
   diagonal_value.push_back(vector_0*vector_1);
   krylov_eigen_value[0] = diagonal_value[0];
   CalculateVectorVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
   
   for (int step = 0; step < params.max_step; ++step) {
      vector_2.Assign(vector_1);
      off_diagonal_value.push_back(vector_2.CalculateL2Norm());
      Orthonormalize(&vector_2, subspace_vectors);
      
      if (params.flag_store_vec) {
         rits_vector.push_back(vector_2.value_list);
      }
      
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
      }
      
      diagonal_value.push_back(vector_1*vector_2);
      
      if (step >= params.min_step) {
         LapackDstev(&krylov_eigen_value[step + 1], &krylov_eigen_vector, diagonal_value, off_diagonal_value);
         const RealType residual_error = std::abs(krylov_eigen_value[step + 1] - krylov_eigen_value[step]);
         
         if (params.flag_display_info) {
            std::cout << "\rLanczos_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
            std::cout << residual_error << std::flush;
         }
         if (residual_error < params.acc) {
            *target_value_out    = krylov_eigen_value[step + 1];
            converge_step_number = step + 1;
            residual_error_final = residual_error;
            break;
         }
      }
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         vector_1.value_list[i] -= diagonal_value[step + 1]*vector_2.value_list[i] + off_diagonal_value[step]*vector_0.value_list[i];
         vector_0.value_list[i]  = vector_2.value_list[i];
      }
   }
   
   if (converge_step_number <= 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Does not converge" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   target_vector_out->value_list.resize(dim);
   
   if (params.flag_store_vec) {
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         RealType temp_val = 0.0;
         for (int j = 0; j <= converge_step_number; ++j) {
            temp_val += krylov_eigen_vector[j]*rits_vector[j][i];
         }
         target_vector_out->value_list[i] = temp_val;
      }
      target_vector_out->Normalize();
   }
   else {
      if (params.flag_use_initial_vec) {
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.value_list[i] = target_vector_out->value_list[i];
            target_vector_out->value_list[i] = 0.0;
         }
      }
      else {
         random_number_engine.seed(seed);
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.value_list[i] = uniform_rand(random_number_engine);
         }
      }
      
      Orthonormalize(&vector_0, subspace_vectors);
      
      CalculateVectorVectorSum(target_vector_out, 1.0, *target_vector_out, krylov_eigen_vector[0], vector_0);
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
      }
      
      CalculateVectorVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
      
      for (int step = 1; step <= converge_step_number; ++step) {
         vector_2.Assign(vector_1);
         Orthonormalize(&vector_2, subspace_vectors);
         CalculateVectorVectorSum(target_vector_out, 1.0, *target_vector_out, krylov_eigen_vector[step], vector_2);
         
         if (params.flag_symmetric_crs) {
            CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
         }
         else {
            CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
         }
         
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_1.value_list[i] -= diagonal_value[step]*vector_2.value_list[i] + off_diagonal_value[step - 1]*vector_0.value_list[i];
            vector_0.value_list[i]  = vector_2.value_list[i];
         }
         if (params.flag_display_info) {
            std::cout << "\rLanczos_Vec_Step:" << step << "/" << converge_step_number << std::string(5, ' ') << std::flush;
         }
      }
      target_vector_out->Normalize();
   }
   if (params.flag_display_info) {
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of diagonalization:" << time_sec << "[sec]" << std::flush;
      std::cout << " (" << residual_error_final << ")" << std::flush;
      std::cout << std::endl;
   }
   
}


} // namespace blas
} // namespace compnel


#endif /* COMPNAL_BLAS_LANCZOS_HPP_ */
