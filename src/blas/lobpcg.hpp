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
//  Created by Kohei Suzuki on 2022/03/08.
//

#ifndef COMPNAL_BLAS_LOBPCG_HPP_
#define COMPNAL_BLAS_LOBPCG_HPP_

#include <chrono>

#include "braket_vector.hpp"
#include "compressed_row_storage.hpp"
#include "lapack.hpp"
#include "orthonormalize.hpp"
#include "parameters.hpp"

namespace compnal {
namespace blas {

template <typename RealType>
void EigendecompositionLOBPCG(RealType *gs_value_out, BraketVector<RealType> *gs_vector_out,
                              const CRS<RealType> &matrix_in, const bool flag_display_info = true,
                              const ParametersLanczos<RealType> &params = ParametersLanczos<RealType>()) {
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in.row_dim << ", col=" << matrix_in.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }

   const auto start = std::chrono::system_clock::now();

   if (matrix_in.row_dim == 0) {
      *gs_value_out = 0.0;
      gs_vector_out->Free();
      if (flag_display_info) {
         const auto time_count =
             std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec = static_cast<double>(time_count) / TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of Eigendecomposition by LOBPCG:" << time_sec << "[sec]" << std::endl;
      }
      return;
   }

   if (matrix_in.row_dim == 1) {
      *gs_value_out = matrix_in.val[0];
      gs_vector_out->value_list.resize(1);
      gs_vector_out->value_list[0] = 1.0;
      if (flag_display_info) {
         const auto time_count =
             std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec = static_cast<double>(time_count) / TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of Eigendecomposition by LOBPCG:" << time_sec << "[sec]" << std::endl;
      }
      return;
   }

   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(gs_value_out, &temp_gs_vector_out, matrix_in);
      gs_vector_out->Assign(temp_gs_vector_out);
      if (flag_display_info) {
         const auto time_count =
             std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec = static_cast<double>(time_count) / TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of Eigendecomposition by LOBPCG:" << time_sec << "[sec]" << std::endl;
      }
      return;
   }

   const std::int64_t dim = matrix_in.row_dim;
   BraketVector<RealType> v0(dim);
   BraketVector<RealType> v1(dim);
   BraketVector<RealType> v2(dim);
   BraketVector<RealType> w0(dim);
   BraketVector<RealType> w1(dim);
   BraketVector<RealType> w2(dim);
   BraketVector<RealType> r(dim);

   std::vector<std::vector<RealType>> vectors_work;

   if (params.flag_symmetric_crs) {
#ifdef _OPENMP
      vectors_work = std::vector<std::vector<RealType>>(omp_get_max_threads(), std::vector<RealType>(dim));
#endif
   }

   if (params.flag_use_initial_vec) {
      v0.Assign(*gs_vector_out);
   } else {
      std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
      std::mt19937 random_number_engine;
      random_number_engine.seed(std::random_device()());
      for (std::int64_t i = 0; i < dim; ++i) {
         v0.value_list[i] = uniform_rand(random_number_engine);
      }
   }

   v0.Normalize();

   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&w0, 1.0, matrix_in, v0, &vectors_work);
   } else {
      CalculateMatrixVectorProduct(&w0, 1.0, matrix_in, v0);
   }

   for (int step = 0; step < params.max_step; ++step) {
      RealType a0 = 0.0;
      RealType a1 = 0.0;
      RealType a2 = 0.0;
      RealType a3 = 0.0;
      RealType a4 = 0.0;
      RealType a5 = 0.0;

      RealType b0 = 1.0;
      RealType b1 = 0.0;
      RealType b2 = 1.0;
      RealType b3 = 0.0;
      RealType b4 = 0.0;
      RealType b5 = 1.0;

#pragma omp parallel for reduction(+ : a0, a1, a2, a3, a4, a5, b1, b3, b4)
      for (std::int64_t i = 0; i < dim; ++i) {
         a0 += v0.value_list[i] * w0.value_list[i];
         a1 += v0.value_list[i] * w1.value_list[i];
         a2 += v1.value_list[i] * w1.value_list[i];
         a3 += v0.value_list[i] * w2.value_list[i];
         a4 += v1.value_list[i] * w2.value_list[i];
         a5 += v2.value_list[i] * w2.value_list[i];

         b1 += v0.value_list[i] * v1.value_list[i];
         b3 += v0.value_list[i] * v2.value_list[i];
         b4 += v1.value_list[i] * v2.value_list[i];
      }

      std::vector<RealType> k_e_val(3);
      std::vector<RealType> k_e_vec(9);

      int krylov_dim = 3;
      std::vector<RealType> mat_a = {a0, a1, a2, a3, a4, a5};
      std::vector<RealType> mat_b = {b0, b1, b2, b3, b4, b5};

      if (step <= 1) {
         krylov_dim = step + 1;
      }

      LapackDspgv(&k_e_val, &k_e_vec, 1, krylov_dim, mat_a, mat_b);

      RealType residual_error = 0.0;

#pragma omp parallel for reduction(+ : residual_error)
      for (std::int64_t i = 0; i < dim; ++i) {
         r.value_list[i] =
             k_e_vec[0] * w0.value_list[i] + k_e_vec[1] * w1.value_list[i] + k_e_vec[2] * w2.value_list[i] -
             k_e_val[0] *
                 (k_e_vec[0] * v0.value_list[i] + k_e_vec[1] * v1.value_list[i] + k_e_vec[2] * v2.value_list[i]);
         residual_error += r.value_list[i] * r.value_list[i];
      }

      if (flag_display_info) {
         std::cout << "\rLOBPCG_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
         std::cout << residual_error << std::flush;
      }

      if (residual_error < params.acc) {
         *gs_value_out = k_e_val[0];
         gs_vector_out->value_list.resize(dim);

#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            gs_vector_out->value_list[i] =
                k_e_vec[0] * w0.value_list[i] + k_e_vec[1] * w1.value_list[i] + k_e_vec[2] * w2.value_list[i];
         }
         if (flag_display_info) {
            gs_vector_out->Normalize();
            if (flag_display_info) {
               const auto time_count =
                   std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start)
                       .count();
               const double time_sec = static_cast<double>(time_count) / TIME_UNIT_CONSTANT;
               std::cout << std::defaultfloat << std::setprecision(8)
                         << "\rElapsed time of diagonalization:" << time_sec << "[sec]" << std::flush;
               std::cout << " (" << residual_error << ")" << std::flush;
               std::cout << std::endl;
            }
         }
         return;
      }

#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         v0.value_list[i] =
             k_e_vec[0] * v0.value_list[i] + k_e_vec[1] * v1.value_list[i] + k_e_vec[2] * v2.value_list[i];
         v2.value_list[i] = k_e_vec[1] * v1.value_list[i] + k_e_vec[2] * v2.value_list[i];
         v1.value_list[i] = r.value_list[i];
      }

      RealType alpha = 1.0 / (v0 * v0);
      RealType beta = v2 * v2;

      v0.Normalize();
      v1.Normalize();

      if (beta > 0) {
         beta = 1.0 / std::sqrt(beta);
         v2.Normalize();
      } else {
         beta = 0.0;
      }

#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         w0.value_list[i] =
             alpha * (k_e_vec[0] * w0.value_list[i] + k_e_vec[1] * w1.value_list[i] + k_e_vec[2] * w2.value_list[i]);
         w2.value_list[i] = beta * (k_e_vec[1] * w1.value_list[i] + k_e_vec[2] * w2.value_list[i]);
      }

      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&w1, 1.0, matrix_in, v1, &vectors_work);
      } else {
         CalculateMatrixVectorProduct(&w1, 1.0, matrix_in, v1);
      }
   }

   std::stringstream ss;
   ss << "Error in " << __func__ << std::endl;
   ss << "Does not converge" << std::endl;
   throw std::runtime_error(ss.str());
}

}  // namespace blas
}  // namespace compnal

#endif /* COMPNAL_BLAS_LOBPCG_HPP_ */
