//
//  eigendecomposition.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#ifndef COMPNAL_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_
#define COMPNAL_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_

#include <vector>
#include <random>

#include "lapack.hpp"
#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"
#include "orthonormalize.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


namespace compnal {
namespace sparse_matrix {

const int TIME_UNIT_CONSTANT = 1000*1000;

struct ParametersLanczos {
   int    min_step = 0;
   int    max_step = 1000;
   double acc      = std::pow(10, -14);
   bool   flag_use_initial_vec = false;
   bool   flag_store_vec       = false;
   bool   flag_output_info     = true;
   bool   flag_symmetric_crs   = false;
};

struct ParametersCG {
   int    max_step = 1000;
   double acc      = std::pow(10, -7);
   bool   flag_use_initial_vec = false;
   bool   flag_output_info     = true;
   bool   flag_symmetric_crs   = false;
};

struct ParametersII {
   ParametersII() {
      cg.flag_use_initial_vec = true;
   }
   
   int    max_step = 3;
   double acc      = std::pow(10, -7);
   double diag_add = std::pow(10, -11);
   bool   flag_output_info = true;
   ParametersCG cg;
};

struct ParametersAll {
   ParametersLanczos lanczos;
   ParametersII ii;
};

template <typename RealType>
std::pair<int, double> EigenvalueDecompositionLanczos(RealType                *target_value_out,
                                                      BraketVector<RealType>  *target_vector_out,
                                                      const CRS<RealType>     &matrix_in,
                                                      const int               target_level,
                                                      const std::vector<BraketVector<RealType>> &subspace_vectors,
                                                      const ParametersLanczos &params = ParametersLanczos()
                                                      ) {
   
   if (target_level == 0) {
      return EigenvalueDecompositionLanczos(target_value_out, target_vector_out, matrix_in, params);
   }
   
   if (target_level < 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Invalid target_level" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
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
      target_vector_out->Clear();
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(target_value_out, &temp_gs_vector_out, matrix_in, target_level);
      target_vector_out->Assign(temp_gs_vector_out);
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (target_level != static_cast<int>(subspace_vectors.size())) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Invalid size of subspace_vectors" << std::endl;
      throw std::runtime_error(ss.str());
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
         vector_0.val[i] = uniform_rand(random_number_engine);
      }
   }
   
   vector_0.Normalize();
   Orthonormalize(&vector_0, subspace_vectors);
   
   if (params.flag_store_vec) {
      rits_vector.push_back(vector_0.val);
   }
   
   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
   }
   else {
      CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
   }
   
   diagonal_value.push_back(CalculateInnerProduct(vector_0, vector_1));
   krylov_eigen_value[0] = diagonal_value[0];
   CalculateVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
   
   for (int step = 0; step < params.max_step; ++step) {
      vector_2.Assign(vector_1);
      off_diagonal_value.push_back(vector_2.L2Norm());
      vector_2.Normalize();
      Orthonormalize(&vector_2, subspace_vectors);
      
      if (params.flag_store_vec) {
         rits_vector.push_back(vector_2.val);
      }
      
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
      }
      
      diagonal_value.push_back(CalculateInnerProduct(vector_1, vector_2));
      
      if (step >= params.min_step) {
         LapackDstev(&krylov_eigen_value[step + 1], &krylov_eigen_vector, diagonal_value, off_diagonal_value);
         const RealType residual_error = std::abs(krylov_eigen_value[step + 1] - krylov_eigen_value[step]);
         
         if (params.flag_output_info) {
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
         vector_1.val[i] -= diagonal_value[step + 1]*vector_2.val[i] + off_diagonal_value[step]*vector_0.val[i];
         vector_0.val[i]  = vector_2.val[i];
      }
   }
   
   if (converge_step_number <= 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Does not converge" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   target_vector_out->val.resize(dim);
   
   if (params.flag_store_vec) {
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         RealType temp_val = 0.0;
         for (int j = 0; j <= converge_step_number; ++j) {
            temp_val += krylov_eigen_vector[j]*rits_vector[j][i];
         }
         target_vector_out->val[i] = temp_val;
      }
      target_vector_out->Normalize();
   }
   else {
      if (params.flag_use_initial_vec) {
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.val[i] = target_vector_out->val[i];
            target_vector_out->val[i] = 0.0;
         }
      }
      else {
         random_number_engine.seed(seed);
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.val[i] = uniform_rand(random_number_engine);
         }
      }
      
      vector_0.Normalize();
      Orthonormalize(&vector_0, subspace_vectors);
      
      CalculateVectorSum(target_vector_out, 1.0, *target_vector_out, krylov_eigen_vector[0], vector_0);
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
      }
      
      CalculateVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
      
      for (int step = 1; step <= converge_step_number; ++step) {
         vector_2.Assign(vector_1);
         vector_2.Normalize();
         Orthonormalize(&vector_2, subspace_vectors);
         CalculateVectorSum(target_vector_out, 1.0, *target_vector_out, krylov_eigen_vector[step], vector_2);
         
         if (params.flag_symmetric_crs) {
            CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
         }
         else {
            CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
         }
         
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_1.val[i] -= diagonal_value[step]*vector_2.val[i] + off_diagonal_value[step - 1]*vector_0.val[i];
            vector_0.val[i]  = vector_2.val[i];
         }
         if (params.flag_output_info) {
            std::cout << "\rLanczos_Vec_Step:" << step << "/" << converge_step_number << std::string(5, ' ') << std::flush;
         }
      }
      target_vector_out->Normalize();
   }
   
   const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
   const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
   std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of diagonalization:" << time_sec << "[sec]" << std::flush;
   std::cout << " (" << residual_error_final << ")" << std::flush;
   std::cout << std::endl;
   return {converge_step_number, time_sec};
   
}

template <typename RealType>
std::pair<int, double> EigenvalueDecompositionLanczos(RealType                *gs_value_out,
                                                      BraketVector<RealType>  *gs_vector_out,
                                                      const CRS<RealType>     &matrix_in,
                                                      const ParametersLanczos &params = ParametersLanczos()
                                                      ) {
   
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
      gs_vector_out->Clear();
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim == 1) {
      *gs_value_out = matrix_in.val[0];
      gs_vector_out->val.resize(1);
      gs_vector_out->val[0] = 1.0;
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(gs_value_out, &temp_gs_vector_out, matrix_in);
      gs_vector_out->Assign(temp_gs_vector_out);
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
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
      vector_0.Assign(*gs_vector_out);
   }
   else {
      random_number_engine.seed(seed);
      for (std::int64_t i = 0; i < dim; ++i) {
         vector_0.val[i] = uniform_rand(random_number_engine);
      }
   }
   
   vector_0.Normalize();
   
   if (params.flag_store_vec) {
      rits_vector.push_back(vector_0.val);
   }
   
   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
   }
   else {
      CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
   }
   
   diagonal_value.push_back(CalculateInnerProduct(vector_0, vector_1));
   krylov_eigen_value[0] = diagonal_value[0];
   CalculateVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
      
   for (int step = 0; step < params.max_step; ++step) {
      vector_2.Assign(vector_1);
      off_diagonal_value.push_back(vector_2.L2Norm());
      vector_2.Normalize();
      
      if (params.flag_store_vec) {
         rits_vector.push_back(vector_2.val);
      }
      
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
      }
      
      diagonal_value.push_back(CalculateInnerProduct(vector_1, vector_2));
      
      if (step >= params.min_step) {
         LapackDstev(&krylov_eigen_value[step + 1], &krylov_eigen_vector, diagonal_value, off_diagonal_value);
         const RealType residual_error = std::abs(krylov_eigen_value[step + 1] - krylov_eigen_value[step]);
         
         if (params.flag_output_info) {
            std::cout << "\rLanczos_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
            std::cout << residual_error << std::flush;
         }
         if (residual_error < params.acc) {
            *gs_value_out = krylov_eigen_value[step + 1];
            converge_step_number = step + 1;
            residual_error_final = residual_error;
            break;
         }
      }
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         vector_1.val[i] -= diagonal_value[step + 1]*vector_2.val[i] + off_diagonal_value[step]*vector_0.val[i];
         vector_0.val[i]  = vector_2.val[i];
      }
   }
   
   if (converge_step_number <= 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Does not converge" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   gs_vector_out->val.resize(dim);
   
   if (params.flag_store_vec) {
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         RealType temp_val = 0.0;
         for (int j = 0; j <= converge_step_number; ++j) {
            temp_val += krylov_eigen_vector[j]*rits_vector[j][i];
         }
         gs_vector_out->val[i] = temp_val;
      }
      gs_vector_out->Normalize();
   }
   else {
      if (params.flag_use_initial_vec) {
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.val[i] = gs_vector_out->val[i];
            gs_vector_out->val[i] = 0.0;
         }
      }
      else {
         random_number_engine.seed(seed);
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_0.val[i] = uniform_rand(random_number_engine);
         }
      }
      
      vector_0.Normalize();
      CalculateVectorSum(gs_vector_out, 1.0, *gs_vector_out, krylov_eigen_vector[0], vector_0);
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
      }
      
      CalculateVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
      
      for (int step = 1; step <= converge_step_number; ++step) {
         vector_2.Assign(vector_1);
         vector_2.Normalize();
         CalculateVectorSum(gs_vector_out, 1.0, *gs_vector_out, krylov_eigen_vector[step], vector_2);
         
         if (params.flag_symmetric_crs) {
            CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
         }
         else {
            CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
         }
         
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            vector_1.val[i] -= diagonal_value[step]*vector_2.val[i] + off_diagonal_value[step - 1]*vector_0.val[i];
            vector_0.val[i]  = vector_2.val[i];
         }
         if (params.flag_output_info) {
            std::cout << "\rLanczos_Vec_Step:" << step << "/" << converge_step_number << std::string(5, ' ') << std::flush;
         }
      }
      gs_vector_out->Normalize();
   }
   
   const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
   const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
   std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of diagonalization:" << time_sec << "[sec]" << std::flush;
   std::cout << " (" << residual_error_final << ")" << std::flush;
   std::cout << std::endl;
   return {converge_step_number, time_sec};
}

template <typename RealType>
std::pair<int, double> EigenvalueDecompositionLOBPCG(RealType                *gs_value_out,
                                                     BraketVector<RealType>  *gs_vector_out,
                                                     const CRS<RealType>     &matrix_in,
                                                     const ParametersLanczos &params = ParametersLanczos()) {
   
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
      gs_vector_out->Clear();
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim == 1) {
      *gs_value_out = matrix_in.val[0];
      gs_vector_out->val.resize(1);
      gs_vector_out->val[0] = 1.0;
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(gs_value_out, &temp_gs_vector_out, matrix_in);
      gs_vector_out->Assign(temp_gs_vector_out);
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   const std::int64_t dim = matrix_in.row_dim;
   BraketVector<RealType> v0(dim);
   BraketVector<RealType> v1(dim);
   BraketVector<RealType> v2(dim);
   BraketVector<RealType> w0(dim);
   BraketVector<RealType> w1(dim);
   BraketVector<RealType> w2(dim);
   BraketVector<RealType> r (dim);
   
   std::vector<std::vector<RealType>> vectors_work;

   if (params.flag_symmetric_crs) {
#ifdef _OPENMP
      vectors_work = std::vector<std::vector<RealType>>(omp_get_max_threads(), std::vector<RealType>(dim));
#endif
   }
   
   if (params.flag_use_initial_vec) {
      v0.Assign(*gs_vector_out);
   }
   else {
      std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
      std::mt19937 random_number_engine;
      random_number_engine.seed(std::random_device()());
      for (std::int64_t i = 0; i < dim; ++i) {
         v0.val[i] = uniform_rand(random_number_engine);
      }
   }
   
   v0.Normalize();
   
   if (params.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&w0, 1.0, matrix_in, v0, &vectors_work);
   }
   else {
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

#pragma omp parallel for reduction (+: a0, a1, a2, a3, a4, a5, b1, b3, b4)
      for (std::int64_t i = 0; i < dim; ++i) {
         a0 += v0.val[i]*w0.val[i];
         a1 += v0.val[i]*w1.val[i];
         a2 += v1.val[i]*w1.val[i];
         a3 += v0.val[i]*w2.val[i];
         a4 += v1.val[i]*w2.val[i];
         a5 += v2.val[i]*w2.val[i];
         
         b1 += v0.val[i]*v1.val[i];
         b3 += v0.val[i]*v2.val[i];
         b4 += v1.val[i]*v2.val[i];
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
      
#pragma omp parallel for reduction (+: residual_error)
      for (std::int64_t i = 0; i < dim; ++i) {
         r.val[i] = k_e_vec[0]*w0.val[i] + k_e_vec[1]*w1.val[i] + k_e_vec[2]*w2.val[i]
                  - k_e_val[0]*(k_e_vec[0]*v0.val[i] + k_e_vec[1]*v1.val[i] + k_e_vec[2]*v2.val[i]);
         residual_error += r.val[i]*r.val[i];
      }
      
      if (params.flag_output_info) {
         std::cout << "\rLOBPCG_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
         std::cout << residual_error << std::flush;
      }
      
      if (residual_error < params.acc) {
         *gs_value_out = k_e_val[0];
         gs_vector_out->val.resize(dim);
         
#pragma omp parallel for
         for (std::int64_t i = 0; i < dim; ++i) {
            gs_vector_out->val[i] = k_e_vec[0]*w0.val[i] + k_e_vec[1]*w1.val[i] + k_e_vec[2]*w2.val[i];
         }
         gs_vector_out->Normalize();
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
         std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of diagonalization:" << time_sec << "[sec]" << std::flush;
         std::cout << " (" << residual_error << ")" << std::flush;
         std::cout << std::endl;
         return {0, time_sec};
      }

#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         v0.val[i] = k_e_vec[0]*v0.val[i] + k_e_vec[1]*v1.val[i] + k_e_vec[2]*v2.val[i];
         v2.val[i] = k_e_vec[1]*v1.val[i] + k_e_vec[2]*v2.val[i];
         v1.val[i] = r.val[i];
      }
      
      RealType alpha = 1.0/CalculateInnerProduct(v0, v0);
      RealType beta  = CalculateInnerProduct(v2, v2);
            
      v0.Normalize();
      v1.Normalize();
      
      if (beta > 0) {
         beta = 1.0/std::sqrt(beta);
         v2.Normalize();
      }
      else {
         beta = 0.0;
      }
      
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         w0.val[i] = alpha*(k_e_vec[0]*w0.val[i] + k_e_vec[1]*w1.val[i] + k_e_vec[2]*w2.val[i]);
         w2.val[i] = beta*(k_e_vec[1]*w1.val[i] + k_e_vec[2]*w2.val[i]);
      }
      
      if (params.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&w1, 1.0, matrix_in, v1, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&w1, 1.0, matrix_in, v1);
      }
   }
   
   std::stringstream ss;
   ss << "Error in " << __func__ << std::endl;
   ss << "Does not converge" << std::endl;
   throw std::runtime_error(ss.str());
   
}

} // namespace sparse_matrix
} // namespace compnal


#endif /* eigendecomposition_hpp */
