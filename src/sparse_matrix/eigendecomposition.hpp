//
//  eigendecomposition.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#ifndef eigendecomposition_hpp
#define eigendecomposition_hpp

#include "lapack.hpp"
#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"
#include <vector>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace compnal {
namespace sparse_matrix {

struct ParametersLanczos {
   int    min_step = 0;
   int    max_step = 1000;
   double acc      = std::pow(10, -14);
   bool   flag_use_initial_vec = false;
   bool   flag_store_vec   = false;
   bool   flag_output_info = true;
   bool   flag_symmetric_crs = false;
};

struct ParametersCG {
   int    max_step = 1000;
   double acc      = std::pow(10, -6);
   bool   flag_output_info = false;
};

struct ParametersII {
   int    max_step = 2;
   double acc      = std::pow(10, -9);
   double diag_add = std::pow(10, -11);
   bool   flag_output_info = false;
};

struct ParametersAll {
   ParametersLanczos lanczos;
   ParametersCG cg;
   ParametersII ii;
};

template <typename RealType>
std::pair<int, double> EigenvalueDecompositionLanczos(RealType                *gs_value_out,
                                                      BraketVector<RealType>  *gs_vector_out,
                                                      const CRS<RealType>     &matrix_in,
                                                      const ParametersLanczos &param = ParametersLanczos()
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
      const double time_sec   = static_cast<double>(time_count)/utility::TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim == 1) {
      *gs_value_out = matrix_in.val[0];
      gs_vector_out->val.resize(1);
      gs_vector_out->val[0] = 1.0;
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/utility::TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   if (matrix_in.row_dim <= 1000) {
      std::vector<RealType> temp_gs_vector_out;
      LapackDsyev(gs_value_out, &temp_gs_vector_out, matrix_in);
      gs_vector_out->Assign(temp_gs_vector_out);
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/utility::TIME_UNIT_CONSTANT;
      return {0, time_sec};
   }
   
   int converge_step_number = 0;
   const std::size_t dim = matrix_in.row_dim;
   BraketVector<RealType> vector_0(dim);
   BraketVector<RealType> vector_1(dim);
   BraketVector<RealType> vector_2(dim);

   std::vector<std::vector<RealType>> rits_vector;
   std::vector<RealType> krylov_eigen_vector;
   std::vector<RealType> krylov_eigen_value(param.max_step + 1);
   std::vector<RealType> diagonal_value;
   std::vector<RealType> off_diagonal_value;
   std::vector<std::vector<RealType>> vectors_work;
   
   if (param.flag_symmetric_crs) {
#ifdef _OPENMP
      vectors_work.resize(omp_get_max_threads());
#pragma omp parallel for
      for (std::size_t i = 0; i < vectors_work.size(); ++i) {
         vectors_work[i].resize(dim, 0.0);
      }
#endif
   }

      
   std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
   const unsigned int seed = std::random_device()();
   std::mt19937 random_number_engine;
   
   if (param.flag_use_initial_vec) {
      vector_0.Assign(*gs_vector_out);
   }
   else {
      random_number_engine.seed(seed);
      for (std::size_t i = 0; i < dim; ++i) {
         vector_0.val[i] = uniform_rand(random_number_engine);
      }
   }
   
   vector_0.Normalize();
   
   if (param.flag_store_vec) {
      rits_vector.push_back(vector_0.val);
   }
   
   if (param.flag_symmetric_crs) {
      CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0, &vectors_work);
   }
   else {
      CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_0);
   }
   
   diagonal_value.push_back(CalculateInnerProduct(vector_0, vector_1));
   krylov_eigen_value[0] = diagonal_value[0];
   CalculateVectorSum(&vector_1, 1.0, vector_1, -krylov_eigen_value[0], vector_0);
      
   for (int step = 0; step < param.max_step; ++step) {
      vector_2.Assign(vector_1);
      off_diagonal_value.push_back(vector_2.L2Norm());
      vector_2.Normalize();
      
      if (param.flag_store_vec) {
         rits_vector.push_back(vector_2.val);
      }
      
      if (param.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
      }
      else {
         CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
      }
      
      diagonal_value.push_back(CalculateInnerProduct(vector_1, vector_2));
      
      if (step >= param.min_step) {
         LapackDstev(&krylov_eigen_value[step + 1], &krylov_eigen_vector, diagonal_value, off_diagonal_value);
         const RealType residual_error = std::abs(krylov_eigen_value[step + 1] - krylov_eigen_value[step]);
         
         if (param.flag_output_info) {
            std::cout << "\rLanczos_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
            std::cout << residual_error << std::flush;
         }
         if (residual_error < param.acc) {
            *gs_value_out = krylov_eigen_value[step + 1];
            converge_step_number = step + 1;
            break;
         }
      }
#pragma omp parallel for
      for (std::size_t i = 0; i < dim; ++i) {
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
   
   if (param.flag_store_vec) {
#pragma omp parallel for
      for (std::size_t i = 0; i < dim; ++i) {
         RealType temp_val = 0.0;
         for (std::size_t j = 0; j <= converge_step_number; ++j) {
            temp_val += krylov_eigen_vector[j]*rits_vector[j][i];
         }
         gs_vector_out->val[i] = temp_val;
      }
      gs_vector_out->Normalize();
   }
   else {
      if (param.flag_use_initial_vec) {
#pragma omp parallel for
         for (std::size_t i = 0; i < dim; ++i) {
            vector_0.val[i] = gs_vector_out->val[i];
            gs_vector_out->val[i] = 0.0;
         }
      }
      else {
         random_number_engine.seed(seed);
         for (std::size_t i = 0; i < dim; ++i) {
            vector_0.val[i] = uniform_rand(random_number_engine);
         }
      }
      
      vector_0.Normalize();
      CalculateVectorSum(gs_vector_out, 1.0, *gs_vector_out, krylov_eigen_vector[0], vector_0);
      if (param.flag_symmetric_crs) {
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
         
         if (param.flag_symmetric_crs) {
            CalculateSymmetricMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2, &vectors_work);
         }
         else {
            CalculateMatrixVectorProduct(&vector_1, 1.0, matrix_in, vector_2);
         }
         
#pragma omp parallel for
         for (std::size_t i = 0; i < dim; ++i) {
            vector_1.val[i] -= diagonal_value[step]*vector_2.val[i] + off_diagonal_value[step - 1]*vector_0.val[i];
            vector_0.val[i]  = vector_2.val[i];
         }
         if (param.flag_output_info) {
            std::cout << "\rLanczos_Vec_Step:" << step << "/" << converge_step_number << std::string(5, ' ') << std::flush;
         }
      }
      gs_vector_out->Normalize();
   }
   
   std::cout << std::endl;
   const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
   const double time_sec   = static_cast<double>(time_count)/utility::TIME_UNIT_CONSTANT;
   return {converge_step_number, time_sec};
}



} // namespace sparse_matrix
} // namespace compnal


#endif /* eigendecomposition_hpp */
