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

struct DiagonalizationParameters {
   //Lanczos method
   int    lanczos_min_step = 0;
   int    lanczos_max_step = 1000;
   double lanczos_acc      = std::pow(10, -14);
   
   bool   lanczos_flag_initial_vec = false;
   bool   lanczos_flag_store_vec   = false;
   
   //Conjugate gradient method
   int    cg_max_step = 1000;
   double cg_acc      = std::pow(10, -6);
   
   //Inverse iteration method
   int    ii_max_step = 2;
   double ii_acc      = std::pow(10, -9);
   double ii_diag_add = std::pow(10, -11);
   
   unsigned int seed  = 0;
   bool flag_use_seed = false;
   
   bool flag_output_converge_step = false;
   
   bool flag_symmetric_matrix = false;
};

struct DiagonalizationInfomation {
   
   void ClearLanczosRecords() {
      lanczos_actual_steps.clear();
      lanczos_actual_time.clear();
   }
   
   void ClearCGRecords() {
      cg_actual_steps.clear();
      cg_actual_time.clear();
   }
   
   void ClearIIRecords() {
      ii_actual_steps.clear();
      ii_actual_time.clear();
   }
   
   void ClearAllRecords() {
      ClearLanczosRecords();
      ClearCGRecords();
      ClearIIRecords();
   }
   
   //Lanczos method
   std::vector<int>    lanczos_actual_steps;
   std::vector<double> lanczos_actual_time;
   
   //Conjugate gradient method
   std::vector<int>    cg_actual_steps;
   std::vector<double> cg_actual_time;
   
   //Inverse iteration method
   std::vector<int>    ii_actual_steps;
   std::vector<double> ii_actual_time;
   
   const int time_unit_const = 1000*1000;
};

template <typename RealType>
void EigenvalueDecompositionLanczos(RealType               *gs_value_out,
                                    BraketVector<RealType> *gs_vector_out,
                                    const CRS<RealType>    &matrix_in,
                                    const DiagonalizationParameters &param,
                                    DiagonalizationInfomation *info
                                    ) {
   
   if (matrix_in.row_dim != matrix_in.col_dim || matrix_in.row_dim < 1 || matrix_in.col_dim < 1) {
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
      info->lanczos_actual_steps.push_back(0);
      const auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      info->lanczos_actual_time.push_back(static_cast<double>(time)/info->time_unit_const);
      return;
   }
   
   if (matrix_in.row_dim == 1) {
      *gs_value_out = matrix_in.val[0];
      gs_vector_out->Resize(1);
      gs_vector_out->Val(0) = 1;
      info->lanczos_actual_steps.push_back(0);
      const auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      info->lanczos_actual_time.push_back(static_cast<double>(time)/info->time_unit_const);
      return;
   }
   
   if (matrix_in.row_dim <= 1000) {
      info->lanczos_actual_steps.push_back(0);
      LapackDsyev(gs_value_out, gs_vector_out, matrix_in);
      const auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      info->lanczos_actual_time.push_back(static_cast<double>(time)/info->time_unit_const);
      return;
   }
   
   int converged_step_number = 0;
   const int64_t dim = matrix_in.row_dim;
   BraketVector<RealType> vector_0(dim);
   BraketVector<RealType> vector_1(dim);
   BraketVector<RealType> vector_2(dim);
   BraketVector<RealType> krylov_eigen_vector;
   
   std::vector<RealType> krylov_eigen_value(param.lanczos_max_step + 1);
   std::vector<RealType> diagonal_value    (param.lanczos_max_step + 1);
   std::vector<RealType> off_diagonal_value(param.lanczos_max_step + 1);
   
#ifdef _OPENMP
   std::vector<std::vector<RealType>> vector_work(omp_get_num_threads(), std::vector<RealType>(dim, 0.0));
#endif
   
   std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
   const unsigned int seed = param.flag_use_seed ? param.seed : std::random_device{}();
   std::mt19937 random_number_engine;
   
   if (param.lanczos_flag_initial_vec) {
      vector_0 = *gs_vector_out;
   }
   else {
      random_number_engine.seed(seed);
      for (int64_t i = 0; i < dim; ++i) {
         vector_0.Val(i) = uniform_rand(random_number_engine);
      }
   }
   
   vector_0.Normalize();
   
   if (param.flag_symmetric_matrix) {
#ifdef _OPENMP
      CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_0, &vector_work);
#else
      CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_0);
#endif
   }
   else {
      vector_1 = matrix_in*vector_0;
   }
   
   diagonal_value[0] = vector_0*vector_1;
   
   krylov_eigen_value[0] = diagonal_value[0];
   
   vector_1 -= krylov_eigen_value[0]*vector_0;
   
   for (int step = 0; step < param.lanczos_max_step; ++step) {
      vector_2 = vector_1;
      off_diagonal_value[step] = vector_2.L2Norm();
      vector_2.Normalize();
      if (param.flag_symmetric_matrix) {
#ifdef _OPENMP
         CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_2, &vector_work);
#else
         CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_2);
#endif
      }
      else {
         vector_1 = matrix_in*vector_2;
      }
      
      diagonal_value[step + 1] = vector_1*vector_2;
      if (step >= param.lanczos_min_step) {
         LapackDstev(&krylov_eigen_value[step + 1], &krylov_eigen_vector, diagonal_value, off_diagonal_value);
         const RealType residual_error = std::abs(krylov_eigen_value[step + 1] - krylov_eigen_value[step]);
         
         if (param.flag_output_converge_step) {
            std::cout << "\rLanczos_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
            std::cout << residual_error << std::flush;
         }
         
         if (residual_error < param.lanczos_acc) {
            *gs_value_out = krylov_eigen_value[step + 1];
            converged_step_number = step + 1;
            break;
         }
      }
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         vector_1.Val(i) -= diagonal_value[step + 1]*vector_2.Val(i) + off_diagonal_value[step]*vector_0.Val(i);
         vector_0.Val(i)  = vector_2.Val(i);
      }
   }
   
   if (converged_step_number <= 0) {
      std::cout << "Error in " << __func__ << std::endl;
      std::cout << "Does not converge" << std::endl;
      std::exit(0);
   }
   
   gs_vector_out->Resize(dim);
   
   if (param.lanczos_flag_store_vec) {
      //TO DO
      
   }
   else {
      if (param.lanczos_flag_initial_vec) {
#pragma omp parallel for
         for (int64_t i = 0; i < dim; ++i) {
            vector_0.Val(i) = gs_vector_out->Val(i);
            gs_vector_out->Val(i) = 0.0;
         }
      }
      else {
         random_number_engine.seed(seed);
         for (int64_t i = 0; i < dim; ++i) {
            vector_0.Val(i) = uniform_rand(random_number_engine);
         }
      }
      vector_0.Normalize();
      
      *gs_vector_out += krylov_eigen_vector.Val(0)*vector_0;
      
      if (param.flag_symmetric_matrix) {
#ifdef _OPENMP
         CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_0, &vector_work);
#else
         CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_0);
#endif
      }
      else {
         vector_1 = matrix_in*vector_0;
      }
      
      vector_1 -= krylov_eigen_value[0]*vector_0;
      
      for (int step = 1; step <= converged_step_number; ++step) {
         vector_2 = vector_1;
         vector_2.Normalize();
         *gs_vector_out += krylov_eigen_vector.Val(step)*vector_2;
         if (param.flag_symmetric_matrix) {
#ifdef _OPENMP
            CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_2, &vector_work);
#else
            CreateSymmetricMatrixVectorProduct(&vector_1, matrix_in, vector_2);
#endif
         }
         else {
            vector_1 = matrix_in*vector_2;
         }
#pragma omp parallel for
         for (int64_t i = 0; i < dim; ++i) {
            vector_1.Val(i) -= diagonal_value[step]*vector_2.Val(i) + off_diagonal_value[step - 1]*vector_0.Val(i);
            vector_0.Val(i)  = vector_2.Val(i);
         }
         if (param.flag_output_converge_step) {
            std::cout << "\rLanczos_Vec_Step:" << step << "/" << converged_step_number << std::string(5, ' ') << std::flush;
         }
      }
      gs_vector_out->Normalize();
   }
   
   info->lanczos_actual_steps.push_back(converged_step_number);
   const auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
   info->lanczos_actual_time.push_back(static_cast<double>(time)/info->time_unit_const);
   std::cout << std::endl;
}



} // namespace sparse_matrix
} // namespace compnal


#endif /* eigendecomposition_hpp */
