//
//  test_eigendecomposition.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/02.
//

#ifndef TEST_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_
#define TEST_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_

#include "sparse_matrix.hpp"
#include <gtest/gtest.h>

template<typename RealType>
compnal::sparse_matrix::CRS<RealType> GenerateTightBindingHamiltonian(const int dim, const bool flag_symmetric = false) {
   if (dim <= 2) {
      throw std::runtime_error("Too small system size");
   }
   compnal::sparse_matrix::CRS<double> matrix;
   if (flag_symmetric) {
      matrix.val.push_back(0.0);
      matrix.col.push_back(0);
      matrix.val.push_back(-1.0);
      matrix.col.push_back(1);
      matrix.val.push_back(-1.0);
      matrix.col.push_back(dim - 1);
      matrix.UpdateRow();
      for (int i = 1; i < dim; ++i) {
         matrix.val.push_back(0.0);
         matrix.col.push_back(i);
         if (i + 1 < dim) {
            matrix.val.push_back(-1.0);
            matrix.col.push_back(i + 1);
         }
         matrix.UpdateRow();
      }
   }
   else {
      matrix.val.push_back(-1.0);
      matrix.col.push_back(1);
      matrix.val.push_back(-1.0);
      matrix.col.push_back(dim - 1);
      matrix.UpdateRow();
      for (int i = 0; i < dim - 1; ++i) {
         if (i == dim - 2) {
            matrix.val.push_back(-1.0);
            matrix.col.push_back(0);
         }
         matrix.val.push_back(-1.0);
         matrix.col.push_back(i);
         if (i + 2 < dim) {
            matrix.val.push_back(-1.0);
            matrix.col.push_back(i + 2);
         }
         matrix.UpdateRow();
      }
   }
   matrix.row_dim = dim;
   matrix.col_dim = dim;
   return matrix;
}

TEST(Lanczos, matrix1) {
   
   compnal::sparse_matrix::CRS<double> matrix;
   double gs_value = 0.0;
   compnal::sparse_matrix::BraketVector<double> gs_vector;
   compnal::sparse_matrix::ParametersLanczos params;
   compnal::sparse_matrix::EigenvalueDecompositionLanczos(&gs_value, &gs_vector, matrix, params);
   
}

TEST(Lanczos, TB1) {
   const auto matrix = GenerateTightBindingHamiltonian<double>(4);
   matrix.PrintMatrix();
   double gs_value = 0.0;
   compnal::sparse_matrix::BraketVector<double> gs_vector;
   compnal::sparse_matrix::ParametersLanczos params;
   compnal::sparse_matrix::EigenvalueDecompositionLanczos(&gs_value, &gs_vector, matrix, params);
   printf("%.15lf\n", gs_value);
}

TEST(Lanczos, TB2) {
   const auto matrix = GenerateTightBindingHamiltonian<double>(4, true);
   matrix.PrintMatrix();
   double gs_value = 0.0;
   compnal::sparse_matrix::BraketVector<double> gs_vector;
   compnal::sparse_matrix::ParametersLanczos params;
   params.flag_symmetric_crs = true;
   compnal::sparse_matrix::EigenvalueDecompositionLanczos(&gs_value, &gs_vector, matrix, params);
   printf("%.15lf\n", gs_value);
}


#endif /* TEST_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_ */
