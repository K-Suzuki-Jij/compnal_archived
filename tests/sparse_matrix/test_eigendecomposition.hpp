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

TEST(Lanczos, matrix1) {
   
   compnal::sparse_matrix::CRS<double> matrix;
   double gs_value = 0.0;
   compnal::sparse_matrix::BraketVector<double> gs_vector;
   compnal::sparse_matrix::ParametersLanczos params;
   compnal::sparse_matrix::EigenvalueDecompositionLanczos(&gs_value, &gs_vector, matrix, params);
   
}


#endif /* TEST_SPARSE_MATRIX_EIGENDECOMPOSITION_HPP_ */
