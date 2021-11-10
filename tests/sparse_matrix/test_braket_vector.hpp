//
//  test_sparse_matrix.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef TEST_SPARSE_MATRIX_BRAKET_VECTOR_HPP_
#define TEST_SPARSE_MATRIX_BRAKET_VECTOR_HPP_

#include "../../src/sparse_matrix/all.hpp"
#include <gtest/gtest.h>

TEST(MatrixVectorOperatopn, MatrixVectorProduct) {
   compnal::sparse_matrix::CRS<double> matrix(std::vector<std::vector<double>>{
      {1.0, 2.0, 0.0},
      {2.0, 3.0, 4.0},
      {0.0, 4.0, 5.0}
   });
   compnal::sparse_matrix::BraketVector<double> vec(std::vector<double>{1.0, 1.0, 1.0});
   const auto vec_out = compnal::sparse_matrix::CalculateMatrixVectorProduct(1.0, matrix, vec);
   vec_out.Print();
   
}

TEST(MatrixVectorOperatopn, SyymetricMatrixVectorProduct) {
   compnal::sparse_matrix::CRS<double> matrix(std::vector<std::vector<double>>{
      {1.0, 0.0, 0.0},
      {2.0, 3.0, 0.0},
      {0.0, 4.0, 5.0}
   });
   matrix.PrintInfo();
   compnal::sparse_matrix::BraketVector<double> vec(std::vector<double>{1.0, 1.0, 1.0});
   const auto vec_out = compnal::sparse_matrix::CalculateSymmetricMatrixVectorProduct(1.0, matrix, vec);
   vec_out.Print();
   
}


#endif /* TEST_SPARSE_MATRIX_BRAKET_VECTOR_HPP_ */
