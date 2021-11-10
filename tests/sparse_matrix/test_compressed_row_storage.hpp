//
//  test_sparse_matrix.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef TEST_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_
#define TEST_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_

#include "../../src/sparse_matrix/all.hpp"
#include <gtest/gtest.h>

TEST(CRS, ConstructorNoArguments) {
   
   compnal::sparse_matrix::CRS<double> matrix;
   
   EXPECT_EQ(matrix.row_dim, 0);
   EXPECT_EQ(matrix.col_dim, 0);
   
   EXPECT_EQ(matrix.col.size(), 0);
   EXPECT_EQ(matrix.val.size(), 0);
   
   EXPECT_EQ(matrix.row.size(), 0);
   
}

TEST(CRS, ConstructorArguments) {
   int row_dim = 3;
   int col_dim = 5;
   compnal::sparse_matrix::CRS<double> matrix(row_dim, col_dim);
   
   EXPECT_EQ(matrix.row_dim, row_dim);
   EXPECT_EQ(matrix.col_dim, col_dim);
   
   EXPECT_EQ(matrix.col.size(), 0);
   EXPECT_EQ(matrix.val.size(), 0);
   
   EXPECT_EQ(matrix.row.size(), row_dim + 1);
   
   for (int64_t i = 0; i < row_dim + 1; ++i) {
      EXPECT_EQ(matrix.row[i], 0);
   }
   
}

TEST(CRS, ConstructorVecVec) {
   std::vector<std::vector<double>> vecvec = {
      {+0.0, +3.0, +1.0},
      {-2.0            },
      {+0.0, -2.5, +0.0}
   };
   
   compnal::sparse_matrix::CRS<double> matrix(vecvec);
   
   EXPECT_EQ(matrix.row_dim, 3);
   EXPECT_EQ(matrix.col_dim, 3);
   
   EXPECT_EQ(matrix.col.size(), 4);
   EXPECT_EQ(matrix.val.size(), 4);
   
   EXPECT_EQ(matrix.row.size(), 4);
   
   EXPECT_EQ(matrix.col[0], 1);
   EXPECT_EQ(matrix.col[1], 2);
   EXPECT_EQ(matrix.col[2], 0);
   EXPECT_EQ(matrix.col[3], 1);
   
   EXPECT_DOUBLE_EQ(matrix.val[0], +3.0);
   EXPECT_DOUBLE_EQ(matrix.val[1], +1.0);
   EXPECT_DOUBLE_EQ(matrix.val[2], -2.0);
   EXPECT_DOUBLE_EQ(matrix.val[3], -2.5);
   
   EXPECT_EQ(matrix.row[0], 0);
   EXPECT_EQ(matrix.row[1], 2);
   EXPECT_EQ(matrix.row[2], 3);
   EXPECT_EQ(matrix.row[3], 4);

}

#endif /* TEST_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_ */
