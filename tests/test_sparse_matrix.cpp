//
//  test_sparse_matrix.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/22.
//

#include "sparse_matrix.hpp"
#include <gtest/gtest.h>

TEST(CRS, ConstructorNoArguments) {
   
   cmpnal::sparse_matrix::CRS<double> matrix;
   
   EXPECT_EQ(matrix.GetRowDim(), 0);
   EXPECT_EQ(matrix.GetColDim(), 0);
   
   EXPECT_EQ(matrix.GetSizeCol(), 0);
   EXPECT_EQ(matrix.GetSizeVal(), 0);
   
   EXPECT_EQ(matrix.GetSizeRow(), 1);
   EXPECT_EQ(matrix.Row(0), 0);
   
}

TEST(CRS, ConstructorArguments) {
   int row_dim = 3;
   int col_dim = 5;
   cmpnal::sparse_matrix::CRS<double> matrix(row_dim, col_dim);
   
   EXPECT_EQ(matrix.GetRowDim(), row_dim);
   EXPECT_EQ(matrix.GetColDim(), col_dim);
   
   EXPECT_EQ(matrix.GetSizeCol(), 0);
   EXPECT_EQ(matrix.GetSizeVal(), 0);
   
   EXPECT_EQ(matrix.GetSizeRow(), row_dim + 1);
   
   for (int64_t i = 0; i < row_dim + 1; ++i) {
      EXPECT_EQ(matrix.Row(i), 0);
   }
   
}

TEST(CRS, ConstructorVecVec) {
   std::vector<std::vector<double>> vecvec = {
      {+0.0, +3.0, +1.0},
      {-2.0            },
      {+0.0, -2.5, +0.0}
   };
   
   cmpnal::sparse_matrix::CRS<double> matrix(vecvec);
   
   EXPECT_EQ(matrix.GetRowDim(), 3);
   EXPECT_EQ(matrix.GetColDim(), 3);
   
   EXPECT_EQ(matrix.GetSizeCol(), 4);
   EXPECT_EQ(matrix.GetSizeVal(), 4);
   
   EXPECT_EQ(matrix.GetSizeRow(), 4);
   
   EXPECT_EQ(matrix.Col(0), 1);
   EXPECT_EQ(matrix.Col(1), 2);
   EXPECT_EQ(matrix.Col(2), 0);
   EXPECT_EQ(matrix.Col(3), 1);
   
   EXPECT_DOUBLE_EQ(matrix.Val(0), +3.0);
   EXPECT_DOUBLE_EQ(matrix.Val(1), +1.0);
   EXPECT_DOUBLE_EQ(matrix.Val(2), -2.0);
   EXPECT_DOUBLE_EQ(matrix.Val(3), -2.5);
   
   EXPECT_EQ(matrix.Row(0), 0);
   EXPECT_EQ(matrix.Row(1), 2);
   EXPECT_EQ(matrix.Row(2), 3);
   EXPECT_EQ(matrix.Row(3), 4);

}


