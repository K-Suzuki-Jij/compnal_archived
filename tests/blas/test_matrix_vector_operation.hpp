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
//  Created by Kohei Suzuki on 2022/03/07.
//

#ifndef COMPNAL_TEST_MATRIX_VECTOR_OPERATION_HPP_
#define COMPNAL_TEST_MATRIX_VECTOR_OPERATION_HPP_

#include <gtest/gtest.h>

#include "../../src/blas/matrix_vector_operation.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace test {

TEST(MatrixVectorOperation, CalculateMatrixVectorProductType1Basic) {
   // float
   blas::BraketVector<float> out_f;
   blas::CRS<float> mat_f({{1.3F, 2.1F}, {-0.1F, 2.2F}});
   blas::BraketVector<float> vec_f({1.2F, 0.3F});

   blas::CalculateMatrixVectorProduct(&out_f, 2.3F, mat_f, vec_f);
   EXPECT_EQ(out_f.Size(), 2);
   EXPECT_TRUE(std::abs(out_f.value_list[0] - 5.037F) < 4 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[1] - 1.242F) < 4 * std::numeric_limits<float>::epsilon());

   // double
   blas::BraketVector<double> out_d;
   blas::CRS<double> mat_d({{1.3, 2.1}, {-0.1, 2.2}});
   blas::BraketVector<double> vec_d({1.2, 0.3});

   blas::CalculateMatrixVectorProduct(&out_d, 2.3, mat_d, vec_d);
   EXPECT_EQ(out_d.Size(), 2);
   EXPECT_TRUE(std::abs(out_d.value_list[0] - 5.037) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[1] - 1.242) < 4 * std::numeric_limits<double>::epsilon());

   // long double
   blas::BraketVector<long double> out_ld;
   blas::CRS<long double> mat_ld({{1.3L, 2.1L}, {-0.1L, 2.2L}});
   blas::BraketVector<long double> vec_ld({1.2L, 0.3L});

   blas::CalculateMatrixVectorProduct(&out_ld, 2.3L, mat_ld, vec_ld);
   EXPECT_EQ(out_ld.Size(), 2);
   EXPECT_TRUE(std::abs(out_ld.value_list[0] - 5.037L) < 4 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[1] - 1.242L) < 4 * std::numeric_limits<long double>::epsilon());
}

TEST(MatrixVectorOperation, CalculateMatrixVectorProductType1Specific) {
   blas::BraketVector<double> out_p1;
   blas::CRS<double> mat_p1({{0.0, 0.0}, {0.0, 0.0}});
   blas::BraketVector<double> vec_p1({1.2, 0.3});

   blas::CalculateMatrixVectorProduct(&out_p1, 2.3, mat_p1, vec_p1);
   EXPECT_EQ(out_p1.Size(), 2);
   EXPECT_TRUE(std::abs(out_p1.value_list[0] - 0.0) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[1] - 0.0) < 4 * std::numeric_limits<double>::epsilon());

   blas::BraketVector<double> out_p2;
   blas::CRS<double> mat_p2({{1.3, 2.1}, {-0.1, 2.2}});
   blas::BraketVector<double> vec_p2({0.0, 0.0});

   blas::CalculateMatrixVectorProduct(&out_p2, 2.3, mat_p2, vec_p2);
   EXPECT_EQ(out_p2.Size(), 2);
   EXPECT_TRUE(std::abs(out_p2.value_list[0] - 0.0) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[1] - 0.0) < 4 * std::numeric_limits<double>::epsilon());

   blas::BraketVector<double> out_p3;
   blas::CRS<double> mat_p3({{1.3, 2.1, 0.1}, {-0.1, 2.2, 0.0}});
   blas::BraketVector<double> vec_p3({1.0, 1.0});
   EXPECT_THROW(blas::CalculateMatrixVectorProduct(&out_p3, 2.3, mat_p3, vec_p3);, std::runtime_error);
}

TEST(MatrixVectorOperation, CalculateMatrixVectorProductType2Basic) {
   // float
   blas::CRS<float> mat_f({{1.3F, 2.1F}, {-0.1F, 2.2F}});
   blas::BraketVector<float> vec_f({1.2F, 0.3F});

   auto out_f = blas::CalculateMatrixVectorProduct(2.3F, mat_f, vec_f);
   EXPECT_EQ(out_f.Size(), 2);
   EXPECT_TRUE(std::abs(out_f.value_list[0] - 5.037F) < 4 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[1] - 1.242F) < 4 * std::numeric_limits<float>::epsilon());

   // double
   blas::CRS<double> mat_d({{1.3, 2.1}, {-0.1, 2.2}});
   blas::BraketVector<double> vec_d({1.2, 0.3});

   auto out_d = blas::CalculateMatrixVectorProduct(2.3, mat_d, vec_d);
   EXPECT_EQ(out_d.Size(), 2);
   EXPECT_TRUE(std::abs(out_d.value_list[0] - 5.037) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[1] - 1.242) < 4 * std::numeric_limits<double>::epsilon());

   // long double
   blas::CRS<long double> mat_ld({{1.3L, 2.1L}, {-0.1L, 2.2L}});
   blas::BraketVector<long double> vec_ld({1.2L, 0.3L});

   auto out_ld = blas::CalculateMatrixVectorProduct(2.3L, mat_ld, vec_ld);
   EXPECT_EQ(out_ld.Size(), 2);
   EXPECT_TRUE(std::abs(out_ld.value_list[0] - 5.037L) < 4 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[1] - 1.242L) < 4 * std::numeric_limits<long double>::epsilon());
}

TEST(MatrixVectorOperation, CalculateMatrixVectorProductType2Specific) {
   blas::CRS<double> mat_p1({{0.0, 0.0}, {0.0, 0.0}});
   blas::BraketVector<double> vec_p1({1.2, 0.3});

   auto out_p1 = blas::CalculateMatrixVectorProduct(2.3, mat_p1, vec_p1);
   EXPECT_EQ(out_p1.Size(), 2);
   EXPECT_TRUE(std::abs(out_p1.value_list[0] - 0.0) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[1] - 0.0) < 4 * std::numeric_limits<double>::epsilon());

   blas::CRS<double> mat_p2({{1.3, 2.1}, {-0.1, 2.2}});
   blas::BraketVector<double> vec_p2({0.0, 0.0});

   auto out_p2 = blas::CalculateMatrixVectorProduct(2.3, mat_p2, vec_p2);
   EXPECT_EQ(out_p2.Size(), 2);
   EXPECT_TRUE(std::abs(out_p2.value_list[0] - 0.0) < 4 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[1] - 0.0) < 4 * std::numeric_limits<double>::epsilon());

   blas::CRS<double> mat_p3({{1.3, 2.1, 0.1}, {-0.1, 2.2, 0.0}});
   blas::BraketVector<double> vec_p3({1.0, 1.0});
   EXPECT_THROW(blas::CalculateMatrixVectorProduct(2.3, mat_p3, vec_p3);, std::runtime_error);
}

TEST(MatrixVectorOperation, CalculateSymmetricMatrixVectorProductType1Basic) {
   // float
   blas::BraketVector<float> out_f;
   blas::CRS<float> mat_f(3, 3);
   blas::BraketVector<float> vec_f({1.2F, 0.3F, -1.0F});

   mat_f.val.push_back(+1.3F);
   mat_f.val.push_back(-0.1F);
   mat_f.val.push_back(+0.0F);
   mat_f.val.push_back(+2.1F);
   mat_f.val.push_back(+1.0F);
   mat_f.val.push_back(-0.3F);

   mat_f.col.push_back(0);
   mat_f.col.push_back(0);
   mat_f.col.push_back(1);
   mat_f.col.push_back(0);
   mat_f.col.push_back(1);
   mat_f.col.push_back(2);

   mat_f.row[0] = 0;
   mat_f.row[1] = 1;
   mat_f.row[2] = 3;
   mat_f.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<float>> work_f(omp_get_max_threads(), std::vector<float>(3));
   blas::CalculateSymmetricMatrixVectorProduct(&out_f, 2.3F, mat_f, vec_f, &work_f);
#else
   blas::CalculateSymmetricMatrixVectorProduct(&out_f, 2.3F, mat_f, vec_f);
#endif
   EXPECT_EQ(out_f.Size(), 3);
   EXPECT_TRUE(std::abs(out_f.value_list[0] + 1.311F) < 5 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[1] + 2.576F) < 5 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[2] - 7.176F) < 5 * std::numeric_limits<float>::epsilon());

   // double
   blas::BraketVector<double> out_d;
   blas::CRS<double> mat_d(3, 3);
   blas::BraketVector<double> vec_d({1.2, 0.3, -1.0});

   mat_d.val.push_back(+1.3);
   mat_d.val.push_back(-0.1);
   mat_d.val.push_back(+0.0);
   mat_d.val.push_back(+2.1);
   mat_d.val.push_back(+1.0);
   mat_d.val.push_back(-0.3);

   mat_d.col.push_back(0);
   mat_d.col.push_back(0);
   mat_d.col.push_back(1);
   mat_d.col.push_back(0);
   mat_d.col.push_back(1);
   mat_d.col.push_back(2);

   mat_d.row[0] = 0;
   mat_d.row[1] = 1;
   mat_d.row[2] = 3;
   mat_d.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_d(omp_get_max_threads(), std::vector<double>(3));
   blas::CalculateSymmetricMatrixVectorProduct(&out_d, 2.3, mat_d, vec_d, &work_d);
#else
   blas::CalculateSymmetricMatrixVectorProduct(&out_d, 2.3, mat_d, vec_d);
#endif
   EXPECT_EQ(out_d.Size(), 3);
   EXPECT_TRUE(std::abs(out_d.value_list[0] + 1.311) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[1] + 2.576) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[2] - 7.176) < 10 * std::numeric_limits<double>::epsilon());

   // long double
   blas::BraketVector<long double> out_ld;
   blas::CRS<long double> mat_ld(3, 3);
   blas::BraketVector<long double> vec_ld({1.2L, 0.3L, -1.0L});

   mat_ld.val.push_back(+1.3L);
   mat_ld.val.push_back(-0.1L);
   mat_ld.val.push_back(+0.0L);
   mat_ld.val.push_back(+2.1L);
   mat_ld.val.push_back(+1.0L);
   mat_ld.val.push_back(-0.3L);

   mat_ld.col.push_back(0);
   mat_ld.col.push_back(0);
   mat_ld.col.push_back(1);
   mat_ld.col.push_back(0);
   mat_ld.col.push_back(1);
   mat_ld.col.push_back(2);

   mat_ld.row[0] = 0;
   mat_ld.row[1] = 1;
   mat_ld.row[2] = 3;
   mat_ld.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<long double>> work_ld(omp_get_max_threads(), std::vector<long double>(3));
   blas::CalculateSymmetricMatrixVectorProduct(&out_ld, 2.3L, mat_ld, vec_ld, &work_ld);
#else
   blas::CalculateSymmetricMatrixVectorProduct(&out_ld, 2.3L, mat_ld, vec_ld);
#endif
   EXPECT_EQ(out_ld.Size(), 3);
   EXPECT_TRUE(std::abs(out_ld.value_list[0] + 1.311L) < 10 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[1] + 2.576L) < 10 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[2] - 7.176L) < 10 * std::numeric_limits<long double>::epsilon());
}

TEST(MatrixVectorOperation, CalculateSymmetricMatrixVectorProductType1Specific) {
   // M == 0.0
   blas::BraketVector<double> out_p1;
   blas::CRS<double> mat_p1(3, 3);
   blas::BraketVector<double> vec_p1({1.2, 0.3, -1.0});

   mat_p1.val.push_back(0.0);
   mat_p1.val.push_back(0.0);
   mat_p1.val.push_back(0.0);

   mat_p1.col.push_back(0);
   mat_p1.col.push_back(1);
   mat_p1.col.push_back(2);

   mat_p1.row[0] = 0;
   mat_p1.row[1] = 1;
   mat_p1.row[2] = 2;
   mat_p1.row[3] = 3;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p1(omp_get_max_threads(), std::vector<double>(3));
   blas::CalculateSymmetricMatrixVectorProduct(&out_p1, 2.3, mat_p1, vec_p1, &work_p1);
#else
   blas::CalculateSymmetricMatrixVectorProduct(&out_p1, 2.3, mat_p1, vec_p1);
#endif
   EXPECT_EQ(out_p1.Size(), 3);
   EXPECT_TRUE(std::abs(out_p1.value_list[0]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[1]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[2]) < 10 * std::numeric_limits<double>::epsilon());

   // V == 0.0
   blas::BraketVector<double> out_p2;
   blas::CRS<double> mat_p2(3, 3);
   blas::BraketVector<double> vec_p2({0.0, 0.0, 0.0});

   mat_p2.val.push_back(+1.3);
   mat_p2.val.push_back(-0.1);
   mat_p2.val.push_back(+0.0);
   mat_p2.val.push_back(+2.1);
   mat_p2.val.push_back(+1.0);
   mat_p2.val.push_back(-0.3);

   mat_p2.col.push_back(0);
   mat_p2.col.push_back(0);
   mat_p2.col.push_back(1);
   mat_p2.col.push_back(0);
   mat_p2.col.push_back(1);
   mat_p2.col.push_back(2);

   mat_p2.row[0] = 0;
   mat_p2.row[1] = 1;
   mat_p2.row[2] = 3;
   mat_p2.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p2(omp_get_max_threads(), std::vector<double>(3));
   blas::CalculateSymmetricMatrixVectorProduct(&out_p2, 2.3, mat_p2, vec_p2, &work_p2);
#else
   blas::CalculateSymmetricMatrixVectorProduct(&out_p2, 2.3, mat_p2, vec_p2);
#endif
   EXPECT_EQ(out_p2.Size(), 3);
   EXPECT_TRUE(std::abs(out_p2.value_list[0]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[1]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[2]) < 10 * std::numeric_limits<double>::epsilon());

   blas::BraketVector<double> out_p3;
   blas::CRS<double> mat_p3(3, 4);
   blas::BraketVector<double> vec_p3({1.2, 0.3, -1.0});

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p3(omp_get_max_threads(), std::vector<double>(3));
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(&out_p3, 2.3, mat_p3, vec_p3, &work_p3),
                std::runtime_error);
#else
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(&out_p3, 2.3, mat_p3, vec_p3), std::runtime_error);
#endif

   blas::BraketVector<double> out_p4;
   blas::CRS<double> mat_p4(3, 3);
   blas::BraketVector<double> vec_p4({1.2, 0.3, -1.0, 2.1});

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p4(omp_get_max_threads(), std::vector<double>(3));
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(&out_p4, 2.3, mat_p4, vec_p4, &work_p4),
                std::runtime_error);
#else
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(&out_p4, 2.3, mat_p4, vec_p4), std::runtime_error);
#endif
}

TEST(MatrixVectorOperation, CalculateSymmetricMatrixVectorProductType2Basic) {
   // float
   blas::CRS<float> mat_f(3, 3);
   blas::BraketVector<float> vec_f({1.2F, 0.3F, -1.0F});

   mat_f.val.push_back(+1.3F);
   mat_f.val.push_back(-0.1F);
   mat_f.val.push_back(+0.0F);
   mat_f.val.push_back(+2.1F);
   mat_f.val.push_back(+1.0F);
   mat_f.val.push_back(-0.3F);

   mat_f.col.push_back(0);
   mat_f.col.push_back(0);
   mat_f.col.push_back(1);
   mat_f.col.push_back(0);
   mat_f.col.push_back(1);
   mat_f.col.push_back(2);

   mat_f.row[0] = 0;
   mat_f.row[1] = 1;
   mat_f.row[2] = 3;
   mat_f.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<float>> work_f(omp_get_max_threads(), std::vector<float>(3));
   auto out_f = blas::CalculateSymmetricMatrixVectorProduct(2.3F, mat_f, vec_f, &work_f);
#else
   auto out_f = blas::CalculateSymmetricMatrixVectorProduct(2.3F, mat_f, vec_f);
#endif
   EXPECT_EQ(out_f.Size(), 3);
   EXPECT_TRUE(std::abs(out_f.value_list[0] + 1.311F) < 5 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[1] + 2.576F) < 5 * std::numeric_limits<float>::epsilon());
   EXPECT_TRUE(std::abs(out_f.value_list[2] - 7.176F) < 5 * std::numeric_limits<float>::epsilon());

   // double
   blas::CRS<double> mat_d(3, 3);
   blas::BraketVector<double> vec_d({1.2, 0.3, -1.0});

   mat_d.val.push_back(+1.3);
   mat_d.val.push_back(-0.1);
   mat_d.val.push_back(+0.0);
   mat_d.val.push_back(+2.1);
   mat_d.val.push_back(+1.0);
   mat_d.val.push_back(-0.3);

   mat_d.col.push_back(0);
   mat_d.col.push_back(0);
   mat_d.col.push_back(1);
   mat_d.col.push_back(0);
   mat_d.col.push_back(1);
   mat_d.col.push_back(2);

   mat_d.row[0] = 0;
   mat_d.row[1] = 1;
   mat_d.row[2] = 3;
   mat_d.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_d(omp_get_max_threads(), std::vector<double>(3));
   auto out_d = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_d, vec_d, &work_d);
#else
   auto out_d = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_d, vec_d);
#endif
   EXPECT_EQ(out_d.Size(), 3);
   EXPECT_TRUE(std::abs(out_d.value_list[0] + 1.311) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[1] + 2.576) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_d.value_list[2] - 7.176) < 10 * std::numeric_limits<double>::epsilon());

   // long double
   blas::CRS<long double> mat_ld(3, 3);
   blas::BraketVector<long double> vec_ld({1.2L, 0.3L, -1.0L});

   mat_ld.val.push_back(+1.3L);
   mat_ld.val.push_back(-0.1L);
   mat_ld.val.push_back(+0.0L);
   mat_ld.val.push_back(+2.1L);
   mat_ld.val.push_back(+1.0L);
   mat_ld.val.push_back(-0.3L);

   mat_ld.col.push_back(0);
   mat_ld.col.push_back(0);
   mat_ld.col.push_back(1);
   mat_ld.col.push_back(0);
   mat_ld.col.push_back(1);
   mat_ld.col.push_back(2);

   mat_ld.row[0] = 0;
   mat_ld.row[1] = 1;
   mat_ld.row[2] = 3;
   mat_ld.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<long double>> work_ld(omp_get_max_threads(), std::vector<long double>(3));
   auto out_ld = blas::CalculateSymmetricMatrixVectorProduct(2.3L, mat_ld, vec_ld, &work_ld);
#else
   auto out_ld = blas::CalculateSymmetricMatrixVectorProduct(2.3L, mat_ld, vec_ld);
#endif
   EXPECT_EQ(out_ld.Size(), 3);
   EXPECT_TRUE(std::abs(out_ld.value_list[0] + 1.311L) < 10 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[1] + 2.576L) < 10 * std::numeric_limits<long double>::epsilon());
   EXPECT_TRUE(std::abs(out_ld.value_list[2] - 7.176L) < 10 * std::numeric_limits<long double>::epsilon());
}

TEST(MatrixVectorOperation, CalculateSymmetricMatrixVectorProductType2Specific) {
   // M == 0.0
   blas::CRS<double> mat_p1(3, 3);
   blas::BraketVector<double> vec_p1({1.2, 0.3, -1.0});

   mat_p1.val.push_back(0.0);
   mat_p1.val.push_back(0.0);
   mat_p1.val.push_back(0.0);

   mat_p1.col.push_back(0);
   mat_p1.col.push_back(1);
   mat_p1.col.push_back(2);

   mat_p1.row[0] = 0;
   mat_p1.row[1] = 1;
   mat_p1.row[2] = 2;
   mat_p1.row[3] = 3;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p1(omp_get_max_threads(), std::vector<double>(3));
   auto out_p1 = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p1, vec_p1, &work_p1);
#else
   auto out_p1 = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p1, vec_p1);
#endif
   EXPECT_EQ(out_p1.Size(), 3);
   EXPECT_TRUE(std::abs(out_p1.value_list[0]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[1]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p1.value_list[2]) < 10 * std::numeric_limits<double>::epsilon());

   // V == 0.0
   blas::CRS<double> mat_p2(3, 3);
   blas::BraketVector<double> vec_p2({0.0, 0.0, 0.0});

   mat_p2.val.push_back(+1.3);
   mat_p2.val.push_back(-0.1);
   mat_p2.val.push_back(+0.0);
   mat_p2.val.push_back(+2.1);
   mat_p2.val.push_back(+1.0);
   mat_p2.val.push_back(-0.3);

   mat_p2.col.push_back(0);
   mat_p2.col.push_back(0);
   mat_p2.col.push_back(1);
   mat_p2.col.push_back(0);
   mat_p2.col.push_back(1);
   mat_p2.col.push_back(2);

   mat_p2.row[0] = 0;
   mat_p2.row[1] = 1;
   mat_p2.row[2] = 3;
   mat_p2.row[3] = 6;

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p2(omp_get_max_threads(), std::vector<double>(3));
   auto out_p2 = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p2, vec_p2, &work_p2);
#else
   auto out_p2 = blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p2, vec_p2);
#endif
   EXPECT_EQ(out_p2.Size(), 3);
   EXPECT_TRUE(std::abs(out_p2.value_list[0]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[1]) < 10 * std::numeric_limits<double>::epsilon());
   EXPECT_TRUE(std::abs(out_p2.value_list[2]) < 10 * std::numeric_limits<double>::epsilon());

   blas::CRS<double> mat_p3(3, 4);
   blas::BraketVector<double> vec_p3({1.2, 0.3, -1.0});

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p3(omp_get_max_threads(), std::vector<double>(3));
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p3, vec_p3, &work_p3), std::runtime_error);
#else
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p3, vec_p3), std::runtime_error);
#endif

   blas::CRS<double> mat_p4(3, 3);
   blas::BraketVector<double> vec_p4({1.2, 0.3, -1.0, 2.1});

#ifdef _OPENMP
   std::vector<std::vector<double>> work_p4(omp_get_max_threads(), std::vector<double>(3));
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p4, vec_p4, &work_p4), std::runtime_error);
#else
   EXPECT_THROW(blas::CalculateSymmetricMatrixVectorProduct(2.3, mat_p4, vec_p4), std::runtime_error);
#endif
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_MATRIX_VECTOR_OPERATION_HPP_ */
