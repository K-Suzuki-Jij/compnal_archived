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
//  Created by Kohei Suzuki on 2022/02/14.
//

#ifndef COMPNAL_TEST_COMPRESSED_ROW_STORAGE_HPP_
#define COMPNAL_TEST_COMPRESSED_ROW_STORAGE_HPP_

#include "../../src/type/compressed_row_storage.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(CRS, Addition) {
   auto expect_crs_eq = [](const auto &lhs, const auto &rhs, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(lhs)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(rhs)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(lhs, rhs);
      }
      else {
         EXPECT_EQ(lhs.row_dim, rhs.row_dim);
         EXPECT_EQ(lhs.col_dim, rhs.col_dim);
         EXPECT_EQ(lhs.row.size(), rhs.row.size());
         EXPECT_EQ(lhs.col.size(), rhs.col.size());
         EXPECT_EQ(lhs.val.size(), rhs.val.size());
         EXPECT_EQ(lhs.tag, rhs.tag);
         EXPECT_EQ(lhs.name, rhs.name);
         for (std::size_t i = 0; i < lhs.row.size(); ++i) {
            EXPECT_EQ(lhs.row.at(i), rhs.row.at(i));
         }
         for (std::size_t i = 0; i < lhs.col.size(); ++i) {
            EXPECT_EQ(lhs.col.at(i), rhs.col.at(i));
         }
         for (std::size_t i = 0; i < lhs.val.size(); ++i) {
            const auto lhs_val = std::abs(lhs.val.at(i) - rhs.val.at(i));
            const auto rhs_val = 4*prec;
            EXPECT_TRUE(lhs_val <= rhs_val) << lhs << " !<= " << rhs << " from " << line;
         }
      }
   };
   using type::CRS;
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.7L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.7 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.7F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<long long  >({{1LL }}), CRS<long double>({{2.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<long       >({{1L  }}), CRS<long double>({{2.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<int        >({{1   }}), CRS<long double>({{2.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<short      >({{1   }}), CRS<long double>({{2.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) + CRS<char       >({{1   }}), CRS<long double>({{2.3L}}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<long double>({{1.4L}}), CRS<double>({{2.7 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<double     >({{1.4 }}), CRS<double>({{2.7 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<float      >({{1.4F}}), CRS<float >({{2.7F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<long long  >({{1LL }}), CRS<double>({{2.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<long       >({{1L  }}), CRS<double>({{2.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<int        >({{1   }}), CRS<double>({{2.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<short      >({{1   }}), CRS<double>({{2.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) + CRS<char       >({{1   }}), CRS<double>({{2.3 }}), DBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<long double>({{1.4L}}), CRS<float>({{2.7F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<double     >({{1.4 }}), CRS<float>({{2.7F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<float      >({{1.4F}}), CRS<float>({{2.7F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<long long  >({{1LL }}), CRS<float>({{2.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<long       >({{1L  }}), CRS<float>({{2.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<int        >({{1   }}), CRS<float>({{2.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<short      >({{1   }}), CRS<float>({{2.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) + CRS<char       >({{1   }}), CRS<float>({{2.3F}}), FLT_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<long long  >({{1LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<long       >({{1L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<int        >({{1   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<short      >({{1   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) + CRS<char       >({{1   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);

   expect_crs_eq(CRS<long>({{1L}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<long long  >({{1LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<long       >({{1L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<int        >({{1   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<short      >({{1   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) + CRS<char       >({{1   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<int>({{1}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<long long  >({{1LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<long       >({{1L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<int        >({{1   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<short      >({{1   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) + CRS<char       >({{1   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<short>({{1}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<long long  >({{1LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<long       >({{1L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<int        >({{1   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<short      >({{1   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) + CRS<char       >({{1   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<char>({{1}}) + CRS<long double>({{1.4L}}), CRS<long double>({{2.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<double     >({{1.4 }}), CRS<double     >({{2.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<float      >({{1.4F}}), CRS<float      >({{2.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<long long  >({{1LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<long       >({{1L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<int        >({{1   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<short      >({{1   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) + CRS<char       >({{1   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);

}

TEST(CRS, Subtraction) {
   auto expect_crs_eq = [](const auto &lhs, const auto &rhs, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(lhs)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(rhs)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(lhs, rhs);
      }
      else {
         EXPECT_EQ(lhs.row_dim, rhs.row_dim);
         EXPECT_EQ(lhs.col_dim, rhs.col_dim);
         EXPECT_EQ(lhs.row.size(), rhs.row.size());
         EXPECT_EQ(lhs.col.size(), rhs.col.size());
         EXPECT_EQ(lhs.val.size(), rhs.val.size());
         EXPECT_EQ(lhs.tag, rhs.tag);
         EXPECT_EQ(lhs.name, rhs.name);
         for (std::size_t i = 0; i < lhs.row.size(); ++i) {
            EXPECT_EQ(lhs.row.at(i), rhs.row.at(i));
         }
         for (std::size_t i = 0; i < lhs.col.size(); ++i) {
            EXPECT_EQ(lhs.col.at(i), rhs.col.at(i));
         }
         for (std::size_t i = 0; i < lhs.val.size(); ++i) {
            const auto lhs_val = std::abs(lhs.val.at(i) - rhs.val.at(i));
            const auto rhs_val = 4*prec;
            EXPECT_TRUE(lhs_val <= rhs_val) << lhs << " !<= " << rhs << " from " << line;
         }
      }
   };
   using type::CRS;
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.1L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.1 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.1F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<long long  >({{1LL }}), CRS<long double>({{+0.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<long       >({{1L  }}), CRS<long double>({{+0.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<int        >({{1   }}), CRS<long double>({{+0.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<short      >({{1   }}), CRS<long double>({{+0.3L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) - CRS<char       >({{1   }}), CRS<long double>({{+0.3L}}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<long double>({{1.4L}}), CRS<double>({{-0.1 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<double     >({{1.4 }}), CRS<double>({{-0.1 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<float      >({{1.4F}}), CRS<float >({{-0.1F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<long long  >({{1LL }}), CRS<double>({{+0.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<long       >({{1L  }}), CRS<double>({{+0.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<int        >({{1   }}), CRS<double>({{+0.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<short      >({{1   }}), CRS<double>({{+0.3 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) - CRS<char       >({{1   }}), CRS<double>({{+0.3 }}), DBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<long double>({{1.4L}}), CRS<float>({{-0.1F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<double     >({{1.4 }}), CRS<float>({{-0.1F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<float      >({{1.4F}}), CRS<float>({{-0.1F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<long long  >({{1LL }}), CRS<float>({{+0.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<long       >({{1L  }}), CRS<float>({{+0.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<int        >({{1   }}), CRS<float>({{+0.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<short      >({{1   }}), CRS<float>({{+0.3F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) - CRS<char       >({{1   }}), CRS<float>({{+0.3F}}), FLT_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<long long  >({{2LL }}), CRS<long long  >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<long       >({{2L  }}), CRS<long       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<int        >({{2   }}), CRS<int        >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<short      >({{2   }}), CRS<short      >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) - CRS<char       >({{2   }}), CRS<char       >({{-1LL }}), LDBL_EPSILON, __LINE__);

   expect_crs_eq(CRS<long>({{1L}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<long long  >({{2LL }}), CRS<long long  >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<long       >({{2L  }}), CRS<long       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<int        >({{2   }}), CRS<int        >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<short      >({{2   }}), CRS<short      >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) - CRS<char       >({{2   }}), CRS<char       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<int>({{1}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<long long  >({{2LL }}), CRS<long long  >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<long       >({{2L  }}), CRS<long       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<int        >({{2   }}), CRS<int        >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<short      >({{2   }}), CRS<short      >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) - CRS<char       >({{2   }}), CRS<char       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<short>({{1}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<long long  >({{2LL }}), CRS<long long  >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<long       >({{2L  }}), CRS<long       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<int        >({{2   }}), CRS<int        >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<short      >({{2   }}), CRS<short      >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) - CRS<char       >({{2   }}), CRS<char       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<char>({{1}}) - CRS<long double>({{1.4L}}), CRS<long double>({{-0.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<double     >({{1.4 }}), CRS<double     >({{-0.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<float      >({{1.4F}}), CRS<float      >({{-0.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<long long  >({{2LL }}), CRS<long long  >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<long       >({{2L  }}), CRS<long       >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<int        >({{2   }}), CRS<int        >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<short      >({{2   }}), CRS<short      >({{-1LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) - CRS<char       >({{2   }}), CRS<char       >({{-1LL }}), LDBL_EPSILON, __LINE__);
}

TEST(CRS, MultiplicationMatrixMatrix) {
   auto expect_crs_eq = [](const auto &lhs, const auto &rhs, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(lhs)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(rhs)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(lhs, rhs);
      }
      else {
         EXPECT_EQ(lhs.row_dim, rhs.row_dim);
         EXPECT_EQ(lhs.col_dim, rhs.col_dim);
         EXPECT_EQ(lhs.row.size(), rhs.row.size());
         EXPECT_EQ(lhs.col.size(), rhs.col.size());
         EXPECT_EQ(lhs.val.size(), rhs.val.size());
         EXPECT_EQ(lhs.tag, rhs.tag);
         EXPECT_EQ(lhs.name, rhs.name);
         for (std::size_t i = 0; i < lhs.row.size(); ++i) {
            EXPECT_EQ(lhs.row.at(i), rhs.row.at(i));
         }
         for (std::size_t i = 0; i < lhs.col.size(); ++i) {
            EXPECT_EQ(lhs.col.at(i), rhs.col.at(i));
         }
         for (std::size_t i = 0; i < lhs.val.size(); ++i) {
            const auto lhs_val = std::abs(lhs.val.at(i) - rhs.val.at(i));
            const auto rhs_val = 4*prec;
            EXPECT_TRUE(lhs_val <= rhs_val) << lhs << " !<= " << rhs << " from " << line;
         }
      }
   };
   using type::CRS;
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.82L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.82 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.82F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<long long  >({{2LL }}), CRS<long double>({{2.6L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<long       >({{2L  }}), CRS<long double>({{2.6L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<int        >({{2   }}), CRS<long double>({{2.6L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<short      >({{2   }}), CRS<long double>({{2.6L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * CRS<char       >({{2   }}), CRS<long double>({{2.6L }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<long double>({{1.4L}}), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<double     >({{1.4 }}), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<float      >({{1.4F}}), CRS<float >({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<long long  >({{2LL }}), CRS<double>({{2.6  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<long       >({{2L  }}), CRS<double>({{2.6  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<int        >({{2   }}), CRS<double>({{2.6  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<short      >({{2   }}), CRS<double>({{2.6  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * CRS<char       >({{2   }}), CRS<double>({{2.6  }}), DBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<long double>({{1.4L}}), CRS<float>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<double     >({{1.4 }}), CRS<float>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<float      >({{1.4F}}), CRS<float>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<long long  >({{2LL }}), CRS<float>({{2.6F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<long       >({{2L  }}), CRS<float>({{2.6F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<int        >({{2   }}), CRS<float>({{2.6F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<short      >({{2   }}), CRS<float>({{2.6F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * CRS<char       >({{2   }}), CRS<float>({{2.6F }}), FLT_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<long long  >({{2LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<long       >({{2L  }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<int        >({{2   }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<short      >({{2   }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * CRS<char       >({{2   }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);

   expect_crs_eq(CRS<long>({{1L}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<long long  >({{2LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<long       >({{2L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<int        >({{2   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<short      >({{2   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * CRS<char       >({{2   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<int>({{1}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<long long  >({{2LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<long       >({{2L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<int        >({{2   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<short      >({{2   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * CRS<char       >({{2   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<short>({{1}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<long long  >({{2LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<long       >({{2L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<int        >({{2   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<short      >({{2   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * CRS<char       >({{2   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<char>({{1}}) * CRS<long double>({{1.4L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<double     >({{1.4 }}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<float      >({{1.4F}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<long long  >({{2LL }}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<long       >({{2L  }}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<int        >({{2   }}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<short      >({{2   }}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * CRS<char       >({{2   }}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
}

TEST(CRS, MultiplicationMatrixScalar) {
   auto expect_crs_eq = [](const auto &lhs, const auto &rhs, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(lhs)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(rhs)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(lhs, rhs);
      }
      else {
         EXPECT_EQ(lhs.row_dim, rhs.row_dim);
         EXPECT_EQ(lhs.col_dim, rhs.col_dim);
         EXPECT_EQ(lhs.row.size(), rhs.row.size());
         EXPECT_EQ(lhs.col.size(), rhs.col.size());
         EXPECT_EQ(lhs.val.size(), rhs.val.size());
         EXPECT_EQ(lhs.tag, rhs.tag);
         EXPECT_EQ(lhs.name, rhs.name);
         for (std::size_t i = 0; i < lhs.row.size(); ++i) {
            EXPECT_EQ(lhs.row.at(i), rhs.row.at(i));
         }
         for (std::size_t i = 0; i < lhs.col.size(); ++i) {
            EXPECT_EQ(lhs.col.at(i), rhs.col.at(i));
         }
         for (std::size_t i = 0; i < lhs.val.size(); ++i) {
            const auto lhs_val = std::abs(lhs.val.at(i) - rhs.val.at(i));
            const auto rhs_val = 4*prec;
            EXPECT_TRUE(lhs_val <= rhs_val) << lhs << " !<= " << rhs << " from " << line;
         }
      }
   };
   using type::CRS;
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<long double>(1.4L), CRS<long double>({{1.82L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.82 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<float      >(1.4F), CRS<float      >({{1.82F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<long long  >(1LL ), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<long       >(1L  ), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<int        >(1   ), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<short      >(1   ), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long double>({{1.3L}}) * static_cast<char       >(1   ), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<long double>(1.4L), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<double     >(1.4 ), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<float      >(1.4F), CRS<float >({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<long long  >(1LL ), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<long       >(1L  ), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<int        >(1   ), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<short      >(1   ), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<double>({{1.3}}) * static_cast<char       >(1   ), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<long double>(1.4L), CRS<double>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<double     >(1.4 ), CRS<double>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<float      >(1.4F), CRS<float >({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<long long  >(1LL ), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<long       >(1L  ), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<int        >(1   ), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<short      >(1   ), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(CRS<float>({{1.3F}}) * static_cast<char       >(1   ), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<long double>(1.4L), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<float      >(1.4F), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<long long  >(2LL ), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<long       >(2L  ), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<int        >(2   ), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<short      >(2   ), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long long>({{1LL}}) * static_cast<char       >(2   ), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);

   expect_crs_eq(CRS<long>({{1L}}) * static_cast<long double>(1.4L), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<float      >(1.4F), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<long long  >(2LL ), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<long       >(2L  ), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<int        >(2   ), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<short      >(2   ), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<long>({{1L}}) * static_cast<char       >(2   ), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<int>({{1}}) * static_cast<long double>(1.4L), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<float      >(1.4F), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<long long  >(2LL ), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<long       >(2L  ), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<int        >(2   ), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<short      >(2   ), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<int>({{1}}) * static_cast<char       >(2   ), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<short>({{1}}) * static_cast<long double>(1.4L), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<float      >(1.4F), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<long long  >(2LL ), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<long       >(2L  ), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<int        >(2   ), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<short      >(2   ), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<short>({{1}}) * static_cast<char       >(2   ), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(CRS<char>({{1}}) * static_cast<long double>(1.4L), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<double     >(1.4 ), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<float      >(1.4F), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<long long  >(2LL ), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<long       >(2L  ), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<int        >(2   ), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<short      >(2   ), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(CRS<char>({{1}}) * static_cast<char       >(2   ), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
}

TEST(CRS, MultiplicationScalarMatrix) {
   auto expect_crs_eq = [](const auto &lhs, const auto &rhs, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(lhs)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(rhs)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(lhs, rhs);
      }
      else {
         EXPECT_EQ(lhs.row_dim, rhs.row_dim);
         EXPECT_EQ(lhs.col_dim, rhs.col_dim);
         EXPECT_EQ(lhs.row.size(), rhs.row.size());
         EXPECT_EQ(lhs.col.size(), rhs.col.size());
         EXPECT_EQ(lhs.val.size(), rhs.val.size());
         EXPECT_EQ(lhs.tag, rhs.tag);
         EXPECT_EQ(lhs.name, rhs.name);
         for (std::size_t i = 0; i < lhs.row.size(); ++i) {
            EXPECT_EQ(lhs.row.at(i), rhs.row.at(i));
         }
         for (std::size_t i = 0; i < lhs.col.size(); ++i) {
            EXPECT_EQ(lhs.col.at(i), rhs.col.at(i));
         }
         for (std::size_t i = 0; i < lhs.val.size(); ++i) {
            const auto lhs_val = std::abs(lhs.val.at(i) - rhs.val.at(i));
            const auto rhs_val = 4*prec;
            EXPECT_TRUE(lhs_val <= rhs_val) << lhs << " !<= " << rhs << " from " << line;
         }
      }
   };
   using type::CRS;
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<long double>({{1.3L}}), CRS<long double>({{1.82L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<long double>({{1.3L}}), CRS<double     >({{1.82 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<long double>({{1.3L}}), CRS<float      >({{1.82F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(1LL ) * CRS<long double>({{1.3L}}), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(1L  ) * CRS<long double>({{1.3L}}), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(1   ) * CRS<long double>({{1.3L}}), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(1   ) * CRS<long double>({{1.3L}}), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(1   ) * CRS<long double>({{1.3L}}), CRS<long double>({{1.3L }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<double>({{1.3}}), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<double>({{1.3}}), CRS<double>({{1.82 }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<double>({{1.3}}), CRS<float >({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long long  >(1LL ) * CRS<double>({{1.3}}), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(1L  ) * CRS<double>({{1.3}}), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(1   ) * CRS<double>({{1.3}}), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(1   ) * CRS<double>({{1.3}}), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(1   ) * CRS<double>({{1.3}}), CRS<double>({{1.3  }}), DBL_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<float>({{1.3F}}), CRS<double>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<float>({{1.3F}}), CRS<double>({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<float>({{1.3F}}), CRS<float >({{1.82F}}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long long  >(1LL ) * CRS<float>({{1.3F}}), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(1L  ) * CRS<float>({{1.3F}}), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(1   ) * CRS<float>({{1.3F}}), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(1   ) * CRS<float>({{1.3F}}), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(1   ) * CRS<float>({{1.3F}}), CRS<double>({{1.3F }}), FLT_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<long long>({{1LL}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<long long>({{1LL}}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<long long>({{1LL}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(2LL ) * CRS<long long>({{1LL}}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(2L  ) * CRS<long long>({{1LL}}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(2   ) * CRS<long long>({{1LL}}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(2   ) * CRS<long long>({{1LL}}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(2   ) * CRS<long long>({{1LL}}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);

   expect_crs_eq(static_cast<long double>(1.4L) * CRS<long>({{1L}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<long>({{1L}}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<long>({{1L}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(2LL ) * CRS<long>({{1L}}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(2L  ) * CRS<long>({{1L}}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(2   ) * CRS<long>({{1L}}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(2   ) * CRS<long>({{1L}}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(2   ) * CRS<long>({{1L}}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<int>({{1}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<int>({{1}}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<int>({{1}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(2LL ) * CRS<int>({{1}}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(2L  ) * CRS<int>({{1}}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(2   ) * CRS<int>({{1}}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(2   ) * CRS<int>({{1}}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(2   ) * CRS<int>({{1}}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<short>({{1}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<short>({{1}}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<short>({{1}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(2LL ) * CRS<short>({{1}}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(2L  ) * CRS<short>({{1}}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(2   ) * CRS<short>({{1}}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(2   ) * CRS<short>({{1}}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(2   ) * CRS<short>({{1}}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
   
   expect_crs_eq(static_cast<long double>(1.4L) * CRS<char>({{1}}), CRS<long double>({{1.4L}}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<double     >(1.4 ) * CRS<char>({{1}}), CRS<double     >({{1.4 }}), DBL_EPSILON , __LINE__);
   expect_crs_eq(static_cast<float      >(1.4F) * CRS<char>({{1}}), CRS<float      >({{1.4F}}), FLT_EPSILON , __LINE__);
   expect_crs_eq(static_cast<long long  >(2LL ) * CRS<char>({{1}}), CRS<long long  >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<long       >(2L  ) * CRS<char>({{1}}), CRS<long       >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<int        >(2   ) * CRS<char>({{1}}), CRS<int        >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<short      >(2   ) * CRS<char>({{1}}), CRS<short      >({{2LL }}), LDBL_EPSILON, __LINE__);
   expect_crs_eq(static_cast<char       >(2   ) * CRS<char>({{1}}), CRS<char       >({{2LL }}), LDBL_EPSILON, __LINE__);
}

TEST(CRS, MultiplicationEquality) {
   using type::CRS;
   EXPECT_TRUE(CRS<long long>({{1LL}}) == CRS<long long>({{1LL}}));
   EXPECT_TRUE(CRS<long long>({{1LL}}) == CRS<long     >({{1L }}));
   EXPECT_TRUE(CRS<long long>({{1LL}}) == CRS<int      >({{1  }}));
   EXPECT_TRUE(CRS<long long>({{1LL}}) == CRS<short    >({{1  }}));
   EXPECT_TRUE(CRS<long long>({{1LL}}) == CRS<char     >({{1  }}));

   EXPECT_TRUE(CRS<long>({{1L}}) == CRS<long long>({{1LL}}));
   EXPECT_TRUE(CRS<long>({{1L}}) == CRS<long     >({{1L }}));
   EXPECT_TRUE(CRS<long>({{1L}}) == CRS<int      >({{1  }}));
   EXPECT_TRUE(CRS<long>({{1L}}) == CRS<short    >({{1  }}));
   EXPECT_TRUE(CRS<long>({{1L}}) == CRS<char     >({{1  }}));
   
   EXPECT_TRUE(CRS<int>({{1}}) == CRS<long long>({{1LL}}));
   EXPECT_TRUE(CRS<int>({{1}}) == CRS<long     >({{1L }}));
   EXPECT_TRUE(CRS<int>({{1}}) == CRS<int      >({{1  }}));
   EXPECT_TRUE(CRS<int>({{1}}) == CRS<short    >({{1  }}));
   EXPECT_TRUE(CRS<int>({{1}}) == CRS<char     >({{1  }}));
   
   EXPECT_TRUE(CRS<short>({{1}}) == CRS<long long>({{1LL}}));
   EXPECT_TRUE(CRS<short>({{1}}) == CRS<long     >({{1L }}));
   EXPECT_TRUE(CRS<short>({{1}}) == CRS<int      >({{1  }}));
   EXPECT_TRUE(CRS<short>({{1}}) == CRS<short    >({{1  }}));
   EXPECT_TRUE(CRS<short>({{1}}) == CRS<char     >({{1  }}));
   
   EXPECT_TRUE(CRS<char>({{1}}) == CRS<long long>({{1LL}}));
   EXPECT_TRUE(CRS<char>({{1}}) == CRS<long     >({{1L }}));
   EXPECT_TRUE(CRS<char>({{1}}) == CRS<int      >({{1  }}));
   EXPECT_TRUE(CRS<char>({{1}}) == CRS<short    >({{1  }}));
   EXPECT_TRUE(CRS<char>({{1}}) == CRS<char     >({{1  }}));
   
   
   EXPECT_FALSE(CRS<long long>({{1LL, 2LL}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<long long>({{1LL, 2LL}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<long long>({{1LL, 2LL}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<long long>({{1LL, 2LL}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<long long>({{1LL, 2LL}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<long>({{1L, 2L}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<long>({{1L, 2L}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<long>({{1L, 2L}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<long>({{1L, 2L}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<long>({{1L, 2L}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<int>({{1, -1}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<int>({{1, -1}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<int>({{1, -1}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<int>({{1, -1}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<int>({{1, -1}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<short>({{1, -1}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<short>({{1, -1}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<short>({{1, -1}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<short>({{1, -1}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<short>({{1, -1}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<char>({{1, -1}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<char>({{1, -1}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<char>({{1, -1}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<char>({{1, -1}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<char>({{1, -1}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<long long>({{2LL}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<long long>({{2LL}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<long long>({{2LL}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<long long>({{2LL}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<long long>({{2LL}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<long>({{2L}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<long>({{2L}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<long>({{2L}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<long>({{2L}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<long>({{2L}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<int>({{2}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<int>({{2}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<int>({{2}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<int>({{2}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<int>({{2}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<short>({{2}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<short>({{2}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<short>({{2}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<short>({{2}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<short>({{2}}) == CRS<char     >({{1  }}));
   
   EXPECT_FALSE(CRS<char>({{2}}) == CRS<long long>({{1LL}}));
   EXPECT_FALSE(CRS<char>({{2}}) == CRS<long     >({{1L }}));
   EXPECT_FALSE(CRS<char>({{2}}) == CRS<int      >({{1  }}));
   EXPECT_FALSE(CRS<char>({{2}}) == CRS<short    >({{1  }}));
   EXPECT_FALSE(CRS<char>({{2}}) == CRS<char     >({{1  }}));
}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_COMPRESSED_ROW_STORAGE_HPP_ */
