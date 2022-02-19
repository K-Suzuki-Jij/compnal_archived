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
//  Created by Kohei Suzuki on 2022/02/15.
//

#ifndef COMPNAL_TEST_BRAKET_VECTOR_HPP_
#define COMPNAL_TEST_BRAKET_VECTOR_HPP_

#include "../../src/type/braket_vector.hpp"
#include <cmath>
#include <type_traits>
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(BraketVector, Fill) {
   
}

TEST(BraketVector, Free) {
   
}

TEST(BraketVector, Clear) {
   
}

TEST(BraketVector, Assign) {
   
}

TEST(BraketVector, Normalize) {
   
}

TEST(BraketVector, MultiplyByScalar) {
   
}

TEST(BraketVector, L1Norm) {
   
}

TEST(BraketVector, L2Norm) {
   
}

TEST(BraketVector, CalculateL1Norm) {
   
}

TEST(BraketVector, Addition) {
   auto check = [](const auto &v1, const auto &v2, const auto &prec, const auto &line) {
      EXPECT_EQ(v1.value_list.size(), v2.value_list.size());
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(v1)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(v2)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(v1, v2);
      }
      else {
         for (std::int64_t i = 0; i < v1.value_list.size(); ++i) {
            const auto lhs = std::abs(v1.value_list[i] - v2.value_list[i]);
            const auto rhs = 4*prec;
            EXPECT_TRUE(lhs <= rhs) << lhs << " !<= " << rhs << " at " << line;
         }
      }
   };
   using type::BraketVector;
   check(BraketVector<long double>({1.3L}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.7L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.7 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.7F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<long long  >({1LL}) , BraketVector<long double>({2.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<long       >({1L }) , BraketVector<long double>({2.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<int        >({1  }) , BraketVector<long double>({2.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<short      >({1  }) , BraketVector<long double>({2.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) + BraketVector<char       >({1  }) , BraketVector<long double>({2.3L}), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<double>({1.3}) + BraketVector<long double>({1.4L}), BraketVector<double>({2.7 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<double     >({1.4 }), BraketVector<double>({2.7 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<float      >({1.4F}), BraketVector<float >({2.7F}), FLT_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<long long  >({1LL}) , BraketVector<double>({2.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<long       >({1L }) , BraketVector<double>({2.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<int        >({1  }) , BraketVector<double>({2.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<short      >({1  }) , BraketVector<double>({2.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) + BraketVector<char       >({1  }) , BraketVector<double>({2.3 }), DBL_EPSILON, __LINE__);
   
   check(BraketVector<float>({1.3F}) + BraketVector<long double>({1.4L}), BraketVector<float>({2.7F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<double     >({1.4 }), BraketVector<float>({2.7F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<float      >({1.4F}), BraketVector<float>({2.7F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<long long  >({1LL}) , BraketVector<float>({2.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<long       >({1L }) , BraketVector<float>({2.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<int        >({1  }) , BraketVector<float>({2.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<short      >({1  }) , BraketVector<float>({2.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) + BraketVector<char       >({1  }) , BraketVector<float>({2.3F}), FLT_EPSILON, __LINE__);
   
   check(BraketVector<long long>({1LL}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<long long  >({1LL}) , BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<long       >({1L }) , BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<int        >({1  }) , BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<short      >({1  }) , BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) + BraketVector<char       >({1  }) , BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);

   check(BraketVector<long>({1L}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<long long  >({1LL}) , BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<long       >({1L }) , BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<int        >({1  }) , BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<short      >({1  }) , BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) + BraketVector<char       >({1  }) , BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<int>({1}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<int>({1}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<int>({1}) + BraketVector<long long  >({1LL}) , BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) + BraketVector<long       >({1L }) , BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) + BraketVector<int        >({1  }) , BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) + BraketVector<short      >({1  }) , BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) + BraketVector<char       >({1  }) , BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<short>({1}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<short>({1}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<short>({1}) + BraketVector<long long  >({1LL}) , BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) + BraketVector<long       >({1L }) , BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) + BraketVector<int        >({1  }) , BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) + BraketVector<short      >({1  }) , BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) + BraketVector<char       >({1  }) , BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<char>({1}) + BraketVector<long double>({1.4L}), BraketVector<long double>({2.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) + BraketVector<double     >({1.4 }), BraketVector<double     >({2.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<char>({1}) + BraketVector<float      >({1.4F}), BraketVector<float      >({2.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<char>({1}) + BraketVector<long long  >({1LL}) , BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) + BraketVector<long       >({1L }) , BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) + BraketVector<int        >({1  }) , BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) + BraketVector<short      >({1  }) , BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) + BraketVector<char       >({1  }) , BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);

}

TEST(BraketVector, Subtraction) {
   auto check = [](const auto &v1, const auto &v2, const auto &prec, const auto &line) {
      EXPECT_EQ(v1.value_list.size(), v2.value_list.size());
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(v1)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(v2)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(v1, v2);
      }
      else {
         for (std::int64_t i = 0; i < v1.value_list.size(); ++i) {
            const auto lhs = std::abs(v1.value_list[i] - v2.value_list[i]);
            const auto rhs = 4*prec;
            EXPECT_TRUE(lhs <= rhs) << lhs << " !<= " << rhs << " at " << line;
         }
      }
   };
   using type::BraketVector;
   check(BraketVector<long double>({1.3L}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.1L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.1 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.1F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<long long  >({1LL}) , BraketVector<long double>({+0.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<long       >({1L }) , BraketVector<long double>({+0.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<int        >({1  }) , BraketVector<long double>({+0.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<short      >({1  }) , BraketVector<long double>({+0.3L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) - BraketVector<char       >({1  }) , BraketVector<long double>({+0.3L}), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<double>({1.3}) - BraketVector<long double>({1.4L}), BraketVector<double>({-0.1 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<double     >({1.4 }), BraketVector<double>({-0.1 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<float      >({1.4F}), BraketVector<float >({-0.1F}), FLT_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<long long  >({1LL}) , BraketVector<double>({+0.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<long       >({1L }) , BraketVector<double>({+0.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<int        >({1  }) , BraketVector<double>({+0.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<short      >({1  }) , BraketVector<double>({+0.3 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) - BraketVector<char       >({1  }) , BraketVector<double>({+0.3 }), DBL_EPSILON, __LINE__);
   
   check(BraketVector<float>({1.3F}) - BraketVector<long double>({1.4L}), BraketVector<float>({-0.1F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<double     >({1.4 }), BraketVector<float>({-0.1F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<float      >({1.4F}), BraketVector<float>({-0.1F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<long long  >({1LL}) , BraketVector<float>({+0.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<long       >({1L }) , BraketVector<float>({+0.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<int        >({1  }) , BraketVector<float>({+0.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<short      >({1  }) , BraketVector<float>({+0.3F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) - BraketVector<char       >({1  }) , BraketVector<float>({+0.3F}), FLT_EPSILON, __LINE__);
   
   check(BraketVector<long long>({1LL}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<long long  >({2LL}) , BraketVector<long long  >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<long       >({2L }) , BraketVector<long       >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<int        >({2  }) , BraketVector<int        >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<short      >({2  }) , BraketVector<short      >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) - BraketVector<char       >({2  }) , BraketVector<char       >({-1LL }), LDBL_EPSILON, __LINE__);

   check(BraketVector<long>({1L}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<long long  >({2LL}) , BraketVector<long long  >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<long       >({2L }) , BraketVector<long       >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<int        >({2  }) , BraketVector<int        >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<short      >({2  }) , BraketVector<short      >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) - BraketVector<char       >({2  }) , BraketVector<char       >({-1LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<int>({1}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<int>({1}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<int>({1}) - BraketVector<long long  >({2LL}) , BraketVector<long long  >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) - BraketVector<long       >({2L }) , BraketVector<long       >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) - BraketVector<int        >({2  }) , BraketVector<int        >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) - BraketVector<short      >({2  }) , BraketVector<short      >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) - BraketVector<char       >({2  }) , BraketVector<char       >({-1LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<short>({1}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<short>({1}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<short>({1}) - BraketVector<long long  >({2LL}) , BraketVector<long long  >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) - BraketVector<long       >({2L }) , BraketVector<long       >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) - BraketVector<int        >({2  }) , BraketVector<int        >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) - BraketVector<short      >({2  }) , BraketVector<short      >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) - BraketVector<char       >({2  }) , BraketVector<char       >({-1LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<char>({1}) - BraketVector<long double>({1.4L}), BraketVector<long double>({-0.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) - BraketVector<double     >({1.4 }), BraketVector<double     >({-0.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<char>({1}) - BraketVector<float      >({1.4F}), BraketVector<float      >({-0.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<char>({1}) - BraketVector<long long  >({2LL}) , BraketVector<long long  >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) - BraketVector<long       >({2L }) , BraketVector<long       >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) - BraketVector<int        >({2  }) , BraketVector<int        >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) - BraketVector<short      >({2  }) , BraketVector<short      >({-1LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) - BraketVector<char       >({2  }) , BraketVector<char       >({-1LL }), LDBL_EPSILON, __LINE__);
}

TEST(BraketVector, MultiplicationVectorVector) {
   auto check = [](const auto &v1, const auto &v2, const auto &prec, const auto &line) {
      const bool c1 = std::is_integral<decltype(v1)>::value;
      const bool c2 = std::is_integral<decltype(v2)>::value;
      if (c1 && c2) {
         EXPECT_EQ(v1, v2);
      }
      else {
         const auto lhs = std::abs(v1 - v2);
         const auto rhs = 4*prec;
         EXPECT_TRUE(lhs <= rhs) << lhs << " !<= " << rhs << " at " << line;
      }
   };
   
   using type::BraketVector;
   check(BraketVector<long double>({1.3L}) * BraketVector<long double>({1.4L}), 1.82L, LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<double     >({1.4 }), 1.82 , DBL_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<float      >({1.4F}), 1.82F, FLT_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<long long  >({2LL}) , 2.6L , LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<long       >({2L }) , 2.6L , LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<int        >({2  }) , 2.6L , LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<short      >({2  }) , 2.6L , LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * BraketVector<char       >({2  }) , 2.6L , LDBL_EPSILON, __LINE__);
   
   check(BraketVector<double>({1.3}) * BraketVector<long double>({1.4L}), 1.82 , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<double     >({1.4 }), 1.82 , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<float      >({1.4F}), 1.82F, FLT_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<long long  >({2LL}) , 2.6  , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<long       >({2L }) , 2.6  , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<int        >({2  }) , 2.6  , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<short      >({2  }) , 2.6  , DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * BraketVector<char       >({2  }) , 2.6  , DBL_EPSILON, __LINE__);
   
   check(BraketVector<float>({1.3F}) * BraketVector<long double>({1.4L}), 1.82F, FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<double     >({1.4 }), 1.82F, FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<float      >({1.4F}), 1.82F, FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<long long  >({2LL}) , 2.6F , FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<long       >({2L }) , 2.6F , FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<int        >({2  }) , 2.6F , FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<short      >({2  }) , 2.6F , FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * BraketVector<char       >({2  }) , 2.6F , FLT_EPSILON, __LINE__);
   
   check(BraketVector<long long>({1LL}) * BraketVector<long double>({1.4L}), 1.4L, LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<double     >({1.4 }), 1.4 , DBL_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<float      >({1.4F}), 1.4F, FLT_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<long long  >({2LL}) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<long       >({2L }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<int        >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<short      >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * BraketVector<char       >({2  }) , 2LL , LDBL_EPSILON, __LINE__);

   check(BraketVector<long>({1L}) * BraketVector<long double>({1.4L}), 1.4L, LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<double     >({1.4 }), 1.4 , DBL_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<float      >({1.4F}), 1.4F, FLT_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<long long  >({2LL}) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<long       >({2L }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<int        >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<short      >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * BraketVector<char       >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   
   check(BraketVector<int>({1}) * BraketVector<long double>({1.4L}), 1.4L, LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * BraketVector<double     >({1.4 }), 1.4 , DBL_EPSILON , __LINE__);
   check(BraketVector<int>({1}) * BraketVector<float      >({1.4F}), 1.4F, FLT_EPSILON , __LINE__);
   check(BraketVector<int>({1}) * BraketVector<long long  >({2LL}) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * BraketVector<long       >({2L }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * BraketVector<int        >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * BraketVector<short      >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * BraketVector<char       >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   
   check(BraketVector<short>({1}) * BraketVector<long double>({1.4L}), 1.4L, LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * BraketVector<double     >({1.4 }), 1.4 , DBL_EPSILON , __LINE__);
   check(BraketVector<short>({1}) * BraketVector<float      >({1.4F}), 1.4F, FLT_EPSILON , __LINE__);
   check(BraketVector<short>({1}) * BraketVector<long long  >({2LL}) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * BraketVector<long       >({2L }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * BraketVector<int        >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * BraketVector<short      >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * BraketVector<char       >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   
   check(BraketVector<char>({1}) * BraketVector<long double>({1.4L}), 1.4L, LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * BraketVector<double     >({1.4 }), 1.4 , DBL_EPSILON , __LINE__);
   check(BraketVector<char>({1}) * BraketVector<float      >({1.4F}), 1.4F, FLT_EPSILON , __LINE__);
   check(BraketVector<char>({1}) * BraketVector<long long  >({2LL}) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * BraketVector<long       >({2L }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * BraketVector<int        >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * BraketVector<short      >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * BraketVector<char       >({2  }) , 2LL , LDBL_EPSILON, __LINE__);
}

TEST(BraketVector, MultiplicationVectorScalar) {
   auto check = [](const auto &v1, const auto &v2, const auto &prec, const auto &line) {
      EXPECT_EQ(v1.value_list.size(), v2.value_list.size());
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(v1)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(v2)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(v1, v2);
      }
      else {
         for (std::int64_t i = 0; i < v1.value_list.size(); ++i) {
            const auto lhs = std::abs(v1.value_list[i] - v2.value_list[i]);
            const auto rhs = 4*prec;
            EXPECT_TRUE(lhs <= rhs) << lhs << " !<= " << rhs << " at " << line;
         }
      }
   };
   using type::BraketVector;
   check(BraketVector<long double>({1.3L}) * static_cast<long double>(1.4L), BraketVector<long double>({1.82L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.82 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<float      >(1.4F), BraketVector<float      >({1.82F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<long long  >(1LL ), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<long       >(1L  ), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<int        >(1   ), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<short      >(1   ), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long double>({1.3L}) * static_cast<char       >(1   ), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<double>({1.3}) * static_cast<long double>(1.4L), BraketVector<double>({1.82 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<double     >(1.4 ), BraketVector<double>({1.82 }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<float      >(1.4F), BraketVector<float >({1.82F}), FLT_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<long long  >(1LL ), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<long       >(1L  ), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<int        >(1   ), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<short      >(1   ), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(BraketVector<double>({1.3}) * static_cast<char       >(1   ), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   
   check(BraketVector<float>({1.3F}) * static_cast<long double>(1.4L), BraketVector<double>({1.82F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<double     >(1.4 ), BraketVector<double>({1.82F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<float      >(1.4F), BraketVector<float >({1.82F}), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<long long  >(1LL ), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<long       >(1L  ), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<int        >(1   ), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<short      >(1   ), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(BraketVector<float>({1.3F}) * static_cast<char       >(1   ), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   
   check(BraketVector<long long>({1LL}) * static_cast<long double>(1.4L), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<float      >(1.4F), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<long long  >(2LL ), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<long       >(2L  ), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<int        >(2   ), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<short      >(2   ), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long long>({1LL}) * static_cast<char       >(2   ), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);

   check(BraketVector<long>({1L}) * static_cast<long double>(1.4L), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) * static_cast<float      >(1.4F), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<long>({1L}) * static_cast<long long  >(2LL ), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * static_cast<long       >(2L  ), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * static_cast<int        >(2   ), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * static_cast<short      >(2   ), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<long>({1L}) * static_cast<char       >(2   ), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<int>({1}) * static_cast<long double>(1.4L), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<int>({1}) * static_cast<float      >(1.4F), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<int>({1}) * static_cast<long long  >(2LL ), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * static_cast<long       >(2L  ), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * static_cast<int        >(2   ), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * static_cast<short      >(2   ), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<int>({1}) * static_cast<char       >(2   ), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<short>({1}) * static_cast<long double>(1.4L), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<short>({1}) * static_cast<float      >(1.4F), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<short>({1}) * static_cast<long long  >(2LL ), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * static_cast<long       >(2L  ), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * static_cast<int        >(2   ), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * static_cast<short      >(2   ), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<short>({1}) * static_cast<char       >(2   ), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(BraketVector<char>({1}) * static_cast<long double>(1.4L), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * static_cast<double     >(1.4 ), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(BraketVector<char>({1}) * static_cast<float      >(1.4F), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(BraketVector<char>({1}) * static_cast<long long  >(2LL ), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * static_cast<long       >(2L  ), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * static_cast<int        >(2   ), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * static_cast<short      >(2   ), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(BraketVector<char>({1}) * static_cast<char       >(2   ), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
}

TEST(BraketVector, MultiplicationScalarVector) {
   auto check = [](const auto &v1, const auto &v2, const auto &prec, const auto &line) {
      EXPECT_EQ(v1.value_list.size(), v2.value_list.size());
      const bool c1 = std::is_integral<typename std::remove_cvref<decltype(v1)>::type::ValueType>::value;
      const bool c2 = std::is_integral<typename std::remove_cvref<decltype(v2)>::type::ValueType>::value;
      if (c1 && c2) {
         EXPECT_EQ(v1, v2);
      }
      else {
         for (std::int64_t i = 0; i < v1.value_list.size(); ++i) {
            const auto lhs = std::abs(v1.value_list[i] - v2.value_list[i]);
            const auto rhs = 4*prec;
            EXPECT_TRUE(lhs <= rhs) << lhs << " !<= " << rhs << " at " << line;
         }
      }
   };
   using type::BraketVector;
   check(static_cast<long double>(1.4L) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.82L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<long double>({1.3L}), BraketVector<double     >({1.82 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<long double>({1.3L}), BraketVector<float      >({1.82F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(1LL ) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(1L  ) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(1   ) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(1   ) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(1   ) * BraketVector<long double>({1.3L}), BraketVector<long double>({1.3L }), LDBL_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<double>({1.3}), BraketVector<double>({1.82 }), DBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<double>({1.3}), BraketVector<double>({1.82 }), DBL_EPSILON, __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<double>({1.3}), BraketVector<float >({1.82F}), FLT_EPSILON, __LINE__);
   check(static_cast<long long  >(1LL ) * BraketVector<double>({1.3}), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(static_cast<long       >(1L  ) * BraketVector<double>({1.3}), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(static_cast<int        >(1   ) * BraketVector<double>({1.3}), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(static_cast<short      >(1   ) * BraketVector<double>({1.3}), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   check(static_cast<char       >(1   ) * BraketVector<double>({1.3}), BraketVector<double>({1.3  }), DBL_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<float>({1.3F}), BraketVector<double>({1.82F}), FLT_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<float>({1.3F}), BraketVector<double>({1.82F}), FLT_EPSILON, __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<float>({1.3F}), BraketVector<float >({1.82F}), FLT_EPSILON, __LINE__);
   check(static_cast<long long  >(1LL ) * BraketVector<float>({1.3F}), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(static_cast<long       >(1L  ) * BraketVector<float>({1.3F}), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(static_cast<int        >(1   ) * BraketVector<float>({1.3F}), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(static_cast<short      >(1   ) * BraketVector<float>({1.3F}), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   check(static_cast<char       >(1   ) * BraketVector<float>({1.3F}), BraketVector<double>({1.3F }), FLT_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<long long>({1LL}), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<long long>({1LL}), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<long long>({1LL}), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(2LL ) * BraketVector<long long>({1LL}), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(2L  ) * BraketVector<long long>({1LL}), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(2   ) * BraketVector<long long>({1LL}), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(2   ) * BraketVector<long long>({1LL}), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(2   ) * BraketVector<long long>({1LL}), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);

   check(static_cast<long double>(1.4L) * BraketVector<long>({1L}), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<long>({1L}), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<long>({1L}), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(2LL ) * BraketVector<long>({1L}), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(2L  ) * BraketVector<long>({1L}), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(2   ) * BraketVector<long>({1L}), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(2   ) * BraketVector<long>({1L}), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(2   ) * BraketVector<long>({1L}), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<int>({1}), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<int>({1}), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<int>({1}), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(2LL ) * BraketVector<int>({1}), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(2L  ) * BraketVector<int>({1}), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(2   ) * BraketVector<int>({1}), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(2   ) * BraketVector<int>({1}), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(2   ) * BraketVector<int>({1}), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<short>({1}), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<short>({1}), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<short>({1}), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(2LL ) * BraketVector<short>({1}), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(2L  ) * BraketVector<short>({1}), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(2   ) * BraketVector<short>({1}), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(2   ) * BraketVector<short>({1}), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(2   ) * BraketVector<short>({1}), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
   
   check(static_cast<long double>(1.4L) * BraketVector<char>({1}), BraketVector<long double>({1.4L}), LDBL_EPSILON, __LINE__);
   check(static_cast<double     >(1.4 ) * BraketVector<char>({1}), BraketVector<double     >({1.4 }), DBL_EPSILON , __LINE__);
   check(static_cast<float      >(1.4F) * BraketVector<char>({1}), BraketVector<float      >({1.4F}), FLT_EPSILON , __LINE__);
   check(static_cast<long long  >(2LL ) * BraketVector<char>({1}), BraketVector<long long  >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<long       >(2L  ) * BraketVector<char>({1}), BraketVector<long       >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<int        >(2   ) * BraketVector<char>({1}), BraketVector<int        >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<short      >(2   ) * BraketVector<char>({1}), BraketVector<short      >({2LL }), LDBL_EPSILON, __LINE__);
   check(static_cast<char       >(2   ) * BraketVector<char>({1}), BraketVector<char       >({2LL }), LDBL_EPSILON, __LINE__);
}

TEST(BraketVector, Equality) {
   using type::BraketVector;
   EXPECT_TRUE(BraketVector<long long>({1LL}) == BraketVector<long long>({1LL}));
   EXPECT_TRUE(BraketVector<long long>({1LL}) == BraketVector<long     >({1L }));
   EXPECT_TRUE(BraketVector<long long>({1LL}) == BraketVector<int      >({1  }));
   EXPECT_TRUE(BraketVector<long long>({1LL}) == BraketVector<short    >({1  }));
   EXPECT_TRUE(BraketVector<long long>({1LL}) == BraketVector<char     >({1  }));

   EXPECT_TRUE(BraketVector<long>({1L}) == BraketVector<long long>({1LL}));
   EXPECT_TRUE(BraketVector<long>({1L}) == BraketVector<long     >({1L }));
   EXPECT_TRUE(BraketVector<long>({1L}) == BraketVector<int      >({1  }));
   EXPECT_TRUE(BraketVector<long>({1L}) == BraketVector<short    >({1  }));
   EXPECT_TRUE(BraketVector<long>({1L}) == BraketVector<char     >({1  }));
   
   EXPECT_TRUE(BraketVector<int>({1}) == BraketVector<long long>({1LL}));
   EXPECT_TRUE(BraketVector<int>({1}) == BraketVector<long     >({1L }));
   EXPECT_TRUE(BraketVector<int>({1}) == BraketVector<int      >({1  }));
   EXPECT_TRUE(BraketVector<int>({1}) == BraketVector<short    >({1  }));
   EXPECT_TRUE(BraketVector<int>({1}) == BraketVector<char     >({1  }));
   
   EXPECT_TRUE(BraketVector<short>({1}) == BraketVector<long long>({1LL}));
   EXPECT_TRUE(BraketVector<short>({1}) == BraketVector<long     >({1L }));
   EXPECT_TRUE(BraketVector<short>({1}) == BraketVector<int      >({1  }));
   EXPECT_TRUE(BraketVector<short>({1}) == BraketVector<short    >({1  }));
   EXPECT_TRUE(BraketVector<short>({1}) == BraketVector<char     >({1  }));
   
   EXPECT_TRUE(BraketVector<char>({1}) == BraketVector<long long>({1LL}));
   EXPECT_TRUE(BraketVector<char>({1}) == BraketVector<long     >({1L }));
   EXPECT_TRUE(BraketVector<char>({1}) == BraketVector<int      >({1  }));
   EXPECT_TRUE(BraketVector<char>({1}) == BraketVector<short    >({1  }));
   EXPECT_TRUE(BraketVector<char>({1}) == BraketVector<char     >({1  }));
   
   
   EXPECT_FALSE(BraketVector<long long>({1LL, 0LL}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<long long>({1LL, 0LL}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<long long>({1LL, 0LL}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<long long>({1LL, 0LL}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<long long>({1LL, 0LL}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<long>({1L, 0L}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<long>({1L, 0L}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<long>({1L, 0L}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<long>({1L, 0L}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<long>({1L, 0L}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<int>({1, 0}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<int>({1, 0}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<int>({1, 0}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<int>({1, 0}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<int>({1, 0}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<short>({1, 0}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<short>({1, 0}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<short>({1, 0}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<short>({1, 0}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<short>({1, 0}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<char>({1, 0}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<char>({1, 0}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<char>({1, 0}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<char>({1, 0}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<char>({1, 0}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<long long>({2LL}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<long long>({2LL}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<long long>({2LL}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<long long>({2LL}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<long long>({2LL}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<long>({2L}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<long>({2L}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<long>({2L}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<long>({2L}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<long>({2L}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<int>({2}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<int>({2}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<int>({2}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<int>({2}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<int>({2}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<short>({2}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<short>({2}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<short>({2}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<short>({2}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<short>({2}) == BraketVector<char     >({1  }));
   
   EXPECT_FALSE(BraketVector<char>({2}) == BraketVector<long long>({1LL}));
   EXPECT_FALSE(BraketVector<char>({2}) == BraketVector<long     >({1L }));
   EXPECT_FALSE(BraketVector<char>({2}) == BraketVector<int      >({1  }));
   EXPECT_FALSE(BraketVector<char>({2}) == BraketVector<short    >({1  }));
   EXPECT_FALSE(BraketVector<char>({2}) == BraketVector<char     >({1  }));
}

}
}

#endif /* COMPNAL_TEST_BRAKET_VECTOR_HPP_ */
