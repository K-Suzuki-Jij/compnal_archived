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
//  Created by Kohei Suzuki on 2022/02/13.
//

#ifndef COMPNAL_TEST_HALF_INT_HPP_
#define COMPNAL_TEST_HALF_INT_HPP_

#include "../../src/type/half_int.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(TypeHalfInt, Addition) {
   using type::HalfInt;
   EXPECT_EQ(HalfInt{0.5} + static_cast<HalfInt>    (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<long double>(1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<double>     (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<float>      (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<long long>  (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<long>       (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<int>        (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<short>      (1.0), 1.5);
   EXPECT_EQ(HalfInt{0.5} + static_cast<char>       (1.0), 1.5);
   
   EXPECT_EQ(static_cast<HalfInt>    (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<long double>(1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<double>     (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<float>      (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<long long>  (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<long>       (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<int>        (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<short>      (1.0) + HalfInt{0.5}, 1.5);
   EXPECT_EQ(static_cast<char>       (1.0) + HalfInt{0.5}, 1.5);
   
   EXPECT_EQ(HalfInt{0.5} + static_cast<long double>(1.3), 1.8);
   EXPECT_EQ(HalfInt{0.5} + static_cast<double>     (1.3), 1.8);
   EXPECT_EQ(HalfInt{0.5} + static_cast<float>      (1.3), 1.8f);
   
   EXPECT_EQ(static_cast<long double>(1.3) + HalfInt{0.5}, 1.8);
   EXPECT_EQ(static_cast<double>     (1.3) + HalfInt{0.5}, 1.8);
   EXPECT_EQ(static_cast<float>      (1.3) + HalfInt{0.5}, 1.8f);
}

TEST(TypeHalfInt, Subtraction) {
   using type::HalfInt;
   EXPECT_EQ(HalfInt{0.5} - static_cast<HalfInt>    (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<long double>(1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<double>     (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<float>      (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<long long>  (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<long>       (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<int>        (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<short>      (1.0), -0.5);
   EXPECT_EQ(HalfInt{0.5} - static_cast<char>       (1.0), -0.5);
   
   EXPECT_EQ(static_cast<HalfInt>    (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long double>(1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<double>     (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<float>      (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long long>  (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long>       (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<int>        (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<short>      (1.0) - HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<char>       (1.0) - HalfInt{0.5}, 0.5);
   
   EXPECT_EQ(HalfInt{0.5} - static_cast<long double>(1.3), -0.8);
   EXPECT_EQ(HalfInt{0.5} - static_cast<double>     (1.3), -0.8);
   EXPECT_NEAR(HalfInt{0.5} - static_cast<float>    (1.3), -0.8, std::pow(10, -7));
   
   EXPECT_EQ(static_cast<long double>(1.3) - HalfInt{0.5}, 0.8);
   EXPECT_EQ(static_cast<double>     (1.3) - HalfInt{0.5}, 0.8);
   EXPECT_NEAR(static_cast<float>    (1.3) - HalfInt{0.5}, 0.8, std::pow(10, -7));
}

TEST(TypeHalfInt, Multiplication) {
   using type::HalfInt;
   EXPECT_EQ(HalfInt{0.5} * static_cast<HalfInt>    (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<long double>(1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<double>     (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<float>      (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<long long>  (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<long>       (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<int>        (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<short>      (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} * static_cast<char>       (1.0), 0.5);
   
   EXPECT_EQ(static_cast<HalfInt>    (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long double>(1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<double>     (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<float>      (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long long>  (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<long>       (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<int>        (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<short>      (1.0) * HalfInt{0.5}, 0.5);
   EXPECT_EQ(static_cast<char>       (1.0) * HalfInt{0.5}, 0.5);
   
   EXPECT_EQ(HalfInt{0.5} * static_cast<long double>(1.3), 0.65);
   EXPECT_EQ(HalfInt{0.5} * static_cast<double>     (1.3), 0.65);
   EXPECT_NEAR(HalfInt{0.5} * static_cast<float>    (1.3), 0.65, std::pow(10, -7));
   
   EXPECT_EQ(static_cast<long double>(1.3) * HalfInt{0.5}, 0.65);
   EXPECT_EQ(static_cast<double>     (1.3) * HalfInt{0.5}, 0.65);
   EXPECT_NEAR(static_cast<float>    (1.3) * HalfInt{0.5}, 0.65, std::pow(10, -7));
}

TEST(TypeHalfInt, Division) {
   using type::HalfInt;
   EXPECT_EQ(HalfInt{0.5} / static_cast<HalfInt>    (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<long double>(1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<double>     (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<float>      (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<long long>  (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<long>       (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<int>        (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<short>      (1.0), 0.5);
   EXPECT_EQ(HalfInt{0.5} / static_cast<char>       (1.0), 0.5);
   
   EXPECT_EQ(static_cast<HalfInt>    (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<long double>(1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<double>     (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<float>      (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<long long>  (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<long>       (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<int>        (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<short>      (1.0) / HalfInt{0.5}, 2.0);
   EXPECT_EQ(static_cast<char>       (1.0) / HalfInt{0.5}, 2.0);
   
   EXPECT_EQ(  HalfInt{0.5} / static_cast<long double>(1.3), static_cast<long double>(0.5)/1.3);
   EXPECT_EQ(  HalfInt{0.5} / static_cast<double>     (1.3), 0.5/1.3);
   EXPECT_NEAR(HalfInt{0.5} / static_cast<float>      (1.3), 0.5/1.3, std::pow(10, -7));
   
   EXPECT_EQ(static_cast<long double>(1.3) / HalfInt{0.5}, 2.6);
   EXPECT_EQ(static_cast<double>     (1.3) / HalfInt{0.5}, 2.6);
   EXPECT_NEAR(static_cast<float>    (1.3) / HalfInt{0.5}, 2.6, std::pow(10, -7));
}

TEST(TypeHalfInt, Equality) {
   using type::HalfInt;
   EXPECT_TRUE(HalfInt{0.5} == static_cast<HalfInt>    (0.5));
   EXPECT_TRUE(HalfInt{0.5} == static_cast<long double>(0.5));
   EXPECT_TRUE(HalfInt{0.5} == static_cast<double>     (0.5));
   EXPECT_TRUE(HalfInt{0.5} == static_cast<float>      (0.5));
   EXPECT_TRUE(HalfInt{1.0} == static_cast<long long>  (1.0));
   EXPECT_TRUE(HalfInt{1.0} == static_cast<long>       (1.0));
   EXPECT_TRUE(HalfInt{1.0} == static_cast<int>        (1.0));
   EXPECT_TRUE(HalfInt{1.0} == static_cast<short>      (1.0));
   EXPECT_TRUE(HalfInt{1.0} == static_cast<char>       (1.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (0.5) == HalfInt{0.5});
   EXPECT_TRUE(static_cast<long double>(0.5) == HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (0.5) == HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (0.5) == HalfInt{0.5});
   EXPECT_TRUE(static_cast<long long>  (1.0) == HalfInt{1.0});
   EXPECT_TRUE(static_cast<long>       (1.0) == HalfInt{1.0});
   EXPECT_TRUE(static_cast<int>        (1.0) == HalfInt{1.0});
   EXPECT_TRUE(static_cast<short>      (1.0) == HalfInt{1.0});
   EXPECT_TRUE(static_cast<char>       (1.0) == HalfInt{1.0});
   
   EXPECT_FALSE(HalfInt{0.5} == static_cast<long double>(1.3));
   EXPECT_FALSE(HalfInt{0.5} == static_cast<double>     (1.3));
   EXPECT_FALSE(HalfInt{0.5} == static_cast<float>      (1.3));
   
   EXPECT_FALSE(static_cast<long double>(1.3) == HalfInt{0.5});
   EXPECT_FALSE(static_cast<double>     (1.3) == HalfInt{0.5});
   EXPECT_FALSE(static_cast<float>      (1.3) == HalfInt{0.5});
}

TEST(TypeHalfInt, Inequality) {
   using type::HalfInt;
   EXPECT_FALSE(HalfInt{0.5} != static_cast<HalfInt>    (0.5));
   EXPECT_FALSE(HalfInt{0.5} != static_cast<long double>(0.5));
   EXPECT_FALSE(HalfInt{0.5} != static_cast<double>     (0.5));
   EXPECT_FALSE(HalfInt{0.5} != static_cast<float>      (0.5));
   EXPECT_FALSE(HalfInt{1.0} != static_cast<long long>  (1.0));
   EXPECT_FALSE(HalfInt{1.0} != static_cast<long>       (1.0));
   EXPECT_FALSE(HalfInt{1.0} != static_cast<int>        (1.0));
   EXPECT_FALSE(HalfInt{1.0} != static_cast<short>      (1.0));
   EXPECT_FALSE(HalfInt{1.0} != static_cast<char>       (1.0));
   
   EXPECT_FALSE(static_cast<HalfInt>    (0.5) != HalfInt{0.5});
   EXPECT_FALSE(static_cast<long double>(0.5) != HalfInt{0.5});
   EXPECT_FALSE(static_cast<double>     (0.5) != HalfInt{0.5});
   EXPECT_FALSE(static_cast<float>      (0.5) != HalfInt{0.5});
   EXPECT_FALSE(static_cast<long long>  (1.0) != HalfInt{1.0});
   EXPECT_FALSE(static_cast<long>       (1.0) != HalfInt{1.0});
   EXPECT_FALSE(static_cast<int>        (1.0) != HalfInt{1.0});
   EXPECT_FALSE(static_cast<short>      (1.0) != HalfInt{1.0});
   EXPECT_FALSE(static_cast<char>       (1.0) != HalfInt{1.0});
   
   EXPECT_TRUE(HalfInt{0.5} != static_cast<long double>(1.3));
   EXPECT_TRUE(HalfInt{0.5} != static_cast<double>     (1.3));
   EXPECT_TRUE(HalfInt{0.5} != static_cast<float>      (1.3));
   
   EXPECT_TRUE(static_cast<long double>(1.3) != HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (1.3) != HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (1.3) != HalfInt{0.5});
}

TEST(TypeHalfInt, ComparisonLessThan) {
   using type::HalfInt;
   EXPECT_TRUE(HalfInt{0.5} < static_cast<HalfInt>    (1.5));
   EXPECT_TRUE(HalfInt{0.5} < static_cast<long double>(1.5));
   EXPECT_TRUE(HalfInt{0.5} < static_cast<double>     (1.5));
   EXPECT_TRUE(HalfInt{0.5} < static_cast<float>      (1.5));
   EXPECT_TRUE(HalfInt{1.0} < static_cast<long long>  (2.0));
   EXPECT_TRUE(HalfInt{1.0} < static_cast<long>       (2.0));
   EXPECT_TRUE(HalfInt{1.0} < static_cast<int>        (2.0));
   EXPECT_TRUE(HalfInt{1.0} < static_cast<short>      (2.0));
   EXPECT_TRUE(HalfInt{1.0} < static_cast<char>       (2.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (0.5) < HalfInt{1.5});
   EXPECT_TRUE(static_cast<long double>(0.5) < HalfInt{1.5});
   EXPECT_TRUE(static_cast<double>     (0.5) < HalfInt{1.5});
   EXPECT_TRUE(static_cast<float>      (0.5) < HalfInt{1.5});
   EXPECT_TRUE(static_cast<long long>  (1.0) < HalfInt{2.0});
   EXPECT_TRUE(static_cast<long>       (1.0) < HalfInt{2.0});
   EXPECT_TRUE(static_cast<int>        (1.0) < HalfInt{2.0});
   EXPECT_TRUE(static_cast<short>      (1.0) < HalfInt{2.0});
   EXPECT_TRUE(static_cast<char>       (1.0) < HalfInt{2.0});
   
   EXPECT_FALSE(HalfInt{0.5} < static_cast<HalfInt>    (0.5));
   EXPECT_FALSE(HalfInt{0.5} < static_cast<long double>(0.5));
   EXPECT_FALSE(HalfInt{0.5} < static_cast<double>     (0.5));
   EXPECT_FALSE(HalfInt{0.5} < static_cast<float>      (0.5));
   EXPECT_FALSE(HalfInt{1.0} < static_cast<long long>  (1.0));
   EXPECT_FALSE(HalfInt{1.0} < static_cast<long>       (1.0));
   EXPECT_FALSE(HalfInt{1.0} < static_cast<int>        (1.0));
   EXPECT_FALSE(HalfInt{1.0} < static_cast<short>      (1.0));
   EXPECT_FALSE(HalfInt{1.0} < static_cast<char>       (1.0));
   
   EXPECT_FALSE(static_cast<HalfInt>    (0.5) < HalfInt{0.5});
   EXPECT_FALSE(static_cast<long double>(0.5) < HalfInt{0.5});
   EXPECT_FALSE(static_cast<double>     (0.5) < HalfInt{0.5});
   EXPECT_FALSE(static_cast<float>      (0.5) < HalfInt{0.5});
   EXPECT_FALSE(static_cast<long long>  (1.0) < HalfInt{1.0});
   EXPECT_FALSE(static_cast<long>       (1.0) < HalfInt{1.0});
   EXPECT_FALSE(static_cast<int>        (1.0) < HalfInt{1.0});
   EXPECT_FALSE(static_cast<short>      (1.0) < HalfInt{1.0});
   EXPECT_FALSE(static_cast<char>       (1.0) < HalfInt{1.0});
}

TEST(TypeHalfInt, ComparisonGreaterThan) {
   using type::HalfInt;
   EXPECT_TRUE(HalfInt{1.5} > static_cast<HalfInt>    (0.5));
   EXPECT_TRUE(HalfInt{1.5} > static_cast<long double>(0.5));
   EXPECT_TRUE(HalfInt{1.5} > static_cast<double>     (0.5));
   EXPECT_TRUE(HalfInt{1.5} > static_cast<float>      (0.5));
   EXPECT_TRUE(HalfInt{2.0} > static_cast<long long>  (1.0));
   EXPECT_TRUE(HalfInt{2.0} > static_cast<long>       (1.0));
   EXPECT_TRUE(HalfInt{2.0} > static_cast<int>        (1.0));
   EXPECT_TRUE(HalfInt{2.0} > static_cast<short>      (1.0));
   EXPECT_TRUE(HalfInt{2.0} > static_cast<char>       (1.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (1.5) > HalfInt{0.5});
   EXPECT_TRUE(static_cast<long double>(1.5) > HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (1.5) > HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (1.5) > HalfInt{0.5});
   EXPECT_TRUE(static_cast<long long>  (2.0) > HalfInt{1.0});
   EXPECT_TRUE(static_cast<long>       (2.0) > HalfInt{1.0});
   EXPECT_TRUE(static_cast<int>        (2.0) > HalfInt{1.0});
   EXPECT_TRUE(static_cast<short>      (2.0) > HalfInt{1.0});
   EXPECT_TRUE(static_cast<char>       (2.0) > HalfInt{1.0});
   
   EXPECT_FALSE(HalfInt{0.5} > static_cast<HalfInt>    (0.5));
   EXPECT_FALSE(HalfInt{0.5} > static_cast<long double>(0.5));
   EXPECT_FALSE(HalfInt{0.5} > static_cast<double>     (0.5));
   EXPECT_FALSE(HalfInt{0.5} > static_cast<float>      (0.5));
   EXPECT_FALSE(HalfInt{1.0} > static_cast<long long>  (1.0));
   EXPECT_FALSE(HalfInt{1.0} > static_cast<long>       (1.0));
   EXPECT_FALSE(HalfInt{1.0} > static_cast<int>        (1.0));
   EXPECT_FALSE(HalfInt{1.0} > static_cast<short>      (1.0));
   EXPECT_FALSE(HalfInt{1.0} > static_cast<char>       (1.0));
   
   EXPECT_FALSE(static_cast<HalfInt>    (0.5) > HalfInt{0.5});
   EXPECT_FALSE(static_cast<long double>(0.5) > HalfInt{0.5});
   EXPECT_FALSE(static_cast<double>     (0.5) > HalfInt{0.5});
   EXPECT_FALSE(static_cast<float>      (0.5) > HalfInt{0.5});
   EXPECT_FALSE(static_cast<long long>  (1.0) > HalfInt{1.0});
   EXPECT_FALSE(static_cast<long>       (1.0) > HalfInt{1.0});
   EXPECT_FALSE(static_cast<int>        (1.0) > HalfInt{1.0});
   EXPECT_FALSE(static_cast<short>      (1.0) > HalfInt{1.0});
   EXPECT_FALSE(static_cast<char>       (1.0) > HalfInt{1.0});
}

TEST(TypeHalfInt, ComparisonLessThanOrEqual) {
   using type::HalfInt;
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<HalfInt>    (1.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<long double>(1.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<double>     (1.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<float>      (1.5));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<long long>  (2.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<long>       (2.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<int>        (2.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<short>      (2.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<char>       (2.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (0.5) <= HalfInt{1.5});
   EXPECT_TRUE(static_cast<long double>(0.5) <= HalfInt{1.5});
   EXPECT_TRUE(static_cast<double>     (0.5) <= HalfInt{1.5});
   EXPECT_TRUE(static_cast<float>      (0.5) <= HalfInt{1.5});
   EXPECT_TRUE(static_cast<long long>  (1.0) <= HalfInt{2.0});
   EXPECT_TRUE(static_cast<long>       (1.0) <= HalfInt{2.0});
   EXPECT_TRUE(static_cast<int>        (1.0) <= HalfInt{2.0});
   EXPECT_TRUE(static_cast<short>      (1.0) <= HalfInt{2.0});
   EXPECT_TRUE(static_cast<char>       (1.0) <= HalfInt{2.0});
   
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<HalfInt>    (0.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<long double>(0.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<double>     (0.5));
   EXPECT_TRUE(HalfInt{0.5} <= static_cast<float>      (0.5));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<long long>  (1.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<long>       (1.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<int>        (1.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<short>      (1.0));
   EXPECT_TRUE(HalfInt{1.0} <= static_cast<char>       (1.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (0.5) <= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long double>(0.5) <= HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (0.5) <= HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (0.5) <= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long long>  (1.0) <= HalfInt{1.0});
   EXPECT_TRUE(static_cast<long>       (1.0) <= HalfInt{1.0});
   EXPECT_TRUE(static_cast<int>        (1.0) <= HalfInt{1.0});
   EXPECT_TRUE(static_cast<short>      (1.0) <= HalfInt{1.0});
   EXPECT_TRUE(static_cast<char>       (1.0) <= HalfInt{1.0});
}

TEST(TypeHalfInt, ComparisonGreaterThanOrEqual) {
   using type::HalfInt;
   EXPECT_TRUE(HalfInt{1.5} >= static_cast<HalfInt>    (0.5));
   EXPECT_TRUE(HalfInt{1.5} >= static_cast<long double>(0.5));
   EXPECT_TRUE(HalfInt{1.5} >= static_cast<double>     (0.5));
   EXPECT_TRUE(HalfInt{1.5} >= static_cast<float>      (0.5));
   EXPECT_TRUE(HalfInt{2.0} >= static_cast<long long>  (1.0));
   EXPECT_TRUE(HalfInt{2.0} >= static_cast<long>       (1.0));
   EXPECT_TRUE(HalfInt{2.0} >= static_cast<int>        (1.0));
   EXPECT_TRUE(HalfInt{2.0} >= static_cast<short>      (1.0));
   EXPECT_TRUE(HalfInt{2.0} >= static_cast<char>       (1.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (1.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long double>(1.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (1.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (1.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long long>  (2.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<long>       (2.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<int>        (2.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<short>      (2.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<char>       (2.0) >= HalfInt{1.0});
   
   EXPECT_TRUE(HalfInt{0.5} >= static_cast<HalfInt>    (0.5));
   EXPECT_TRUE(HalfInt{0.5} >= static_cast<long double>(0.5));
   EXPECT_TRUE(HalfInt{0.5} >= static_cast<double>     (0.5));
   EXPECT_TRUE(HalfInt{0.5} >= static_cast<float>      (0.5));
   EXPECT_TRUE(HalfInt{1.0} >= static_cast<long long>  (1.0));
   EXPECT_TRUE(HalfInt{1.0} >= static_cast<long>       (1.0));
   EXPECT_TRUE(HalfInt{1.0} >= static_cast<int>        (1.0));
   EXPECT_TRUE(HalfInt{1.0} >= static_cast<short>      (1.0));
   EXPECT_TRUE(HalfInt{1.0} >= static_cast<char>       (1.0));
   
   EXPECT_TRUE(static_cast<HalfInt>    (0.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long double>(0.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<double>     (0.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<float>      (0.5) >= HalfInt{0.5});
   EXPECT_TRUE(static_cast<long long>  (1.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<long>       (1.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<int>        (1.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<short>      (1.0) >= HalfInt{1.0});
   EXPECT_TRUE(static_cast<char>       (1.0) >= HalfInt{1.0});
}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_HALF_INT_HPP_ */
