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
//  Created by Kohei Suzuki on 2022/02/11.
//

#ifndef COMPNAL_TEST_UTILITY_INTEGER_HPP_
#define COMPNAL_TEST_UTILITY_INTEGER_HPP_

#include <gtest/gtest.h>

#include "../../src/utility/integer.hpp"

namespace compnal {
namespace test {

TEST(UtilityInteger, GenerateIntegerPartition) {
   using VecType = std::vector<std::vector<int>>;

   VecType vec_5_5 = {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}};

   VecType vec_5_4 = {{4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}};

   VecType vec_5_3 = {{3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}};

   VecType vec_5_2 = {{2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}};

   VecType vec_5_1 = {{1, 1, 1, 1, 1}};

   EXPECT_EQ(utility::GenerateIntegerPartition(5, 1), vec_5_1);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 2), vec_5_2);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 3), vec_5_3);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4), vec_5_4);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 5), vec_5_5);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 6), vec_5_5);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 10), vec_5_5);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 100000), vec_5_5);

   EXPECT_EQ(utility::GenerateIntegerPartition(1, 1), VecType{{1}});

   EXPECT_THROW(utility::GenerateIntegerPartition(+5, -2), std::runtime_error);
   EXPECT_THROW(utility::GenerateIntegerPartition(+5, +0), std::runtime_error);
   EXPECT_THROW(utility::GenerateIntegerPartition(-4, -2), std::runtime_error);
   EXPECT_THROW(utility::GenerateIntegerPartition(-4, +2), std::runtime_error);
}

TEST(UtilityInteger, GenerateIntegerPartitionMaxSize) {
   using VecType = std::vector<std::vector<int>>;

   VecType vec_5_4_5 = {{4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}};

   VecType vec_5_4_4 = {{4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}};

   VecType vec_5_4_3 = {{4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}};

   VecType vec_5_4_2 = {{4, 1}, {3, 2}};

   VecType vec_5_4_1 = {};

   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4, 5), vec_5_4_5);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4, 4), vec_5_4_4);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4, 3), vec_5_4_3);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4, 2), vec_5_4_2);
   EXPECT_EQ(utility::GenerateIntegerPartition(5, 4, 1), vec_5_4_1);

   EXPECT_EQ(utility::GenerateIntegerPartition(5, 5, 1), VecType{{5}});
   EXPECT_THROW(utility::GenerateIntegerPartition(5, 4, 0), std::runtime_error);
}

TEST(UtilityInteger, CalculateBinomialCoefficient) {
   EXPECT_THROW(utility::CalculateBinomialCoefficient(9, -1), std::runtime_error);
   EXPECT_THROW(utility::CalculateBinomialCoefficient(-9, 3), std::runtime_error);
   EXPECT_THROW(utility::CalculateBinomialCoefficient(-9, -2), std::runtime_error);
   EXPECT_THROW(utility::CalculateBinomialCoefficient(9, 11), std::runtime_error);

   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 0), 1);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 1), 9);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 2), 36);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 3), 84);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 4), 126);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 5), 126);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 6), 84);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 7), 36);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 8), 9);
   EXPECT_EQ(utility::CalculateBinomialCoefficient(9, 9), 1);
}

TEST(UtilityInteger, GenerateBinomialTable) {
   std::vector<std::vector<std::int64_t>> ans(4);
   ans[0].resize(1);
   ans[0][0] = 1;

   ans[1].resize(2);
   ans[1][0] = 1;
   ans[1][1] = 1;

   ans[2].resize(3);
   ans[2][0] = 1;
   ans[2][1] = 2;
   ans[2][2] = 1;

   ans[3].resize(4);
   ans[3][0] = 1;
   ans[3][1] = 3;
   ans[3][2] = 3;
   ans[3][3] = 1;

   EXPECT_EQ(utility::GenerateBinomialTable(3), ans);
   EXPECT_EQ(utility::GenerateBinomialTable(0), std::vector<std::vector<std::int64_t>>{{1}});
   EXPECT_THROW(utility::GenerateBinomialTable(-1), std::runtime_error);
}

TEST(UtilityInteger, CalculateNumPermutation) {
   EXPECT_EQ(utility::CalculateNumPermutation<int>({1, 2, 3}), 6);
   EXPECT_EQ(utility::CalculateNumPermutation<int>({1, 2, 2, 3}), 12);
   EXPECT_EQ(utility::CalculateNumPermutation<std::string>({"a", "b", "b", "c"}), 12);
}

TEST(UtilityInteger, GenerateNthPermutationInt) {
   std::vector<int> vec_int_p1 = {1, 2, 3};
   std::vector<int> vec_ans_p1_0 = {1, 2, 3};
   std::vector<int> vec_ans_p1_1 = {1, 3, 2};
   std::vector<int> vec_ans_p1_2 = {2, 1, 3};
   std::vector<int> vec_ans_p1_3 = {2, 3, 1};
   std::vector<int> vec_ans_p1_4 = {3, 1, 2};
   std::vector<int> vec_ans_p1_5 = {3, 2, 1};

   EXPECT_THROW(utility::GenerateNthPermutation(vec_int_p1, -1), std::runtime_error);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 0), vec_ans_p1_0);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 1), vec_ans_p1_1);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 2), vec_ans_p1_2);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 3), vec_ans_p1_3);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 4), vec_ans_p1_4);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p1, 5), vec_ans_p1_5);
   EXPECT_THROW(utility::GenerateNthPermutation(vec_int_p1, 6), std::runtime_error);

   std::vector<int> vec_int_p2 = {1, 2, 2, 3};
   std::vector<int> vec_ans_p2_0 = {1, 2, 2, 3};
   std::vector<int> vec_ans_p2_1 = {1, 2, 3, 2};
   std::vector<int> vec_ans_p2_2 = {1, 3, 2, 2};
   std::vector<int> vec_ans_p2_3 = {2, 1, 2, 3};
   std::vector<int> vec_ans_p2_4 = {2, 1, 3, 2};
   std::vector<int> vec_ans_p2_5 = {2, 2, 1, 3};
   std::vector<int> vec_ans_p2_6 = {2, 2, 3, 1};
   std::vector<int> vec_ans_p2_7 = {2, 3, 1, 2};
   std::vector<int> vec_ans_p2_8 = {2, 3, 2, 1};
   std::vector<int> vec_ans_p2_9 = {3, 1, 2, 2};
   std::vector<int> vec_ans_p2_10 = {3, 2, 1, 2};
   std::vector<int> vec_ans_p2_11 = {3, 2, 2, 1};

   EXPECT_THROW(utility::GenerateNthPermutation(vec_int_p2, -1), std::runtime_error);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 0), vec_ans_p2_0);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 1), vec_ans_p2_1);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 2), vec_ans_p2_2);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 3), vec_ans_p2_3);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 4), vec_ans_p2_4);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 5), vec_ans_p2_5);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 6), vec_ans_p2_6);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 7), vec_ans_p2_7);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 8), vec_ans_p2_8);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 9), vec_ans_p2_9);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 10), vec_ans_p2_10);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_int_p2, 11), vec_ans_p2_11);
   EXPECT_THROW(utility::GenerateNthPermutation(vec_int_p2, 12), std::runtime_error);
}

TEST(UtilityInteger, GenerateNthPermutationString) {
   std::vector<std::string> vec_str_p1 = {"a", "b", "c"};
   std::vector<std::string> vec_ans_p1_0 = {"a", "b", "c"};
   std::vector<std::string> vec_ans_p1_1 = {"a", "c", "b"};
   std::vector<std::string> vec_ans_p1_2 = {"b", "a", "c"};
   std::vector<std::string> vec_ans_p1_3 = {"b", "c", "a"};
   std::vector<std::string> vec_ans_p1_4 = {"c", "a", "b"};
   std::vector<std::string> vec_ans_p1_5 = {"c", "b", "a"};

   EXPECT_THROW(utility::GenerateNthPermutation(vec_str_p1, -1), std::runtime_error);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 0), vec_ans_p1_0);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 1), vec_ans_p1_1);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 2), vec_ans_p1_2);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 3), vec_ans_p1_3);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 4), vec_ans_p1_4);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p1, 5), vec_ans_p1_5);
   EXPECT_THROW(utility::GenerateNthPermutation(vec_str_p1, 6), std::runtime_error);

   std::vector<std::string> vec_str_p2 = {"a", "b", "b", "c"};
   std::vector<std::string> vec_ans_p2_0 = {"a", "b", "b", "c"};
   std::vector<std::string> vec_ans_p2_1 = {"a", "b", "c", "b"};
   std::vector<std::string> vec_ans_p2_2 = {"a", "c", "b", "b"};
   std::vector<std::string> vec_ans_p2_3 = {"b", "a", "b", "c"};
   std::vector<std::string> vec_ans_p2_4 = {"b", "a", "c", "b"};
   std::vector<std::string> vec_ans_p2_5 = {"b", "b", "a", "c"};
   std::vector<std::string> vec_ans_p2_6 = {"b", "b", "c", "a"};
   std::vector<std::string> vec_ans_p2_7 = {"b", "c", "a", "b"};
   std::vector<std::string> vec_ans_p2_8 = {"b", "c", "b", "a"};
   std::vector<std::string> vec_ans_p2_9 = {"c", "a", "b", "b"};
   std::vector<std::string> vec_ans_p2_10 = {"c", "b", "a", "b"};
   std::vector<std::string> vec_ans_p2_11 = {"c", "b", "b", "a"};

   EXPECT_THROW(utility::GenerateNthPermutation(vec_str_p2, -1), std::runtime_error);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 0), vec_ans_p2_0);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 1), vec_ans_p2_1);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 2), vec_ans_p2_2);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 3), vec_ans_p2_3);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 4), vec_ans_p2_4);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 5), vec_ans_p2_5);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 6), vec_ans_p2_6);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 7), vec_ans_p2_7);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 8), vec_ans_p2_8);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 9), vec_ans_p2_9);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 10), vec_ans_p2_10);
   EXPECT_EQ(utility::GenerateNthPermutation(vec_str_p2, 11), vec_ans_p2_11);
   EXPECT_THROW(utility::GenerateNthPermutation(vec_str_p2, 12), std::runtime_error);
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_UTILITY_INTEGER_HPP_ */
