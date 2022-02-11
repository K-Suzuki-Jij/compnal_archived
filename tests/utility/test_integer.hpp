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

#include "../../src/utility/integer.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(UtilityInteger, GenerateIntegerPartition) {
   using VecType = std::vector<std::vector<int>>;
   
   VecType vec_5_5 = {
      {5},
      {4, 1},
      {3, 2},
      {3, 1, 1},
      {2, 2, 1},
      {2, 1, 1, 1},
      {1, 1, 1, 1, 1}
   };
   
   VecType vec_5_4 = {
      {4, 1},
      {3, 2},
      {3, 1, 1},
      {2, 2, 1},
      {2, 1, 1, 1},
      {1, 1, 1, 1, 1}
   };
   
   VecType vec_5_3 = {
      {3, 2},
      {3, 1, 1},
      {2, 2, 1},
      {2, 1, 1, 1},
      {1, 1, 1, 1, 1}
   };
   
   VecType vec_5_2 = {
      {2, 2, 1},
      {2, 1, 1, 1},
      {1, 1, 1, 1, 1}
   };
   
   VecType vec_5_1 = {
      {1, 1, 1, 1, 1}
   };
   
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


} //namespace test
} //namespace compnal


#endif /* COMPNAL_TEST_UTILITY_INTEGER_HPP_ */
