//
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
//  test_quadratic_any_interaction.hpp
//  compnal
//
//  Created by kohei on 2022/10/18.
//  
//

#ifndef COMPNAL_TEST_INTERACTION_QUADRATIC_ANY_HPP_
#define COMPNAL_TEST_INTERACTION_QUADRATIC_ANY_HPP_

#include "../../src/interaction/quadratic_any_interaction.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(InteractionQuadraticAny, Basic) {
   using QuadClass = interaction::QuadraticAnyInteraction<double>;
   const QuadClass::LinearType linear = {
      {1, 1.0},
      {"a", 1.0},
      {utility::AnyTupleType{2, "b"}, 1.5}
   };
   
   const QuadClass::QuadraticType quadratic = {
      {{1, 1}, +3.0},
      {{1, 2}, -1.0},
      {{"a", 1}, -1.5},
      {{utility::AnyTupleType{2, "b"}, utility::AnyTupleType{2, "a"}}, -1.5}
   };
   
   const auto quad_interaction = QuadClass{linear, quadratic};
   auto index_list = std::vector<QuadClass::IndexType>{1, 2, "a", utility::AnyTupleType{2, "b"}, utility::AnyTupleType{2, "a"}};
   std::sort(index_list.begin(), index_list.end());
   
   EXPECT_EQ(quad_interaction.GetSystemSize(), 5);
   EXPECT_EQ(quad_interaction.GetDegree(), 2);
   EXPECT_EQ(quad_interaction.GetIndexList(), index_list);
   for (const auto &it: quad_interaction.GetIndexMap()) {
      EXPECT_EQ(it.first, index_list[it.second]);
   }
   EXPECT_EQ(quad_interaction.GetConstant(), 3.0);
   for (const auto &it: linear) {
      EXPECT_DOUBLE_EQ(quad_interaction.GetLinear().at(quad_interaction.GetIndexMap().at(it.first)), it.second);
   }
   
   //TODO: Add quadratic CRS tests
   
}

} // namespace test
} // namespace compnal


#endif /* COMPNAL_TEST_INTERACTION_QUADRATIC_ANY_HPP_ */
