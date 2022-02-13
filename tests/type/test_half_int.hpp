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

TEST(TypeHalfInt, Basic) {
   type::HalfInt x = 1.5;
   EXPECT_THROW(x += 2.1, std::runtime_error);
}


} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_HALF_INT_HPP_ */
