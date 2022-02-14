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
   using type::CRS;
   CRS<int> i ({
      {0, 1},
      {1, 0}
   });
   
   CRS<long double> j ({
      {0, 1.0},
      {1.0, 0}
   });
   
   
   EXPECT_EQ(i + j, CRS<long double>({{0, 2.0}, {2.0, 0.0}}));
   
   
   
   
}




} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_COMPRESSED_ROW_STORAGE_HPP_ */
