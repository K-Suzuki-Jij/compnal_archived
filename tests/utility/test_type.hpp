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
//  Created by Kohei Suzuki on 2022/01/30.
//

#ifndef COMPNAL_TEST_UTILITY_TYPE_HPP_
#define COMPNAL_TEST_UTILITY_TYPE_HPP_

#include "../../src/type.hpp"
#include "../include/all.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(HalfInt, Arithmetic) {
   EXPECT_EQ(0.5_hi + 0.5_hi, 1.0_hi);
   EXPECT_EQ(0.5_hi + 0.5_hi, 1.0   );
   EXPECT_EQ(0.5_hi + 0.5_hi, 1     );

   EXPECT_EQ(0.5_hi + 0.5, 1.0_hi);
   EXPECT_EQ(0.5_hi + 0.5, 1.0   );
   EXPECT_EQ(0.5_hi + 0.5, 1     );
   
   EXPECT_EQ(0.5_hi + 0.75, 1.25);
   
   
   EXPECT_EQ(0.5_hi - 1.5_hi, -1.0_hi);
   EXPECT_EQ(0.5_hi - 1.5_hi, -1.0   );
   EXPECT_EQ(0.5_hi - 1.5_hi, -1     );

   EXPECT_EQ(0.5_hi - 1.5, -1.0_hi);
   EXPECT_EQ(0.5_hi - 1.5, -1.0   );
   EXPECT_EQ(0.5_hi - 1.5, -1     );
   
   EXPECT_EQ(0.5_hi - 0.75, -0.25);
   
   EXPECT_EQ(0.5_hi*1.5_hi, 0.75);
   EXPECT_EQ(0.5_hi*1.5, 0.75);
   EXPECT_EQ(0.5_hi*1.0_hi, 0.5_hi);
   EXPECT_EQ(0.5_hi*1.0_hi, 0.5);
   EXPECT_EQ(0.5_hi*1.0, 0.5_hi);
   EXPECT_EQ(0.5_hi*1.0, 0.5);
   EXPECT_EQ(0.5_hi*1, 0.5_hi);
   EXPECT_EQ(0.5_hi*1, 0.5);

   
}

} // namespace test
} // namespace compnal


#endif /* COMPNAL_TEST_UTILITY_TYPE_HPP_ */
