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
//  Created by Kohei Suzuki on 2022/01/08.
//

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_1D_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_1D_HPP_

#include "../../src/model/base_u1_spin_1d.hpp"
#include "../include/all.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(ModelBaseU1Spin1D, ConstructorDefault) {
   model::BaseU1Spin_1D<double> model;
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 0);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{});
   ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{});
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSize) {
   model::BaseU1Spin_1D<double> model(10);
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{});
   ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{});
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSizeSpin) {
   model::BaseU1Spin_1D<double> model(10, 1);
   TestSpinOne(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{});
   ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{});
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSizeSpinTotalSz) {
   model::BaseU1Spin_1D<double> model(10, 1, 1);
   TestSpinOne(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 1);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{});
   ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{});
}

TEST(ModelBaseU1Spin1D, SetSystemSize) {
   model::BaseU1Spin_1D<double> model;
   model.SetSystemSize(5);
   EXPECT_EQ(model.GetSystemSize(), 5);
}

TEST(ModelBaseU1Spin1D, SetTotalSz) {
   model::BaseU1Spin_1D<double> model;
   model.SetTotalSz(2);
   EXPECT_EQ(model.GetTotalSz(), 2);
   EXPECT_THROW(model.SetTotalSz(1.9), std::runtime_error);
}

TEST(ModelBaseU1Spin1D, SetMagnitudeSpin) {
   model::BaseU1Spin_1D<double> model;
   model.SetMagnitudeSpin(1.5);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.5);
   EXPECT_THROW(model.SetMagnitudeSpin(1.3), std::runtime_error);
}

TEST(ModelBaseU1Spin1D, isValidQNumber) {
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.0));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.0));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.0));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -4.5));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.5));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.5));

   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +10.0));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +0.0 ));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, -10.0));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -9.0 ));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -1.0 ));
   EXPECT_TRUE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, +9.0 ));
   
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.5));
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.5));
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.5));
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -5.0));
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.0));
   EXPECT_FALSE(model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.0));
}

TEST(ModelBaseU1Spin1D, CalculateTargetDim) {
   //Spin-1/2
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(0, 0.5, +0.0), 0);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(1, 0.5, +0.5), 1);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(2, 0.5, +0.0), 2);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(3, 0.5, +0.5), 3);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, +0.0), 6);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, +2.0), 1);
   EXPECT_EQ(model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, -2.0), 1);
}

}
}

#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_1D_HPP_ */
