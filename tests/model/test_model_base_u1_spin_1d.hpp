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
#include "../test.hpp"
#include <gtest/gtest.h>

namespace {

using compnal::LInt;
using compnal::Map;
using compnal::test::ExpectEQ;
using compnal::test::ExpectNear;
using compnal::model::BaseU1Spin_1D;
using compnal::sparse_matrix::CRS;

}

namespace {

template<typename RealType>
void TestSpinOneHalf(const BaseU1Spin_1D<RealType> &model) {
   
   EXPECT_EQ(model.GetDimOnsite(), 2);
   EXPECT_EQ(model.GetMagnitudeSpin(), 0.5);
   
   const CRS<RealType> ref_sp ({{+0.0, +1.0}, {+0.0, +0.0}});
   const CRS<RealType> ref_sm ({{+0.0, +0.0}, {+1.0, +0.0}});
   const CRS<RealType> ref_sx ({{+0.0, +0.5}, {+0.5, +0.0}});
   const CRS<RealType> ref_isy({{+0.0, +0.5}, {-0.5, +0.0}});
   const CRS<RealType> ref_sz ({{+0.5, +0.0}, {+0.0, -0.5}});
   
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , std::pow(10, -15)));
}

template<typename RealType>
void TestSpinOne(const BaseU1Spin_1D<RealType> &model) {
   
   EXPECT_EQ(model.GetDimOnsite(), 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.0);
   
   const CRS<RealType> ref_sp ({
      {+0.0, +std::sqrt(2), +0.0},
      {+0.0, +0.0, +std::sqrt(2)},
      {+0.0, +0.0, +0.0}
   });
   
   const CRS<RealType> ref_sm ({
      {+0.0, +0.0, +0.0},
      {+std::sqrt(2), +0.0, +0.0},
      {+0.0, +std::sqrt(2), +0.0}
   });
   
   const CRS<RealType> ref_sx ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {+1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, +1.0/std::sqrt(2), +0.0}
   });
   
   const CRS<RealType> ref_isy ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {-1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, -1.0/std::sqrt(2), +0.0}
   });
   
   const CRS<RealType> ref_sz ({
      {+1.0, +0.0, +0.0},
      {+0.0, +0.0, +0.0},
      {+0.0, +0.0, -1.0}
   });
   
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, std::pow(10, -15)));
   EXPECT_TRUE(ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , std::pow(10, -15)));
}

} // namespace

TEST(ModelBaseU1Spin1D, ConstructorDefault) {
   BaseU1Spin_1D<double> model;
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 0);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(ExpectEQ(model.GetBases(), Map<int, std::vector<LInt>>{}));
   EXPECT_TRUE(ExpectEQ(model.GetBasesInv(), Map<int, Map<LInt, LInt>>{}));
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSize) {
   BaseU1Spin_1D<double> model(10);
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(ExpectEQ(model.GetBases(), Map<int, std::vector<LInt>>{}));
   EXPECT_TRUE(ExpectEQ(model.GetBasesInv(), Map<int, Map<LInt, LInt>>{}));
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSizeSpin) {
   BaseU1Spin_1D<double> model(10, 1);
   TestSpinOne(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(ExpectEQ(model.GetBases(), Map<int, std::vector<LInt>>{}));
   EXPECT_TRUE(ExpectEQ(model.GetBasesInv(), Map<int, Map<LInt, LInt>>{}));
}

TEST(ModelBaseU1Spin1D, ConstructorSystemSizeSpinTotalSz) {
   BaseU1Spin_1D<double> model(10, 1, 1);
   TestSpinOne(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 1);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(ExpectEQ(model.GetBases(), Map<int, std::vector<LInt>>{}));
   EXPECT_TRUE(ExpectEQ(model.GetBasesInv(), Map<int, Map<LInt, LInt>>{}));
}

TEST(ModelBaseU1Spin1D, SetSystemSize) {
   BaseU1Spin_1D<double> model;
   model.SetSystemSize(5);
   EXPECT_EQ(model.GetSystemSize(), 5);
}

TEST(ModelBaseU1Spin1D, SetTotalSz) {
   BaseU1Spin_1D<double> model;
   model.SetTotalSz(2);
   EXPECT_EQ(model.GetTotalSz(), 2);
   EXPECT_THROW(model.SetTotalSz(1.9), std::runtime_error);
}

TEST(ModelBaseU1Spin1D, SetMagnitudeSpin) {
   BaseU1Spin_1D<double> model;
   model.SetMagnitudeSpin(1.5);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.5);
   EXPECT_THROW(model.SetMagnitudeSpin(1.3), std::runtime_error);
}

TEST(ModelBaseU1Spin1D, isValidQNumber) {
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.0));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.0));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.0));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -4.5));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.5));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.5));

   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +10.0));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +0.0 ));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, -10.0));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -9.0 ));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -1.0 ));
   EXPECT_TRUE(BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, +9.0 ));
   
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.5));
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.5));
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.5));
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -5.0));
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.0));
   EXPECT_FALSE(BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.0));
}

TEST(ModelBaseU1Spin1D, CalculateTargetDim) {
   //Spin-1/2
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(0, 0.5, +0.0), 0);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(1, 0.5, +0.5), 1);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(2, 0.5, +0.0), 2);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(3, 0.5, +0.5), 3);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, +0.0), 6);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, +2.0), 1);
   EXPECT_EQ(BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, -2.0), 1);
}

#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_1D_HPP_ */
