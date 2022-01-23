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
//  Created by Kohei Suzuki on 2022/01/23.
//

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_ELECTRON_1D_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_ELECTRON_1D_HPP_

#include "../../src/model/base_u1_spin_electron_1d.hpp"
#include "../test.hpp"
#include <gtest/gtest.h>

namespace {

using compnal::LInt;
using compnal::IntPair;
using compnal::Map;
using compnal::PairHash;
using compnal::test::ExpectEQ;
using compnal::test::ExpectNear;
using compnal::model::BaseU1SpinElectron_1D;

}

template<typename RealType>
void TestSpinOneHalf(const BaseU1SpinElectron_1D<RealType> &model) {
   
   const double threshold = std::pow(10, -15);
   
   EXPECT_EQ(model.GetDimOnsite(), 4*2);
   EXPECT_EQ(model.GetMagnitudeLSpin(), 0.5);
   EXPECT_EQ(model.GetDimOnsiteELectron(), 4);
   EXPECT_EQ(model.GetDimOnsiteLSpin(), 2);

   const CRS<RealType> ref_c_up ({
      {0, 0, 1, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 0},
      {0, 0, 0, 0, 0, 0, 0, 1},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}
   }, CRSTag::FERMION);
   
   const CRS<RealType> ref_c_up_d ({
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0}
   }, CRSTag::FERMION);
   
   const CRS<RealType> ref_c_down ({
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0},
      {0, 0, 0, 0, 0, 0, -1, 0},
      {0, 0, 0, 0, 0, 0, 0, -1},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}
   }, CRSTag::FERMION);
   
   const CRS<RealType> ref_c_down_d ({
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, -1, 0, 0, 0, 0, 0},
      {0, 0, 0, -1, 0, 0, 0, 0}
   }, CRSTag::FERMION);
 
   ExpectNear(model.GetOnsiteOperatorCUp()        , ref_c_up    , threshold);
   ExpectNear(model.GetOnsiteOperatorCUpDagger()  , ref_c_up_d  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDown()      , ref_c_down  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDownDagger(), ref_c_down_d, threshold);

}

TEST(ModelBaseU1SpinElectron1D, ConstructorDefault) {
   BaseU1SpinElectron_1D<double> model;
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize()    , 0  );
   EXPECT_EQ(model.GetTotalSz()       , 0.0);
   EXPECT_EQ(model.GetTotalElectron() , 0  );
   EXPECT_EQ(model.GetMagnitudeLSpin(), 0.5);
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   
   ExpectEQ(model.GetBases(), std::unordered_map<IntPair, std::vector<LInt>, PairHash>{});
   ExpectEQ(model.GetBasesInv(), Map<IntPair, Map<LInt, LInt>, PairHash>{});
   
}


#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_1D_HPP_ */
   
