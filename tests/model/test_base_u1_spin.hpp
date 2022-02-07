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

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_HPP_

#include "../../src/model/base_u1_spin.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(ModelBaseU1Spin, SpinOneHalf) {
   using RealType = double;
   model::BaseU1Spin<RealType> model;

   EXPECT_EQ(model.GetDimOnsite(), 2);
   EXPECT_EQ(model.GetMagnitudeSpin(), 0.5);
   EXPECT_EQ(model.GetTotalSz(), 0);
   
   const sparse_matrix::CRS<RealType> ref_sp ({{+0.0, +1.0}, {+0.0, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sm ({{+0.0, +0.0}, {+1.0, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sx ({{+0.0, +0.5}, {+0.5, +0.0}});
   const sparse_matrix::CRS<RealType> ref_isy({{+0.0, +0.5}, {-0.5, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sz ({{+0.5, +0.0}, {+0.0, -0.5}});
   
   EXPECT_EQ(model.GetOnsiteOperatorSp() , ref_sp );
   EXPECT_EQ(model.GetOnsiteOperatorSm() , ref_sm );
   EXPECT_EQ(model.GetOnsiteOperatorSx() , ref_sx );
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), ref_isy);
   EXPECT_EQ(model.GetOnsiteOperatorSz() , ref_sz );
}

TEST(ModelBaseU1Spin, SpinOne) {
   using RealType = long double;
   model::BaseU1Spin<RealType> model(1);
   
   EXPECT_EQ(model.GetDimOnsite(), 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1);
   EXPECT_EQ(model.GetTotalSz(), 0);

   const RealType sqrt2 = std::sqrt(static_cast<RealType>(2));
   const sparse_matrix::CRS<RealType> ref_sp ({
      {+0.0, +sqrt2, +0.0},
      {+0.0, +0.0, +sqrt2},
      {+0.0, +0.0, +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sm ({
      {+0.0, +0.0, +0.0},
      {+sqrt2, +0.0, +0.0},
      {+0.0, +sqrt2, +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sx ({
      {+0.0, +sqrt2/2.0, +0.0},
      {+sqrt2/2.0, +0.0, +sqrt2/2.0},
      {+0.0, +sqrt2/2.0, +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_isy ({
      {+0.0, +sqrt2/2.0, +0.0},
      {-sqrt2/2.0, +0.0, +sqrt2/2.0},
      {+0.0, -sqrt2/2.0, +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sz ({
      {+1.0, +0.0, +0.0},
      {+0.0, +0.0, +0.0},
      {+0.0, +0.0, -1.0}
   });
   
   EXPECT_EQ(model.GetOnsiteOperatorSp() , ref_sp );
   EXPECT_EQ(model.GetOnsiteOperatorSm() , ref_sm );
   EXPECT_EQ(model.GetOnsiteOperatorSx() , ref_sx );
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), ref_isy);
   EXPECT_EQ(model.GetOnsiteOperatorSz() , ref_sz );
}

TEST(ModelBaseU1Spin, ValidateQNumber) {
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, +5.0));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, +0.0));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, -5.0));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, -4.5));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, -0.5));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, +4.5));

   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 1.0, +10.0));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 1.0, +0.0 ));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(10, 1.0, -10.0));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 1.0, -9.0 ));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 1.0, -1.0 ));
   EXPECT_TRUE(model::BaseU1Spin<double>::ValidateQNumber(9 , 1.0, +9.0 ));
   
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, +5.5));
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, +0.5));
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(10, 0.5, -5.5));
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, -5.0));
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, -0.0));
   EXPECT_FALSE(model::BaseU1Spin<double>::ValidateQNumber(9 , 0.5, +4.0));
}

TEST(ModelBaseU1Spin, CalculateTargetDim) {
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(0, 0.5, +0.0), 0);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(1, 0.5, +0.5), 1);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(2, 0.5, +0.0), 2);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(3, 0.5, +0.5), 3);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(4, 0.5, +0.0), 6);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(4, 0.5, +2.0), 1);
   EXPECT_EQ(model::BaseU1Spin<double>::CalculateTargetDim(4, 0.5, -2.0), 1);
}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_HPP_ */
