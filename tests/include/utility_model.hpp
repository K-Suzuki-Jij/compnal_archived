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
//  Created by Kohei Suzuki on 2022/01/29.
//

#ifndef COMPNAL_TEST_UTILITY_UTILITY_MODEL_HPP_
#define COMPNAL_TEST_UTILITY_UTILITY_MODEL_HPP_

#include "../../src/model/all.hpp"
#include "../../src/sparse_matrix/all.hpp"
#include "./test.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

template<typename RealType>
void TestElectron(const model::BaseU1Electron_1D<RealType> &model) {
   
   const double threshold = std::pow(10, -15);
   EXPECT_EQ(model.GetDimOnsite(), 4);
   
   const sparse_matrix::CRS<RealType> ref_c_up ({
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0},
      {0.0, 0.0, 0.0, 0.0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_up_d ({
      {0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_down ({
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, -1.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_down_d ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {0.0, -1.0, 0.0, 0.0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_nc_up ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_nc_down ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 1.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_nc ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 2.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sz ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.5, 0.0, 0.0},
      {0.0, 0.0, -0.5, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sx ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.5, 0.0},
      {0.0, 0.5, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_isy ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.5, 0.0},
      {0.0, -0.5, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sp ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sm ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   ExpectNear(model.GetOnsiteOperatorCUp()        , ref_c_up    , threshold);
   ExpectNear(model.GetOnsiteOperatorCUpDagger()  , ref_c_up_d  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDown()      , ref_c_down  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDownDagger(), ref_c_down_d, threshold);
   
   ExpectNear(model.GetOnsiteOperatorNCUp()  , ref_nc_up  , threshold);
   ExpectNear(model.GetOnsiteOperatorNCDown(), ref_nc_down, threshold);
   ExpectNear(model.GetOnsiteOperatorNC()    , ref_nc     , threshold);
   ExpectNear(model.GetOnsiteOperatorSz()    , ref_sz     , threshold);
   ExpectNear(model.GetOnsiteOperatorSx()    , ref_sx     , threshold);
   ExpectNear(model.GetOnsiteOperatoriSy()   , ref_isy    , threshold);
   ExpectNear(model.GetOnsiteOperatorSp()    , ref_sp     , threshold);
   ExpectNear(model.GetOnsiteOperatorSm()    , ref_sm     , threshold);
   
}

template<typename RealType>
void TestSpinOneHalf(const model::BaseU1Spin_1D<RealType> &model) {
   
   const double threshold = std::pow(10, -15);
   
   EXPECT_EQ(model.GetDimOnsite(), 2);
   EXPECT_EQ(model.GetMagnitudeSpin(), 0.5_hi);
   
   const sparse_matrix::CRS<RealType> ref_sp ({{+0.0, +1.0}, {+0.0, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sm ({{+0.0, +0.0}, {+1.0, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sx ({{+0.0, +0.5}, {+0.5, +0.0}});
   const sparse_matrix::CRS<RealType> ref_isy({{+0.0, +0.5}, {-0.5, +0.0}});
   const sparse_matrix::CRS<RealType> ref_sz ({{+0.5, +0.0}, {+0.0, -0.5}});
   
   ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , threshold);
   ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , threshold);
   ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , threshold);
   ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, threshold);
   ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , threshold);
}

template<typename RealType>
void TestSpinOne(const model::BaseU1Spin_1D<RealType> &model) {
   
   const double threshold = std::pow(10, -15);
   
   EXPECT_EQ(model.GetDimOnsite(), 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1_hi);
   
   const sparse_matrix::CRS<RealType> ref_sp ({
      {+0.0, +std::sqrt(2), +0.0},
      {+0.0, +0.0, +std::sqrt(2)},
      {+0.0, +0.0, +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sm ({
      {+0.0, +0.0, +0.0},
      {+std::sqrt(2), +0.0, +0.0},
      {+0.0, +std::sqrt(2), +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sx ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {+1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, +1.0/std::sqrt(2), +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_isy ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {-1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, -1.0/std::sqrt(2), +0.0}
   });
   
   const sparse_matrix::CRS<RealType> ref_sz ({
      {+1.0, +0.0, +0.0},
      {+0.0, +0.0, +0.0},
      {+0.0, +0.0, -1.0}
   });
   
   ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , threshold);
   ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , threshold);
   ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , threshold);
   ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, threshold);
   ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , threshold);
}

template<typename RealType>
void TestSpinOneHalf(const model::BaseU1SpinElectron_1D<RealType> &model) {
   
   const double threshold = std::pow(10, -15);
   
   EXPECT_EQ(model.GetDimOnsite(), 4*2);
   EXPECT_EQ(model.GetMagnitudeLSpin(), 0.5);
   EXPECT_EQ(model.GetDimOnsiteELectron(), 4);
   EXPECT_EQ(model.GetDimOnsiteLSpin(), 2);

   const sparse_matrix::CRS<RealType> ref_c_up ({
      {0, 0, 1, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 0},
      {0, 0, 0, 0, 0, 0, 0, 1},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_up_d ({
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_down ({
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0},
      {0, 0, 0, 0, 0, 0, -1, 0},
      {0, 0, 0, 0, 0, 0, 0, -1},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}
   }, sparse_matrix::CRSTag::FERMION);
   
   const sparse_matrix::CRS<RealType> ref_c_down_d ({
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},
      {1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, -1, 0, 0, 0, 0, 0},
      {0, 0, 0, -1, 0, 0, 0, 0}
   }, sparse_matrix::CRSTag::FERMION);
 
   ExpectNear(model.GetOnsiteOperatorCUp()        , ref_c_up    , threshold);
   ExpectNear(model.GetOnsiteOperatorCUpDagger()  , ref_c_up_d  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDown()      , ref_c_down  , threshold);
   ExpectNear(model.GetOnsiteOperatorCDownDagger(), ref_c_down_d, threshold);

}

}
}


#endif /* COMPNAL_TEST_UTILITY_UTILITY_MODEL_HPP_ */
