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

#ifndef COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_1D_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_1D_HPP_

#include "../../src/model/base_u1_electron_1d.hpp"
#include "../test.hpp"
#include <gtest/gtest.h>

template<typename RealType>
void TestElectron(const compnal::model::BaseU1Electron_1D<RealType> &model) {
   
   EXPECT_EQ(model.GetDimOnsite(), 4);
   
   const compnal::sparse_matrix::CRS<RealType> ref_c_up ({
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0},
      {0.0, 0.0, 0.0, 0.0}
   }, compnal::sparse_matrix::CRSTag::FERMION);
   
   const compnal::sparse_matrix::CRS<RealType> ref_c_up_d ({
      {0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0}
   }, compnal::sparse_matrix::CRSTag::FERMION);
   
   const compnal::sparse_matrix::CRS<RealType> ref_c_down ({
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, -1.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   }, compnal::sparse_matrix::CRSTag::FERMION);
   
   const compnal::sparse_matrix::CRS<RealType> ref_c_down_d ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {0.0, -1.0, 0.0, 0.0}
   }, compnal::sparse_matrix::CRSTag::FERMION);
   
   const compnal::sparse_matrix::CRS<RealType> ref_nc_up ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_nc_down ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 1.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_nc ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 2.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sz ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.5, 0.0, 0.0},
      {0.0, 0.0, -0.5, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sx ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.5, 0.0},
      {0.0, 0.5, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_isy ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.5, 0.0},
      {0.0, -0.5, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sp ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sm ({
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0}
   });
      
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorCUp()        , ref_c_up    , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorCUpDagger()  , ref_c_up_d  , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorCDown()      , ref_c_down  , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorCDownDagger(), ref_c_down_d, std::pow(10, -15)));
   
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorNCUp()  , ref_nc_up  , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorNCDown(), ref_nc_down, std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorNC()    , ref_nc     , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSz()    , ref_sz     , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSx()    , ref_sx     , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatoriSy()   , ref_isy    , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSp()    , ref_sp     , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSm()    , ref_sm     , std::pow(10, -15)));

}

TEST(ModelBaseU1Electron1D, ConstructorDefault) {
   compnal::model::BaseU1Electron_1D<double> model;
   TestElectron(model);
   
   EXPECT_EQ(model.GetSystemSize(), 0);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   EXPECT_EQ(model.GetTotalElectron(), 0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(),
                                       std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, compnal::utility::PairHash>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(),
                                       std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, compnal::utility::PairHash>{}));
}

TEST(ModelBaseU1Electron1D, ConstructorSystemSize) {
   compnal::model::BaseU1Electron_1D<double> model(8);
   TestElectron(model);
   
   EXPECT_EQ(model.GetSystemSize(), 8);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   EXPECT_EQ(model.GetTotalElectron(), 0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(),
                                       std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, compnal::utility::PairHash>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(),
                                       std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, compnal::utility::PairHash>{}));
}

TEST(ModelBaseU1Electron1D, ConstructorSystemSizeElectron) {
   compnal::model::BaseU1Electron_1D<double> model(8, 2);
   TestElectron(model);
   
   EXPECT_EQ(model.GetSystemSize(), 8);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   EXPECT_EQ(model.GetTotalElectron(), 2);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(),
                                       std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, compnal::utility::PairHash>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(),
                                       std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, compnal::utility::PairHash>{}));
}

TEST(ModelBaseU1Electron1D, SetSystemSize) {
   compnal::model::BaseU1Electron_1D<double> model;
   model.SetSystemSize(3);
   EXPECT_EQ(model.GetSystemSize(), 3);
}

TEST(ModelBaseU1Electron1D, SetTotalSz) {
   compnal::model::BaseU1Electron_1D<double> model;
   model.SetTotalSz(2.5);
   EXPECT_EQ(model.GetTotalSz(), 2.5);
}

TEST(ModelBaseU1Electron1D, SetTotalElectron) {
   compnal::model::BaseU1Electron_1D<double> model;
   model.SetTotalElectron(100);
   EXPECT_EQ(model.GetTotalElectron(), 100);
}

TEST(ModelBaseU1Electron1D, CalculateNumElectron) {
   compnal::model::BaseU1Electron_1D<double> model;
   EXPECT_EQ(model.CalculateNumElectron(0), 0);
   EXPECT_EQ(model.CalculateNumElectron(1), 1);
   EXPECT_EQ(model.CalculateNumElectron(2), 1);
   EXPECT_EQ(model.CalculateNumElectron(3), 2);
   EXPECT_THROW(model.CalculateNumElectron(4), std::runtime_error);
}

#endif /* COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_1D_HPP_ */
