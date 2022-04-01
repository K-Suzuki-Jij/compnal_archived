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

TEST(ModelBaseU1Spin, Constructors) {
   using RealType = double;
   
   model::BaseU1Spin<RealType> model(0.5);
   EXPECT_EQ(model.GetDimOnsite()    , 2  );
   EXPECT_EQ(model.GetMagnitudeSpin(), 0.5);
   
   EXPECT_EQ(model.GetOnsiteOperatorSp (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSp (0.5));
   EXPECT_EQ(model.GetOnsiteOperatorSm (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSm (0.5));
   EXPECT_EQ(model.GetOnsiteOperatorSx (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSx (0.5));
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Spin<RealType>::CreateOnsiteOperatoriSy(0.5));
   EXPECT_EQ(model.GetOnsiteOperatorSz (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSz (0.5));
}

TEST(ModelBaseU1Spin, ConstructorsSpin) {
   using RealType = long double;
   
   model::BaseU1Spin<RealType> model(1);
   EXPECT_EQ(model.GetDimOnsite()    , 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1);
   
   EXPECT_EQ(model.GetOnsiteOperatorSp (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSp (1));
   EXPECT_EQ(model.GetOnsiteOperatorSm (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSm (1));
   EXPECT_EQ(model.GetOnsiteOperatorSx (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSx (1));
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Spin<RealType>::CreateOnsiteOperatoriSy(1));
   EXPECT_EQ(model.GetOnsiteOperatorSz (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSz (1));
}

TEST(ModelBaseU1Spin, ConstructorsSpinSz) {
   model::BaseU1Spin<double> model(1);
   EXPECT_EQ(model.GetDimOnsite()    , 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1);
}

TEST(ModelBaseU1Spin, SetMagnitudeSpin) {
   using RealType = long double;
   
   model::BaseU1Spin<RealType> model;
   model.SetMagnitudeSpin(1.5);
   EXPECT_EQ(model.GetDimOnsite()    , 4);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.5);
   EXPECT_EQ(model.GetOnsiteOperatorSp (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSp (1.5));
   EXPECT_EQ(model.GetOnsiteOperatorSm (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSm (1.5));
   EXPECT_EQ(model.GetOnsiteOperatorSx (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSx (1.5));
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Spin<RealType>::CreateOnsiteOperatoriSy(1.5));
   EXPECT_EQ(model.GetOnsiteOperatorSz (), model::BaseU1Spin<RealType>::CreateOnsiteOperatorSz (1.5));
   
   EXPECT_THROW(model.SetMagnitudeSpin(0) , std::runtime_error);
   EXPECT_THROW(model.SetMagnitudeSpin(-1), std::runtime_error);
}

TEST(ModelBaseU1Spin, CalculateNumElectron) {
   model::BaseU1Spin<double> model;
   EXPECT_THROW(model.CalculateNumElectron(-1), std::runtime_error);
   EXPECT_EQ(model.CalculateNumElectron(0), 0);
   EXPECT_EQ(model.CalculateNumElectron(1), 0);
   EXPECT_THROW(model.CalculateNumElectron(2) , std::runtime_error);
}

TEST(ModelBaseU1Spin, CalculateQNumber) {
   model::BaseU1Spin<double> model;
   model::BaseU1Spin<double>::QType total_sz = 3;
   EXPECT_THROW(model.CalculateQNumber(-1,  0, total_sz), std::runtime_error);
   EXPECT_THROW(model.CalculateQNumber( 0, -1, total_sz), std::runtime_error);
   EXPECT_THROW(model.CalculateQNumber(-2, -1, total_sz), std::runtime_error);
   
   EXPECT_THROW(model.CalculateQNumber(2, 0, total_sz), std::runtime_error);
   EXPECT_THROW(model.CalculateQNumber(0, 2, total_sz), std::runtime_error);
   EXPECT_THROW(model.CalculateQNumber(2, 2, total_sz), std::runtime_error);

   EXPECT_EQ(model.CalculateQNumber(0, 0, total_sz), 3);
   EXPECT_EQ(model.CalculateQNumber(0, 1, total_sz), 4);
   EXPECT_EQ(model.CalculateQNumber(1, 0, total_sz), 2);
   EXPECT_EQ(model.CalculateQNumber(1, 1, total_sz), 3);
}

TEST(ModelBaseU1Spin, GenerateBasisSpin05) {
   model::BaseU1Spin<double> model(0.5);
   std::vector<std::int64_t> basis_spin05_sz_p1 = {0};
   std::vector<std::int64_t> basis_spin05_sz_0  = {1, 2};
   std::vector<std::int64_t> basis_spin05_sz_m1 = {3};
   EXPECT_THROW(model.GenerateBasis(2, -2.0, false), std::runtime_error);
   EXPECT_THROW(model.GenerateBasis(2, -1.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2, -1.0, false), basis_spin05_sz_m1);
   EXPECT_THROW(model.GenerateBasis(2, -0.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2,  0.0, false), basis_spin05_sz_0 );
   EXPECT_THROW(model.GenerateBasis(2, +0.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2, +1.0, false), basis_spin05_sz_p1);
   EXPECT_THROW(model.GenerateBasis(2, +1.5, false), std::runtime_error);
   EXPECT_THROW(model.GenerateBasis(2, +2.0, false), std::runtime_error);
}

TEST(ModelBaseU1Spin, GenerateBasisSpin1) {
   model::BaseU1Spin<double> model(1);
   std::vector<std::int64_t> basis_spin1_sz_p2 = {0};
   std::vector<std::int64_t> basis_spin1_sz_p1 = {1, 3};
   std::vector<std::int64_t> basis_spin1_sz_0  = {2, 4, 6};
   std::vector<std::int64_t> basis_spin1_sz_m1 = {5, 7};
   std::vector<std::int64_t> basis_spin1_sz_m2 = {8};
   EXPECT_EQ   (model.GenerateBasis(2, -2.0, false), basis_spin1_sz_m2);
   EXPECT_THROW(model.GenerateBasis(2, -1.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2, -1.0, false), basis_spin1_sz_m1);
   EXPECT_THROW(model.GenerateBasis(2, -0.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2,  0.0, false), basis_spin1_sz_0);
   EXPECT_THROW(model.GenerateBasis(2, +0.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2, +1.0, false), basis_spin1_sz_p1);
   EXPECT_THROW(model.GenerateBasis(2, +1.5, false), std::runtime_error);
   EXPECT_EQ   (model.GenerateBasis(2, +2.0, false), basis_spin1_sz_p2);
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

TEST(ModelBaseU1Spin, SpinOneHalf) {
   using RealType = double;
   
   const blas::CRS<RealType> ref_sp ({{+0.0, +1.0}, {+0.0, +0.0}}, blas::CRSTag::BOSON);
   const blas::CRS<RealType> ref_sm ({{+0.0, +0.0}, {+1.0, +0.0}}, blas::CRSTag::BOSON);
   const blas::CRS<RealType> ref_sx ({{+0.0, +0.5}, {+0.5, +0.0}}, blas::CRSTag::BOSON);
   const blas::CRS<RealType> ref_isy({{+0.0, +0.5}, {-0.5, +0.0}}, blas::CRSTag::BOSON);
   const blas::CRS<RealType> ref_sz ({{+0.5, +0.0}, {+0.0, -0.5}}, blas::CRSTag::BOSON);
   
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSp (0.5) , ref_sp );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSm (0.5) , ref_sm );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSx (0.5) , ref_sx );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatoriSy(0.5) , ref_isy);
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSz (0.5) , ref_sz );
}

TEST(ModelBaseU1Spin, SpinOne) {
   using RealType = double;
   const RealType sqrt2   = std::sqrt(RealType{2.0});
   const RealType sqrt2_2 = std::sqrt(RealType{2.0})/RealType{2.0};

   const blas::CRS<RealType> ref_sp ({
      {+0.0, +sqrt2, +0.0},
      {+0.0, +0.0, +sqrt2},
      {+0.0, +0.0, +0.0}
   }, blas::CRSTag::BOSON);
   
   const blas::CRS<RealType> ref_sm ({
      {+0.0, +0.0, +0.0},
      {+sqrt2, +0.0, +0.0},
      {+0.0, +sqrt2, +0.0}
   }, blas::CRSTag::BOSON);
   
   const blas::CRS<RealType> ref_sx ({
      {+0.0, +sqrt2_2, +0.0},
      {+sqrt2_2, +0.0, +sqrt2_2},
      {+0.0, +sqrt2_2, +0.0}
   }, blas::CRSTag::BOSON);
   
   const blas::CRS<RealType> ref_isy ({
      {+0.0, +sqrt2_2, +0.0},
      {-sqrt2_2, +0.0, +sqrt2_2},
      {+0.0, -sqrt2_2, +0.0}
   }, blas::CRSTag::BOSON);
   
   const blas::CRS<RealType> ref_sz ({
      {+1.0, +0.0, +0.0},
      {+0.0, +0.0, +0.0},
      {+0.0, +0.0, -1.0}
   }, blas::CRSTag::BOSON);
   
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSp (1) , ref_sp );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSm (1) , ref_sm );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSx (1) , ref_sx );
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatoriSy(1) , ref_isy);
   EXPECT_EQ(model::BaseU1Spin<RealType>::CreateOnsiteOperatorSz (1) , ref_sz );
}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_HPP_ */
