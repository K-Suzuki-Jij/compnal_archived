//
//  test_model_base_u1_spin_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2022/01/08.
//

#ifndef TEST_MODEL_BASE_U1_SPIN_1D_HPP_
#define TEST_MODEL_BASE_U1_SPIN_1D_HPP_

#include "../../src/model/base_u1_spin_1d.hpp"
#include "../test.hpp"
#include <gtest/gtest.h>

template<typename RealType>
void TestSpinOneHalf(const compnal::model::BaseU1Spin_1D<RealType> &model) {
   
   EXPECT_EQ(model.GetDimOnsite(), 2);
   EXPECT_EQ(model.GetMagnitudeSpin(), 0.5);
   
   const compnal::sparse_matrix::CRS<RealType> ref_sp ({{+0.0, +1.0}, {+0.0, +0.0}});
   const compnal::sparse_matrix::CRS<RealType> ref_sm ({{+0.0, +0.0}, {+1.0, +0.0}});
   const compnal::sparse_matrix::CRS<RealType> ref_sx ({{+0.0, +0.5}, {+0.5, +0.0}});
   const compnal::sparse_matrix::CRS<RealType> ref_isy({{+0.0, +0.5}, {-0.5, +0.0}});
   const compnal::sparse_matrix::CRS<RealType> ref_sz ({{+0.5, +0.0}, {+0.0, -0.5}});

   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , std::pow(10, -15)));
}

template<typename RealType>
void TestSpinOne(const compnal::model::BaseU1Spin_1D<RealType> &model) {
   
   EXPECT_EQ(model.GetDimOnsite(), 3);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.0);
   
   const compnal::sparse_matrix::CRS<RealType> ref_sp ({
      {+0.0, +std::sqrt(2), +0.0},
      {+0.0, +0.0, +std::sqrt(2)},
      {+0.0, +0.0, +0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sm ({
      {+0.0, +0.0, +0.0},
      {+std::sqrt(2), +0.0, +0.0},
      {+0.0, +std::sqrt(2), +0.0}
   });

   const compnal::sparse_matrix::CRS<RealType> ref_sx ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {+1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, +1.0/std::sqrt(2), +0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_isy ({
      {+0.0, +1.0/std::sqrt(2), +0.0},
      {-1.0/std::sqrt(2), +0.0, +1.0/std::sqrt(2)},
      {+0.0, -1.0/std::sqrt(2), +0.0}
   });
   
   const compnal::sparse_matrix::CRS<RealType> ref_sz ({
      {+1.0, +0.0, +0.0},
      {+0.0, +0.0, +0.0},
      {+0.0, +0.0, -1.0}
   });

   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSp() , ref_sp , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSm() , ref_sm , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSx() , ref_sx , std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatoriSy(), ref_isy, std::pow(10, -15)));
   EXPECT_TRUE(compnal::test::ExpectNear(model.GetOnsiteOperatorSz() , ref_sz , std::pow(10, -15)));
}


TEST(BaseU1Spin1D, ConstructorDefault) {
   compnal::model::BaseU1Spin_1D<double> model;
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 0);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{}));
}

TEST(BaseU1Spin1D, ConstructorSystemSize) {
   compnal::model::BaseU1Spin_1D<double> model(10);
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{}));
}

TEST(BaseU1Spin1D, ConstructorSystemSizeSpin) {
   compnal::model::BaseU1Spin_1D<double> model(10, 1);
   TestSpinOne(model);
   
   EXPECT_EQ(model.GetSystemSize(), 10);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{}));
}

TEST(BaseU1Spin1D, SetSystemSize) {
   compnal::model::BaseU1Spin_1D<double> model;
   model.SetSystemSize(5);
   EXPECT_EQ(model.GetSystemSize(), 5);
}

TEST(BaseU1Spin1D, SetTotalSz) {
   compnal::model::BaseU1Spin_1D<double> model;
   model.SetTotalSz(2);
   EXPECT_EQ(model.GetTotalSz(), 2);
   EXPECT_THROW(model.SetTotalSz(1.9), std::runtime_error);
}

TEST(BaseU1Spin1D, SetMagnitudeSpin) {
   compnal::model::BaseU1Spin_1D<double> model;
   model.SetMagnitudeSpin(1.5);
   EXPECT_EQ(model.GetMagnitudeSpin(), 1.5);
   EXPECT_THROW(model.SetMagnitudeSpin(1.3), std::runtime_error);
}

TEST(BaseU1Spin1D, isValidQNumber) {
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -4.5));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.5));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.5));

   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +10.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, +0.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 1.0, -10.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -9.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, -1.0));
   EXPECT_TRUE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 1.0, +9.0));
   
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +5.5));
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, +0.5));
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(10, 0.5, -5.5));
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -5.0));
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, -0.0));
   EXPECT_FALSE(compnal::model::BaseU1Spin_1D<double>::isValidQNumber(9 , 0.5, +4.0));
}

TEST(BaseU1Spin1D, CalculateTargetDim) {
   //Spin-1/2
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(0, 0.5, 0.0), 0);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(1, 0.5, 0.5), 1);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(2, 0.5, 0.0), 2);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(3, 0.5, 0.5), 3);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, 0.0), 6);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, 2.0), 1);
   EXPECT_EQ(compnal::model::BaseU1Spin_1D<double>::CalculateTargetDim(4, 0.5, -2.0), 1);
}

#endif /* TEST_MODEL_BASE_U1_SPIN_1D_HPP_ */
