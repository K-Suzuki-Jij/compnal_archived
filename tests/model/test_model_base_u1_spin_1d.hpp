//
//  test_model_base_u1_spin_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2022/01/08.
//

#ifndef TEST_MODEL_BASE_U1_SPIN_1d_HPP_
#define TEST_MODEL_BASE_U1_SPIN_1d_HPP_

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

TEST(Constructor, Default) {
   
   compnal::model::BaseU1Spin_1D<double> model;
   
   TestSpinOneHalf(model);
   
   EXPECT_EQ(model.GetSystemSize(), 0);
   EXPECT_EQ(model.GetTotalSz(), 0.0);
   
   EXPECT_EQ(model.GetCalculatedEigenvectorSet(), std::unordered_set<int>());
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBases(), std::unordered_map<int, std::vector<std::int64_t>>{}));
   EXPECT_TRUE(compnal::test::ExpectEQ(model.GetBasesInv(), std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>>{}));

}


#endif /* TEST_MODEL_BASE_U1_SPIN_1d_HPP_ */
