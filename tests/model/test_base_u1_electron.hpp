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

#ifndef COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_HPP_

#include <gtest/gtest.h>

#include "../../src/model/base_u1_electron.hpp"

namespace compnal {
namespace test {

TEST(ModelBaseU1Electron, Constructors) {
   using RealType = double;
   model::BaseU1Electron<RealType> model;

   EXPECT_EQ(model.GetDimOnsite(), 4);
   EXPECT_EQ(model.GetTotalSz(), 0);
   EXPECT_EQ(model.GetTotalElectron(), 0);

   EXPECT_EQ(model.GetOnsiteOperatorCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUp());
   EXPECT_EQ(model.GetOnsiteOperatorCUpDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUpDagger());
   EXPECT_EQ(model.GetOnsiteOperatorCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDown());
   EXPECT_EQ(model.GetOnsiteOperatorCDownDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDownDagger());
   EXPECT_EQ(model.GetOnsiteOperatorNCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCUp());
   EXPECT_EQ(model.GetOnsiteOperatorNCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCDown());
   EXPECT_EQ(model.GetOnsiteOperatorNC(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNC());

   EXPECT_EQ(model.GetOnsiteOperatorSp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSp());
   EXPECT_EQ(model.GetOnsiteOperatorSm(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSm());
   EXPECT_EQ(model.GetOnsiteOperatorSx(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSx());
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Electron<RealType>::CreateOnsiteOperatoriSy());
   EXPECT_EQ(model.GetOnsiteOperatorSz(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSz());
}

TEST(ModelBaseU1Electron, ConstructorsElectron) {
   using RealType = long double;
   model::BaseU1Electron<RealType> model(3);
   EXPECT_EQ(model.GetDimOnsite(), 4);
   EXPECT_EQ(model.GetTotalSz(), 0);
   EXPECT_EQ(model.GetTotalElectron(), 3);

   EXPECT_EQ(model.GetOnsiteOperatorCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUp());
   EXPECT_EQ(model.GetOnsiteOperatorCUpDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUpDagger());
   EXPECT_EQ(model.GetOnsiteOperatorCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDown());
   EXPECT_EQ(model.GetOnsiteOperatorCDownDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDownDagger());
   EXPECT_EQ(model.GetOnsiteOperatorNCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCUp());
   EXPECT_EQ(model.GetOnsiteOperatorNCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCDown());
   EXPECT_EQ(model.GetOnsiteOperatorNC(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNC());

   EXPECT_EQ(model.GetOnsiteOperatorSp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSp());
   EXPECT_EQ(model.GetOnsiteOperatorSm(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSm());
   EXPECT_EQ(model.GetOnsiteOperatorSx(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSx());
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Electron<RealType>::CreateOnsiteOperatoriSy());
   EXPECT_EQ(model.GetOnsiteOperatorSz(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSz());
}

TEST(ModelBaseU1Electron, ConstructorsElectronSz) {
   using RealType = long double;
   model::BaseU1Electron<RealType> model(3, 1.5);
   EXPECT_EQ(model.GetDimOnsite(), 4);
   EXPECT_EQ(model.GetTotalSz(), 1.5);
   EXPECT_EQ(model.GetTotalElectron(), 3);

   EXPECT_EQ(model.GetOnsiteOperatorCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUp());
   EXPECT_EQ(model.GetOnsiteOperatorCUpDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCUpDagger());
   EXPECT_EQ(model.GetOnsiteOperatorCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDown());
   EXPECT_EQ(model.GetOnsiteOperatorCDownDagger(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorCDownDagger());
   EXPECT_EQ(model.GetOnsiteOperatorNCUp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCUp());
   EXPECT_EQ(model.GetOnsiteOperatorNCDown(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNCDown());
   EXPECT_EQ(model.GetOnsiteOperatorNC(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorNC());

   EXPECT_EQ(model.GetOnsiteOperatorSp(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSp());
   EXPECT_EQ(model.GetOnsiteOperatorSm(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSm());
   EXPECT_EQ(model.GetOnsiteOperatorSx(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSx());
   EXPECT_EQ(model.GetOnsiteOperatoriSy(), model::BaseU1Electron<RealType>::CreateOnsiteOperatoriSy());
   EXPECT_EQ(model.GetOnsiteOperatorSz(), model::BaseU1Electron<RealType>::CreateOnsiteOperatorSz());
}

TEST(ModelBaseU1Electron, SetTotalSz) {
   model::BaseU1Electron<double> model(3, 2);
   EXPECT_EQ(model.GetTotalSz(), 2);
   model.SetTotalSz(4);
   EXPECT_EQ(model.GetTotalSz(), 4);
}

TEST(ModelBaseU1Electron, SetTotalElectron) {
   model::BaseU1Electron<double> model(3, 2);
   EXPECT_EQ(model.GetTotalElectron(), 3);
   model.SetTotalElectron(2);
   EXPECT_EQ(model.GetTotalElectron(), 2);
   EXPECT_THROW(model.SetTotalElectron(-2), std::runtime_error);
}

TEST(ModelBaseU1Electron, CalculateNumElectron) {
   model::BaseU1Electron<double> model;
   EXPECT_EQ(model.CalculateNumElectron(0), 0);
   EXPECT_EQ(model.CalculateNumElectron(1), 1);
   EXPECT_EQ(model.CalculateNumElectron(2), 1);
   EXPECT_EQ(model.CalculateNumElectron(3), 2);
   EXPECT_THROW(model.CalculateNumElectron(4), std::runtime_error);
   EXPECT_THROW(model.CalculateNumElectron(-1), std::runtime_error);
}

TEST(ModelBaseU1Electron, PrintBasisOnsite) {
   model::BaseU1Electron<double> model;
   EXPECT_NO_THROW(model.PrintBasisOnsite());
}

TEST(ModelBaseU1Electron, CalculateQNumber) {
   using BU1E = model::BaseU1Electron<double>;
   BU1E model(2, 1.5);
   EXPECT_EQ(model.CalculateQNumber(0, 0), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() + 0));
   EXPECT_EQ(model.CalculateQNumber(0, 1), BU1E::QType(model.GetTotalElectron() - 1, model.GetTotalSz() - 0.5));
   EXPECT_EQ(model.CalculateQNumber(0, 2), BU1E::QType(model.GetTotalElectron() - 1, model.GetTotalSz() + 0.5));
   EXPECT_EQ(model.CalculateQNumber(0, 3), BU1E::QType(model.GetTotalElectron() - 2, model.GetTotalSz() + 0));
   EXPECT_EQ(model.CalculateQNumber(1, 0), BU1E::QType(model.GetTotalElectron() + 1, model.GetTotalSz() + 0.5));
   EXPECT_EQ(model.CalculateQNumber(1, 1), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() + 0));
   EXPECT_EQ(model.CalculateQNumber(1, 2), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() + 1));
   EXPECT_EQ(model.CalculateQNumber(1, 3), BU1E::QType(model.GetTotalElectron() - 1, model.GetTotalSz() + 0.5));
   EXPECT_EQ(model.CalculateQNumber(2, 0), BU1E::QType(model.GetTotalElectron() + 1, model.GetTotalSz() - 0.5));
   EXPECT_EQ(model.CalculateQNumber(2, 1), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() - 1));
   EXPECT_EQ(model.CalculateQNumber(2, 2), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() + 0));
   EXPECT_EQ(model.CalculateQNumber(2, 3), BU1E::QType(model.GetTotalElectron() - 1, model.GetTotalSz() - 0.5));
   EXPECT_EQ(model.CalculateQNumber(3, 0), BU1E::QType(model.GetTotalElectron() + 2, model.GetTotalSz() + 0));
   EXPECT_EQ(model.CalculateQNumber(3, 1), BU1E::QType(model.GetTotalElectron() + 1, model.GetTotalSz() - 0.5));
   EXPECT_EQ(model.CalculateQNumber(3, 2), BU1E::QType(model.GetTotalElectron() + 1, model.GetTotalSz() + 0.5));
   EXPECT_EQ(model.CalculateQNumber(3, 3), BU1E::QType(model.GetTotalElectron() + 0, model.GetTotalSz() + 0));
}

TEST(ModelBaseU1Electron, GenerateBasis) {
   model::BaseU1Electron<double> model(4, 2);
   const auto basis_2_0 = model.GenerateBasis(2, {2, 0});
   EXPECT_EQ(basis_2_0.size(), std::size_t{4});
   EXPECT_EQ(basis_2_0.at(0), 3);
   EXPECT_EQ(basis_2_0.at(1), 6);
   EXPECT_EQ(basis_2_0.at(2), 9);
   EXPECT_EQ(basis_2_0.at(3), 12);

   const auto basis_2_p1 = model.GenerateBasis(2, {2, +1});
   EXPECT_EQ(basis_2_p1.size(), std::size_t{1});
   EXPECT_EQ(basis_2_p1.at(0), 5);

   const auto basis_2_m1 = model.GenerateBasis(2, {2, -1});
   EXPECT_EQ(basis_2_m1.size(), std::size_t{1});
   EXPECT_EQ(basis_2_m1.at(0), 10);

   const auto basis_1_p1 = model.GenerateBasis(2, {1, +0.5});
   EXPECT_EQ(basis_1_p1.size(), std::size_t{2});
   EXPECT_EQ(basis_1_p1.at(0), 1);
   EXPECT_EQ(basis_1_p1.at(1), 4);

   const auto basis_1_m1 = model.GenerateBasis(2, {1, -0.5});
   EXPECT_EQ(basis_1_m1.size(), std::size_t{2});
   EXPECT_EQ(basis_1_m1.at(0), 2);
   EXPECT_EQ(basis_1_m1.at(1), 8);

   EXPECT_THROW(model.GenerateBasis(2, {2, 0.5}), std::runtime_error);
}

TEST(ModelBaseU1Electron, ValidateQNumber) {
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 4, 0));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 3, +0.5));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 3, -0.5));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 2, 0));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 2, +1));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 2, -1));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 1, +0.5));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 1, -0.5));
   EXPECT_TRUE(model::BaseU1Electron<double>::ValidateQNumber(2, 0, 0));

   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +4, +1));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +3, +0));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +3, -1.5));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +2, +0.5));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +2, -0.5));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +1, +1));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +1, -1));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, +0, +1));
   EXPECT_FALSE(model::BaseU1Electron<double>::ValidateQNumber(2, -2, +0));
}

TEST(ModelBaseU1Electron, CalculateTargetDim) {
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(1, 0, 0), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(1, 1, -0.5), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(1, 1, +0.5), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(1, 2, 0), 1);

   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 0, 0), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 1, -0.5), 2);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 1, +0.5), 2);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 2, -1), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 2, +0), 4);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 2, +1), 1);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 3, -0.5), 2);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 3, +0.5), 2);
   EXPECT_EQ(model::BaseU1Electron<double>::CalculateTargetDim(2, 4, 0), 1);
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_MODEL_BASE_U1_ELECTRON_HPP_ */
