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
//  Created by Kohei Suzuki on 2022/02/07.
//

#ifndef COMPNAL_TEST_MODEL_GENERAL_MODEL_HPP_
#define COMPNAL_TEST_MODEL_GENERAL_MODEL_HPP_

#include "../../src/model/base_u1_spin.hpp"
#include "../../src/model/base_u1_electron.hpp"
#include "../../src/model/base_u1_spin_electron.hpp"
#include "../../src/model/base_u1_spin_multi_electrons.hpp"
#include "../../src/model/general_model.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(ModelGeneralModel, Basic) {
   
   auto check_model = [](auto model) {
      using VariantVecType = std::vector<std::variant<std::int64_t, std::string>>;
      
      const auto sz = model.GetOnsiteOperatorSz();
      const auto sp = model.GetOnsiteOperatorSp();
      const auto sm = model.GetOnsiteOperatorSm();
      
      model.AddPotential(-1 , sz);
      model.AddPotential(-1 , sm);
      model.AddPotential("a", sp);
      model.AddPotential(VariantVecType{1, "a"}, sz);
      
      model.AddInteraction(1, sp, -2 , sm);
      model.AddInteraction(1, sp, "a", sm);
      model.AddInteraction(1, sm, "a", sp);
      model.AddInteraction(VariantVecType{1, "a"}, sp, 1, sm);
      model.AddInteraction(VariantVecType{1, "a"}, sp, VariantVecType{1, "a"}, sz);
      
      
      EXPECT_EQ(model.GetSystemSize(), 5);
      EXPECT_EQ(model.GetIndexList().count(+1), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetIndexList().count(-1), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetIndexList().count(-2), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetIndexList().count("a"), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetIndexList().count(VariantVecType{1, "a"}), static_cast<unsigned long>(1));
      
      EXPECT_EQ(model.GetPotentialList().size(), static_cast<unsigned long>(3));
      EXPECT_EQ(model.GetPotential(-1), sz + sm);
      EXPECT_EQ(model.GetPotential("a"), sp);
      EXPECT_EQ(model.GetPotential(VariantVecType{1, "a"}), sz + sp*sz);
      
      EXPECT_EQ(model.GetInteraction(1, -2).size(), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetInteraction(1, -2).at(0).first , sp);
      EXPECT_EQ(model.GetInteraction(1, -2).at(0).second, sm);
      
      EXPECT_EQ(model.GetInteraction(1, "a").size(), static_cast<unsigned long>(2));
      EXPECT_EQ(model.GetInteraction(1, "a").at(0).first , sp);
      EXPECT_EQ(model.GetInteraction(1, "a").at(0).second, sm);
      EXPECT_EQ(model.GetInteraction(1, "a").at(1).first , sm);
      EXPECT_EQ(model.GetInteraction(1, "a").at(1).second, sp);
      
      EXPECT_EQ(model.GetInteraction(VariantVecType{1, "a"}, 1).size(), static_cast<unsigned long>(1));
      EXPECT_EQ(model.GetInteraction(VariantVecType{1, "a"}, 1).at(0).first , sp);
      EXPECT_EQ(model.GetInteraction(VariantVecType{1, "a"}, 1).at(0).second, sm);
   };
   
   check_model(model::GeneralModel<model::BaseU1Spin<float>>());
   check_model(model::GeneralModel<model::BaseU1Electron<float>>());
   check_model(model::GeneralModel<model::BaseU1SpinElectron<float>>());
   check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<float>>());
   
   check_model(model::GeneralModel<model::BaseU1Spin<double>>());
   check_model(model::GeneralModel<model::BaseU1Electron<double>>());
   check_model(model::GeneralModel<model::BaseU1SpinElectron<double>>());
   check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<double>>());
   
   check_model(model::GeneralModel<model::BaseU1Spin<long double>>());
   check_model(model::GeneralModel<model::BaseU1Electron<long double>>());
   check_model(model::GeneralModel<model::BaseU1SpinElectron<long double>>());
   check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<long double>>());
   
}

TEST(ModelGeneralModel, BaseU1Spin) {
   
   auto check_model = [](auto model) {
      model.AddPotential(-1, model.GetOnsiteOperatorSx());
      model.AddPotential(-1, model.GetOnsiteOperatoriSy());
      model.AddPotential(-1, model.GetOnsiteOperatorSz());
      model.AddPotential(-1, model.GetOnsiteOperatorSp());
      model.AddPotential(-1, model.GetOnsiteOperatorSm());
   };
   
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Spin<float>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Spin<double>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Spin<long double>>()));
   
}

TEST(ModelGeneralModel, BaseU1Electron) {
   
   auto check_model = [](auto model) {
      model.AddPotential(-1, model.GetOnsiteOperatorCUp());
      model.AddPotential(-1, model.GetOnsiteOperatorCUpDagger());
      model.AddPotential(-1, model.GetOnsiteOperatorCDown());
      model.AddPotential(-1, model.GetOnsiteOperatorCDownDagger());
      model.AddPotential(-1, model.GetOnsiteOperatorNCUp());
      model.AddPotential(-1, model.GetOnsiteOperatorNCDown());
      model.AddPotential(-1, model.GetOnsiteOperatorNC());
      model.AddPotential(-1, model.GetOnsiteOperatorSx());
      model.AddPotential(-1, model.GetOnsiteOperatoriSy());
      model.AddPotential(-1, model.GetOnsiteOperatorSz());
      model.AddPotential(-1, model.GetOnsiteOperatorSp());
      model.AddPotential(-1, model.GetOnsiteOperatorSm());
   };
   
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Electron<float>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Electron<double>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1Electron<long double>>()));
   
}

TEST(ModelGeneralModel, BaseU1SpinElectron) {
   
   auto check_model = [](auto model) {
      model.AddPotential(-1, model.GetOnsiteOperatorCUp());
      model.AddPotential(-1, model.GetOnsiteOperatorCUpDagger());
      model.AddPotential(-1, model.GetOnsiteOperatorCDown());
      model.AddPotential(-1, model.GetOnsiteOperatorCDownDagger());
      model.AddPotential(-1, model.GetOnsiteOperatorNCUp());
      model.AddPotential(-1, model.GetOnsiteOperatorNCDown());
      model.AddPotential(-1, model.GetOnsiteOperatorNC());
      model.AddPotential(-1, model.GetOnsiteOperatorSxC());
      model.AddPotential(-1, model.GetOnsiteOperatoriSyC());
      model.AddPotential(-1, model.GetOnsiteOperatorSzC());
      model.AddPotential(-1, model.GetOnsiteOperatorSpC());
      model.AddPotential(-1, model.GetOnsiteOperatorSmC());
      
      model.AddPotential(-1, model.GetOnsiteOperatorSxL());
      model.AddPotential(-1, model.GetOnsiteOperatoriSyL());
      model.AddPotential(-1, model.GetOnsiteOperatorSzL());
      model.AddPotential(-1, model.GetOnsiteOperatorSpL());
      model.AddPotential(-1, model.GetOnsiteOperatorSmL());
      model.AddPotential(-1, model.GetOnsiteOperatorSCSL());
      
      model.AddPotential(-1, model.GetOnsiteOperatorSx());
      model.AddPotential(-1, model.GetOnsiteOperatoriSy());
      model.AddPotential(-1, model.GetOnsiteOperatorSz());
      model.AddPotential(-1, model.GetOnsiteOperatorSp());
      model.AddPotential(-1, model.GetOnsiteOperatorSm());
   };
   
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinElectron<float>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinElectron<double>>()));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinElectron<long double>>()));
   
}

TEST(ModelGeneralModel, BaseU1SpinMultiElectrons) {
   
   auto check_model = [](auto model) {
      model.AddPotential(-1, model.GetOnsiteOperatorCUp(0));
      model.AddPotential(-1, model.GetOnsiteOperatorCUpDagger(0));
      model.AddPotential(-1, model.GetOnsiteOperatorCDown(0));
      model.AddPotential(-1, model.GetOnsiteOperatorCDownDagger(0));
      model.AddPotential(-1, model.GetOnsiteOperatorNCUp(0));
      model.AddPotential(-1, model.GetOnsiteOperatorNCDown(0));
      model.AddPotential(-1, model.GetOnsiteOperatorNC(0));
      model.AddPotential(-1, model.GetOnsiteOperatorSxC(0));
      model.AddPotential(-1, model.GetOnsiteOperatoriSyC(0));
      model.AddPotential(-1, model.GetOnsiteOperatorSzC(0));
      model.AddPotential(-1, model.GetOnsiteOperatorSpC(0));
      model.AddPotential(-1, model.GetOnsiteOperatorSmC(0));
      model.AddPotential(-1, model.GetOnsiteOperatorSCSL(0));
      
      model.AddPotential(-1, model.GetOnsiteOperatorCUp(1));
      model.AddPotential(-1, model.GetOnsiteOperatorCUpDagger(1));
      model.AddPotential(-1, model.GetOnsiteOperatorCDown(1));
      model.AddPotential(-1, model.GetOnsiteOperatorCDownDagger(1));
      model.AddPotential(-1, model.GetOnsiteOperatorNCUp(1));
      model.AddPotential(-1, model.GetOnsiteOperatorNCDown(1));
      model.AddPotential(-1, model.GetOnsiteOperatorNC(1));
      model.AddPotential(-1, model.GetOnsiteOperatorSxC(1));
      model.AddPotential(-1, model.GetOnsiteOperatoriSyC(1));
      model.AddPotential(-1, model.GetOnsiteOperatorSzC(1));
      model.AddPotential(-1, model.GetOnsiteOperatorSpC(1));
      model.AddPotential(-1, model.GetOnsiteOperatorSmC(1));
      model.AddPotential(-1, model.GetOnsiteOperatorSCSL(1));
      
      model.AddPotential(-1, model.GetOnsiteOperatorNCTot());
      
      model.AddPotential(-1, model.GetOnsiteOperatorSxL());
      model.AddPotential(-1, model.GetOnsiteOperatoriSyL());
      model.AddPotential(-1, model.GetOnsiteOperatorSzL());
      model.AddPotential(-1, model.GetOnsiteOperatorSpL());
      model.AddPotential(-1, model.GetOnsiteOperatorSmL());
      
      model.AddPotential(-1, model.GetOnsiteOperatorSx());
      model.AddPotential(-1, model.GetOnsiteOperatoriSy());
      model.AddPotential(-1, model.GetOnsiteOperatorSz());
      model.AddPotential(-1, model.GetOnsiteOperatorSp());
      model.AddPotential(-1, model.GetOnsiteOperatorSm());
   };
   
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<float>>(0.5, {2, 2})));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<double>>(0.5, {2, 2})));
   EXPECT_NO_THROW(check_model(model::GeneralModel<model::BaseU1SpinMultiElectrons<long double>>(0.5, {2, 2})));
   
}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_MODEL_GENERAL_MODEL_HPP_ */
