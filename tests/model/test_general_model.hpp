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

   check_model(model::GeneralModel<model::BaseU1Spin<double>>());
   check_model(model::GeneralModel<model::BaseU1Electron<double>>());

}

} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_MODEL_GENERAL_MODEL_HPP_ */
