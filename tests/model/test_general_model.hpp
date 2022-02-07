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

TEST(ModelGeneralModel, U1Spin) {
   using RealType = double;
   model::GeneralModel<model::BaseU1Spin<RealType>> model;
   model.AddPotential(-1, model.GetOnsiteOperatorSz());
   model.AddPotential("a", model.GetOnsiteOperatorSz());
   model.AddPotential(std::vector<std::variant<int, std::string>>{1, "a"}, model.GetOnsiteOperatorSz());
   
   model.AddInteraction(1, model.GetOnsiteOperatorSp(), 2, model.GetOnsiteOperatorSm());
   model.AddInteraction(1, model.GetOnsiteOperatorSp(), "a", model.GetOnsiteOperatorSm());
   model.AddInteraction(1, model.GetOnsiteOperatorSp(), std::vector<std::variant<int, std::string>>{1, "a"}, model.GetOnsiteOperatorSm());
   
   EXPECT_EQ(model.GetSystemSize(), 5);


}



} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_MODEL_GENERAL_MODEL_HPP_ */
