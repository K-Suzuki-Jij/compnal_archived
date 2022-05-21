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
//  Created by Kohei Suzuki on 2022/03/03.
//

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_

#include <gtest/gtest.h>

#include "../../src/model/base_u1_spin_multi_electrons.hpp"

namespace compnal {
namespace test {

TEST(BaseU1SpinMultiElectrons, Constructors) {
   using RealType = double;
   model::BaseU1SpinMultiElectrons<RealType> model;
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_ */
