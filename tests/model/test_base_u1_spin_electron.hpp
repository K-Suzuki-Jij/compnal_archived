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
//  Created by Kohei Suzuki on 2022/03/01.
//

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_ELECTRON_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_ELECTRON_HPP_

#include <gtest/gtest.h>

#include "../../src/model/base_u1_spin_electron.hpp"

namespace compnal {
namespace test {

TEST(ModelBaseU1SpinElectron, Constructors) {
   using RealType = double;
   model::BaseU1SpinElectron<RealType> model;
}

}  // namespace test
}  // namespace compnal

#endif /*  COMPNAL_TEST_MODEL_BASE_U1_SPIN_ELECTRON_HPP_ */
