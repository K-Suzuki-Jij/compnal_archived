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
//  Created by Kohei Suzuki on 2022/01/23.
//

#ifndef COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_
#define COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_

#include "../../src/model/base_u1_spin_multi_electrons_1d.hpp"
#include "../test.hpp"
#include <gtest/gtest.h>


TEST(ModelBaseU1SpinMultiElectrons1D, ConstructorDefault) {
   compnal::model::BaseU1SpinMultiElectrons_1D<double> model;
}


#endif /* COMPNAL_TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_ */
   
