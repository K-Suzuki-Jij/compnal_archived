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
//  Created by Kohei Suzuki on 2021/06/14.
//

#include "model/test_base_u1_spin.hpp"
#include "model/test_general_model.hpp"
//#include "model/test_base_u1_electron_1d.hpp"
//#include "model/test_base_u1_spin_electron_1d.hpp"
//#include "model/test_base_u1_spin_multi_electrons_1d.hpp"
//#include "solver/test_exact_diag_general_model.hpp"
#include "utility/test_integer.hpp"
#include "gtest/gtest.h"

int main(int argc, char **argv) {
   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
