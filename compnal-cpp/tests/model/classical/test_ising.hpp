//
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
//  test_ising.hpp
//  compnal
//
//  Created by kohei on 2022/08/28.
//  
//

#ifndef COMPNAL_TEST_MODEL_ISING_HPP_
#define COMPNAL_TEST_MODEL_ISING_HPP_

#include "../../../src/model/classical/ising.hpp"

namespace compnal {
namespace test {

TEST(ModelIsing, AnyLattice) {
   
   auto model = model::make_ising<double>(lattice::AnyLattice{}, {}, {});
   
   
}

TEST(ModelIsing, InfiniteRange) {
   
   lattice::InfiniteRange lattice{10};
   auto model = model::make_ising(lattice, 0.0, 1.0);
   model.GetInteraction();
   
}


} // namespace test
} // namespace compnal

#endif /* COMPNAL_TEST_MODEL_POLYNOMIAL_ISING_HPP_ */
