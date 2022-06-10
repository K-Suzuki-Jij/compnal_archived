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
//  Created by Kohei Suzuki on 2022/03/08.
//

#ifndef COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_
#define COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_

#include <gtest/gtest.h>

#include "../../src/blas/braket_vector.hpp"
#include "../../src/blas/compressed_row_storage.hpp"
#include "../../src/blas/lobpcg.hpp"

namespace compnal {
namespace test {

TEST(LOBPCG, Basic) {
   blas::CRS<double> m_d({{1.0, 2.0}, {2.0, 1.0}});
   double value = 0.0;
   blas::BraketVector<double> vector;
   blas::EigendecompositionLOBPCG(&value, &vector, m_d);
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_ */
