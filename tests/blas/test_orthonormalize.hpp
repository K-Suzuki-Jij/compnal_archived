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
//  Created by Kohei Suzuki on 2022/03/05.
//

#ifndef COMPNAL_TEST_ORTHONORMALIZE_HPP_
#define COMPNAL_TEST_ORTHONORMALIZE_HPP_

#include <gtest/gtest.h>

#include "../../src/blas/braket_vector.hpp"
#include "../../src/blas/orthonormalize.hpp"

namespace compnal {
namespace test {

TEST(Orthonormalize, Type1) {
   auto check = [](auto vectors) {
      blas::Orthonormalize(&vectors);
      for (std::size_t i = 0; i < vectors.size(); ++i) {
         using ElementType = typename std::remove_cvref<decltype(vectors[i])>::type::ValueType;
         EXPECT_TRUE(std::abs(vectors[i] * vectors[i] - 1) < 4 * std::numeric_limits<ElementType>::epsilon());
         for (std::size_t j = 0; j < i; ++j) {
            EXPECT_TRUE(vectors[i] * vectors[j] < 4 * std::numeric_limits<ElementType>::epsilon());
         }
      }
   };

   std::vector<blas::BraketVector<float>> vectors_f = {blas::BraketVector<float>{+0.1F, +0.2F, +0.3F, +0.4F},
                                                       blas::BraketVector<float>{+0.1F, -0.2F, +0.3F, +0.4F},
                                                       blas::BraketVector<float>{+0.1F, +0.2F, -0.3F, +0.4F}};

   std::vector<blas::BraketVector<double>> vectors_d = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                        blas::BraketVector<double>{+0.1, -0.2, +0.3, +0.4},
                                                        blas::BraketVector<double>{+0.1, +0.2, -0.3, +0.4}};

   std::vector<blas::BraketVector<long double>> vectors_ld = {
       blas::BraketVector<long double>{+0.1L, +0.2L, +0.3L, +0.4L},
       blas::BraketVector<long double>{+0.1L, -0.2L, +0.3L, +0.4L},
       blas::BraketVector<long double>{+0.1L, +0.2L, -0.3L, +0.4L}};

   check(vectors_f);
   check(vectors_d);
   check(vectors_ld);

   std::vector<blas::BraketVector<double>> vectors_e1 = {blas::BraketVector<double>{+0.1, +0.2},
                                                         blas::BraketVector<double>{+0.1, -0.2},
                                                         blas::BraketVector<double>{+0.1, +0.3}};

   std::vector<blas::BraketVector<double>> vectors_e2 = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                         blas::BraketVector<double>{+0.1, -0.2, +0.3, +0.4, +0.6},
                                                         blas::BraketVector<double>{+0.1, +0.2, -0.3, +0.4}};

   std::vector<blas::BraketVector<double>> vectors_e3 = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                         blas::BraketVector<double>{+0.0, +0.0, +0.0, +0.0},
                                                         blas::BraketVector<double>{+0.1, +0.2, -0.3, +0.4}};

   std::vector<blas::BraketVector<double>> vectors_e4 = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                         blas::BraketVector<double>{+0.1, +0.2, -0.3, +0.4},
                                                         blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4}};

   // Errors
   EXPECT_THROW(blas::Orthonormalize(&vectors_e1), std::runtime_error);
   EXPECT_THROW(blas::Orthonormalize(&vectors_e2), std::runtime_error);
   EXPECT_THROW(blas::Orthonormalize(&vectors_e3), std::runtime_error);
   EXPECT_THROW(blas::Orthonormalize(&vectors_e4), std::runtime_error);
}

TEST(Orthonormalize, Type2) {
   auto check = [](auto vector, auto vectors) {
      blas::Orthonormalize(&vector, vectors);
      using ElementType = typename std::remove_cvref<decltype(vector)>::type::ValueType;
      EXPECT_TRUE(std::abs(vector * vector - 1) < 4 * std::numeric_limits<ElementType>::epsilon());
      for (std::size_t i = 0; i < vectors.size(); ++i) {
         EXPECT_TRUE(vector * vectors[i] < 4 * std::numeric_limits<ElementType>::epsilon());
      }
   };

   std::vector<blas::BraketVector<float>> vectors_f = {blas::BraketVector<float>{+0.1F, +0.2F, +0.3F, +0.4F},
                                                       blas::BraketVector<float>{+0.1F, -0.2F, +0.3F, +0.4F}};

   std::vector<blas::BraketVector<double>> vectors_d = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                        blas::BraketVector<double>{+0.1, -0.2, +0.3, +0.4}};

   std::vector<blas::BraketVector<long double>> vectors_ld = {
       blas::BraketVector<long double>{+0.1L, +0.2L, +0.3L, +0.4L},
       blas::BraketVector<long double>{+0.1L, -0.2L, +0.3L, +0.4L}};

   blas::Orthonormalize(&vectors_f);
   blas::Orthonormalize(&vectors_d);
   blas::Orthonormalize(&vectors_ld);

   check(blas::BraketVector<float>{+0.1F, +0.2F, -0.3F, +0.4F}, vectors_f);
   check(blas::BraketVector<double>{+0.1, +0.2, -0.3, +0.4}, vectors_d);
   check(blas::BraketVector<long double>{+0.1L, +0.2L, -0.3L, +0.4L}, vectors_ld);

   blas::BraketVector<double> ref_e1({+0.1, +0.3});
   std::vector<blas::BraketVector<double>> vectors_e1 = {blas::BraketVector<double>{+0.1, +0.2},
                                                         blas::BraketVector<double>{+0.1, -0.2}};
   blas::Orthonormalize(&vectors_e1);

   blas::BraketVector<double> ref_e2({+0.1, +0.2, -0.3, +0.4});
   std::vector<blas::BraketVector<double>> vectors_e2 = {blas::BraketVector<double>{+0.1, +0.2, +0.3, +0.4},
                                                         blas::BraketVector<double>{+0.1, -0.2, +0.3, +0.4, +0.6}};

   // Errors
   EXPECT_THROW(blas::Orthonormalize(&ref_e1, vectors_e1), std::runtime_error);
   EXPECT_THROW(blas::Orthonormalize(&ref_e2, vectors_e2), std::runtime_error);
}

}  // namespace test
}  // namespace compnal

#endif /* COMPNAL_TEST_ORTHONORMALIZE_HPP_ */
