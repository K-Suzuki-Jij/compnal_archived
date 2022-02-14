//
//  test_braket_vector.hpp.h
//  compnal-cpp
//
//  Created by Kohei Suzuki on 2022/02/15.
//

#ifndef COMPNAL_TEST_BRAKET_VECTOR_HPP_
#define COMPNAL_TEST_BRAKET_VECTOR_HPP_

#include "../../src/type/braket_vector.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(BraketVector, Addition) {
   using type::BraketVector;
   EXPECT_EQ(BraketVector<long double>({1.0}) + BraketVector<long double>({1.0}), BraketVector<long double>({2.0}));
   
}

}
}

#endif /* COMPNAL_TEST_BRAKET_VECTOR_HPP_ */
