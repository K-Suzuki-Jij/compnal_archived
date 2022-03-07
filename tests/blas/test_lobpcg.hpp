//
//  test_lobpcg.hpp
//  compnal-cpp
//
//  Created by Kohei Suzuki on 2022/03/08.
//

#ifndef COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_
#define COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_

#include "../../src/blas/lobpcg.hpp"
#include "../../src/blas/compressed_row_storage.hpp"
#include "../../src/blas/braket_vector.hpp"

#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(LOBPCG, Basic) {
   
   blas::CRS<double> m_d({{1.0, 2.0}, {2.0, 1.0}});
   double value = 0.0;
   blas::BraketVector<double> vector;
   blas::EigendecompositionLOBPCG(&value, &vector, m_d);
   
}



} //namespace test
} //namespace compnal

#endif /* COMPNAL_TEST_EIGENDECOMPOSITION_LOBPCG_HPP_ */
