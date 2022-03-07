//
//  test_eigendecomposition_lanczos.hpp
//  compnal-cpp
//
//  Created by Kohei Suzuki on 2022/03/08.
//

#ifndef COMPNAL_TEST_EIGENDECOMPOSITION_LANCZOS_HPP_
#define COMPNAL_TEST_EIGENDECOMPOSITION_LANCZOS_HPP_

#include "../../src/blas/lanczos.hpp"
#include "../../src/blas/compressed_row_storage.hpp"
#include "../../src/blas/braket_vector.hpp"

#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(Lanczos, Basic) {
   blas::CRS<double> m_d({{1.0, 2.0}, {2.0, 1.0}});
   double value = 0.0;
   blas::BraketVector<double> vector;
   blas::EigendecompositionLanczos(&value, &vector, m_d);
   
   
   
   
   
   
   
   
}



} //namespace test
} //namespace compnal


#endif /* COMPNAL_TEST_EIGENDECOMPOSITION_LANCZOS_HPP_ */
