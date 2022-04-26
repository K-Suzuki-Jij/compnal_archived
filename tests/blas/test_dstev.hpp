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
//  Created by Kohei Suzuki on 2022/04/27.
//

#ifndef COMPNAL_TEST_DSTEV_HPP_
#define COMPNAL_TEST_DSTEV_HPP_

#include "../../src/blas/compnal_lapack.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(Dstev, Basic) {
   
   std::vector<double> diag = {1, 2, 3};
   std::vector<double> off_diag = {2, 2};
   double eigenvalue;
   std::vector<double> eigenvector;
   
   blas::Dstev(&eigenvalue, &eigenvector, 2, diag, off_diag);
   
   std::cout << eigenvalue << std::endl;
   for (std::size_t i = 0; i < eigenvector.size(); ++i) {
      printf("vec[%ld]=%lf\n", i, eigenvector[i]);
   }
   
}

} //namespace test
} //namespace compnal



#endif /* COMPNAL_TEST_DSTEV_HPP_ */
