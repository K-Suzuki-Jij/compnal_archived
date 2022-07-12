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
//  Created by Kohei Suzuki on 2022/06/11.
//

#ifndef COMPNAL_TEST_SOLVER_CLASSICAL_MONTE_CARLO_HPP_
#define COMPNAL_TEST_SOLVER_CLASSICAL_MONTE_CARLO_HPP_

#include "../../src/solver/classical_monte_carlo.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(SolverClassicalMonteCarlo, InfinitRangePolyIsing) {
   
   model::PolynomialIsing<double> model(10, {0.0, -1.0}, model::Lattice::INFINIT_RANGE);
   solver::ClassicalMonteCarlo solver(model, solver::Updater::METROPOLIS);
   solver.SetNumSweeps(1000);
   solver.SetNumSamples(10);
   solver.SetInverseTemperature(100);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%d, ", it);
      }
      printf("\n");
   }
   
}


} // namespace test
} // namespace compnal

#endif /* COMPNAL_TEST_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
