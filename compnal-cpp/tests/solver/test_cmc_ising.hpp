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
//  test_cmc_ising.hpp
//  compnal
//
//  Created by kohei on 2022/09/25.
//  
//

#ifndef COMPNAL_TEST_SOLVER_CMC_ISING_HPP_
#define COMPNAL_TEST_SOLVER_CMC_ISING_HPP_


#include "../../src/solver/classical_monte_carlo.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {


TEST(SolverClassicalMonteCarlo, IsingChainOBC) {
   const std::int32_t system_size = 10;
   const lattice::Chain lattice(system_size, lattice::BoundaryCondition::OBC);
   model::Ising model(lattice, 0.0, -1.0);
   
   const std::uint64_t seed = 1;
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumThreads(4);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(5);
   solver.SetTemperature(0.01);
   solver.Run(seed);
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      EXPECT_EQ(std::abs(std::accumulate(solver.GetSample(i).begin(), solver.GetSample(i).end(), 0)), system_size);
   }
}

TEST(SolverClassicalMonteCarlo, IsingChainPBC) {
   const std::int32_t system_size = 10;
   const lattice::Chain lattice(system_size, lattice::BoundaryCondition::PBC);
   model::Ising model(lattice, 0.0, -1.0);
   
   const std::uint64_t seed = 1;
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumThreads(4);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(5);
   solver.SetTemperature(0.01);
   solver.Run(seed);
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      EXPECT_EQ(std::abs(std::accumulate(solver.GetSample(i).begin(), solver.GetSample(i).end(), 0)), system_size);
   }
}

TEST(SolverClassicalMonteCarlo, IsingSquareOBC) {
   const std::int32_t x_size = 3;
   const std::int32_t y_size = 4;
   const lattice::Square lattice(x_size, y_size, lattice::BoundaryCondition::OBC);
   model::Ising model(lattice, 0.0, -1.0);
   
   const std::uint64_t seed = 1;
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumThreads(4);
   solver.SetNumSweeps(1000);
   solver.SetNumSamples(8);
   solver.SetTemperature(0.1);
   solver.Run(seed);
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      EXPECT_EQ(std::abs(std::accumulate(solver.GetSample(i).begin(), solver.GetSample(i).end(), 0)), x_size*y_size);
   }
   
}

TEST(SolverClassicalMonteCarlo, IsingSquarePBC) {
   const std::int32_t x_size = 3;
   const std::int32_t y_size = 5;
   const lattice::Square lattice(x_size, y_size, lattice::BoundaryCondition::PBC);
   model::Ising model(lattice, 0.0, -1.0);
   
   const std::uint64_t seed = 1;
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumThreads(4);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(8);
   solver.SetTemperature(0.5);
   solver.Run(seed);
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      EXPECT_EQ(std::abs(std::accumulate(solver.GetSample(i).begin(), solver.GetSample(i).end(), 0)), x_size*y_size);
   }
}



} // namespace test
} // namespace compnal

#endif /* COMPNAL_TEST_SOLVER_CMC_ISING_HPP_ */
