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
#include "../../src/solver/updater/monte_carlo_updater.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(SolverClassicalMonteCarlo, AnyLattice) {
   model::PolynomialIsing<lattice::AnyLattice, double> model{lattice::AnyLattice{}};
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
   
}

TEST(SolverClassicalMonteCarlo, Chain) {
   const lattice::Chain lattice(10, lattice::BoundaryCondition::PBC);
   const std::unordered_map<std::int32_t, double> interaction{{3, -1.0}};
   model::PolynomialIsing model(lattice, interaction);
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
   
}

TEST(SolverClassicalMonteCarlo, Cubic) {
   const lattice::Cubic lattice(3, 3, 3);
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice, interaction);
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
}

TEST(SolverClassicalMonteCarlo, Honeycomb) {
   const lattice::Honeycomb lattice(4, 4);
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice, interaction);
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
}

TEST(SolverClassicalMonteCarlo, Square) {
   const lattice::Square lattice(4, 4);
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice, interaction);
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
}

TEST(SolverClassicalMonteCarlo, Triangle) {
   const lattice::Triangle lattice(4, 4);
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice, interaction);
   solver::ClassicalMonteCarlo solver(model, solver::CMCUpdater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
}

TEST(SolverClassicalMonteCarlo, InfinitRangePolyIsing) {
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice::InfiniteRange{10}, interaction);
   
   solver::ClassicalMonteCarlo solver{model, solver::CMCUpdater::METROPOLIS};
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
   
}

TEST(SolverClassicalMonteCarlo, InfinitRangePolyIsing2) {
   const std::unordered_map<std::int32_t, double> interaction{{3, -0.03}};
   model::PolynomialIsing model(lattice::InfiniteRange{10}, interaction);
   
   solver::ClassicalMonteCarlo solver{model, solver::CMCUpdater::HEAT_BATH};
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(10);
   solver.SetTemperature(0.15);
   solver.Run();
   
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%+d, ", it);
      }
      printf("\n");
   }
   
   printf("%lf\n", solver.CalculateSampleAverage());
   
}


} // namespace test
} // namespace compnal

#endif /* COMPNAL_TEST_SOLVER_CLASSICAL_MONTE_CARLO_HPP_ */
