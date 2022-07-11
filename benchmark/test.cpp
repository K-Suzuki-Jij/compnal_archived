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
//  Created by Kohei Suzuki on 2022/07/12.
//

#include "../src/model/polynomial_ising.hpp"
#include "../src/solver/classical_monte_carlo.hpp"
#include <chrono>

int main(void) {
   
   compnal::model::PolynomialIsing<double> model(1000, {0.0, 0.0, -1.0}, compnal::model::Lattice::INFINIT_RANGE);
   compnal::solver::ClassicalMonteCarlo solver(model, compnal::solver::Updater::METROPOLIS);
   solver.SetNumSweeps(10000);
   solver.SetNumSamples(100);
   solver.SetInverseTemperature(100);
   
   std::chrono::system_clock::time_point start, end;
   
   start = std::chrono::system_clock::now();
   solver.Run();
   end = std::chrono::system_clock::now();
   
   double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
   printf("time %lf[ms]\n", time);
   
   /*
   for (std::size_t i = 0; i < solver.GetSamples().size(); ++i) {
      for (const auto &it: solver.GetSample(i)) {
         printf("%d, ", it);
      }
      printf("\n");
   }
    */
   
   return 0;
}
