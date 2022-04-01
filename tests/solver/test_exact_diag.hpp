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
//  Created by Kohei Suzuki on 2022/03/26.
//

#ifndef COMPNAL_TEST_SOLVER_EXACT_DIAG_HPP_
#define COMPNAL_TEST_SOLVER_EXACT_DIAG_HPP_

#include "../../src/model/all.hpp"
#include "../../src/solver/exact_diag.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(SolverExactDiag, BaseU1Spin) {
   model::GeneralModel<model::BaseU1Spin<long double>> model;
   const int system_size = 10;
   for (int i = 0; i < system_size - 1; ++i)  {
      model.AddInteraction(i, model.GetOnsiteOperatorSp(), i+1, model.GetOnsiteOperatorSm(), 0.5);
      model.AddInteraction(i, model.GetOnsiteOperatorSm(), i+1, model.GetOnsiteOperatorSp(), 0.5);
      model.AddInteraction(i, model.GetOnsiteOperatorSz(), i+1, model.GetOnsiteOperatorSz(), 1.0);
   }
   solver::ExactDiag solver(model);
   solver.SetTargetQNumber(0);
   solver.SetFlagDisplayInfo(false);
   solver.CalculateGroundState();
   printf("%.32Lf\n", solver.GetEigenvalue(0)/system_size);
   
}



} //namespace testc
} //namespace compnal




#endif /* COMPNAL_TEST_SOLVER_EXACT_DIAG_HPP_ */
