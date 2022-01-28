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
//  Created by Kohei Suzuki on 2022/01/29.
//

#ifndef COMPNAL_TEST_SOLVER_EXACT_DIAG_GENERAL_MODEL_HPP_
#define COMPNAL_TEST_SOLVER_EXACT_DIAG_GENERAL_MODEL_HPP_

#include "../../src/model/all.hpp"
#include "../../src/solver/all.hpp"
#include "../utility/all.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(ExactDiagGeneralModel, U1Electron) {
   solver::ExactDiag<model::GeneralModel_1D<model::BaseU1Electron_1D<double>>> solver;
   solver.model.SetSystemSize(4);
   solver.CalculateGroundState();
}

TEST(ExactDiagGeneralModel, U1Spin) {
   solver::ExactDiag<model::GeneralModel_1D<model::BaseU1Spin_1D<double>>> solver;
   solver.model.SetSystemSize(4);
   solver.CalculateGroundState();
}

TEST(ExactDiagGeneralModel, U1SpinElectron) {
   solver::ExactDiag<model::GeneralModel_1D<model::BaseU1SpinElectron_1D<double>>> solver;
   solver.model.SetSystemSize(4);
   solver.CalculateGroundState();
}

TEST(ExactDiagGeneralModel, U1SpinMultiElectron) {
   solver::ExactDiag<model::GeneralModel_1D<model::BaseU1SpinMultiElectrons_1D<double>>> solver;
   solver.model.SetSystemSize(2);
   solver.CalculateGroundState();
}

}
}

#endif /* COMPNAL_TEST_SOLVER_EXACT_DIAG_GENERAL_MODEL_HPP_ */
