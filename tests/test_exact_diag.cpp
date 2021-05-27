//
//  test_exact_diag.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#include "solver.hpp"
#include "model.hpp"
#include <gtest/gtest.h>

TEST(ExactDiag, Constructor) {
   
   compnal::model::Heisenberg1D<double> model(3, 0.5);
   
   compnal::solver::ExactDiag1D exact_diag(model);
   
   exact_diag.model.GetOperatorSz().Print();
   
   exact_diag.CalculateExpectationValue(exact_diag.model.GetOperatorSz(), 0);
}
