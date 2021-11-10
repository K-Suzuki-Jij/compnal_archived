//
//  test_exact_diag.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/07.
//

#ifndef TEST_SOLVER_EXACT_DIAG_HPP_
#define TEST_SOLVER_EXACT_DIAG_HPP_

#include "../../src/solver/all.hpp"
#include "../../src/model/all.hpp"

TEST(XXZ, Basic) {
   compnal::model::XXZ_1D<double> model;
   model.SetSystemSize(18);
   model.SetMagnitudeSpin(0.5);
   compnal::solver::ExactDiag ed(model);
   ed.CalculateGroundState();
   printf("%.30lf\n", ed.GetEigenvalues()[0]);
}


#endif /* TEST_SOLVER_EXACT_DIAG_HPP_*/
