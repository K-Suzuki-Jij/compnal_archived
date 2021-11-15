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
   model.SetSystemSize(4);
   model.SetMagnitudeSpin(0.5);
   compnal::solver::ExactDiag ed(model);
   ed.CalculateGroundState("Lanczos");
   printf("%.30lf\n", ed.GetEigenvalues()[0]);
   printf("%.15lf\n", ed.CalculateCorrelationFunction(model.GetOnsiteOperatorSz(),
                                                      0,
                                                      model.GetOnsiteOperatorSz(),
                                                      1));
   
   printf("%.15lf\n", ed.CalculateCorrelationFunction(model.GetOnsiteOperatorSx(),
                                                      0,
                                                      model.GetOnsiteOperatorSx(),
                                                      1));
}


#endif /* TEST_SOLVER_EXACT_DIAG_HPP_*/
