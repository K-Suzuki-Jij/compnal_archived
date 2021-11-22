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

TEST(U1Spin, Basis) {
   const int N = 4;
   compnal::model::U1Spin_1D<double> model(N, 0.5);

   for (int i = 0; i < N - 1; ++i) {
      model.AddInteraction(model.GetOnsiteOperatorSz(), i, model.GetOnsiteOperatorSz(), i+1);
      model.AddInteraction(model.GetOnsiteOperatorSx(), i, model.GetOnsiteOperatorSx(), i+1);
      model.AddInteraction(-1.0, model.GetOnsiteOperatoriSy(), i, model.GetOnsiteOperatoriSy(), i+1);
   }
      
   compnal::solver::ExactDiag ed(model);
   ed.CalculateGroundState("Lanczos");
   printf("%.30lf\n", ed.GetEigenvalues()[0]);

}

#endif /* TEST_SOLVER_EXACT_DIAG_HPP_*/
