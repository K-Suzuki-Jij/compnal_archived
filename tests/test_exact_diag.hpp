//
//  test_exact_diag.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef test_exact_diag_hpp
#define test_exact_diag_hpp

#include "solver.hpp"
#include "model.hpp"
#include <gtest/gtest.h>

TEST(ExactDiag, Heisenberg1D) {
   
   compnal::model::Heisenberg1D<double> model(16);
   
   model.SetMagnitudeSpin(0.5);
   
   compnal::solver::ExactDiag1D<compnal::model::Heisenberg1D<double>> exact_diag(model);
   
   exact_diag.Diagonalize();
   
   printf("Dim=%lld, e=%.15lf(%.15lf)\n", exact_diag.GetDim(), exact_diag.GetGSValue(), exact_diag.GetGSValue()/exact_diag.model.GetSystemSize());
   
   
}


#endif /* test_exact_diag_hpp */
