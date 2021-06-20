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
   
   compnal::model::Heisenberg1D<double> model(4);
   
   model.SetMagnitudeSpin(0.5);
   
   compnal::solver::ExactDiag1D exact_diag(model);
   
   exact_diag.model.GetOperatorSz().Print();
   
   exact_diag.SetBasis();
   
   exact_diag.PrintBasis();
   
   exact_diag.GenerateHamiltonian().Print();
   
}


#endif /* test_exact_diag_hpp */
