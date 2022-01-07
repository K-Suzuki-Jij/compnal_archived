//
//  test_model_base_u1_spin_multi_electrons.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/31.
//

#ifndef TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_
#define TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_

#include "../../src/model/base_u1_spin_multi_electrons_1d.hpp"
#include <gtest/gtest.h>

TEST(Constructor, Basis) {
   
   compnal::model::BaseU1SpinMultiElectrons_1D<double> model(4, 0.5, {2});
   model.GetOnsiteOperatorCUp(0).PrintMatrix();
   model.GetOnsiteOperatorCUp(0).PrintInfo();
}


#endif /* TEST_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_HPP_ */
