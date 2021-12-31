//
//  test_main.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#include "sparse_matrix/test_compressed_row_storage.hpp"
#include "sparse_matrix/test_eigendecomposition.hpp"
#include "sparse_matrix/test_braket_vector.hpp"
#include "model/test_model_xxz_1d.hpp"
#include "model/test_model_base_u1_spin_multi_electrons.hpp"
#include "solver/test_exact_diag.hpp"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
