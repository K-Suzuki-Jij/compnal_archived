//
//  test_main.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#include "test_exact_diag.hpp"
#include "test_model.hpp"
#include "test_sparse_matrix.hpp"
#include "gtest/gtest.h"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
