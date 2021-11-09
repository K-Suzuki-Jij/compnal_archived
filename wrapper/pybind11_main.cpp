//
//  pybind11_main.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/29.
//

#include "./src/pybind11_model.hpp"
#include "./src/pybind11_sparse_matrix.hpp"

PYBIND11_MODULE(compnal, m) {
   
   py::module_ m_model = m.def_submodule("model");
   pybindModelXXZ1D<double>(m_model);
   
   py::module_ m_sp_mat = m.def_submodule("sparse_matrix");
   pybindSparseMatrixCRS<double>(m_sp_mat);
   
};
