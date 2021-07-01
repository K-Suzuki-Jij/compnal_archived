//
//  pybind11_main.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/29.
//

#include "pybind11_declare.hpp"


PYBIND11_MODULE(compnal, m) {
   
   py::module_ m_solver      = m.def_submodule("solver");
   py::module_ m_solver_ed1d = m_solver.def_submodule("exact_diag_1d");
   pybindSolverExactDiag1D<compnal::model::Heisenberg1D<double>>(m_solver_ed1d);
   
   py::module_ m_sparse_matrix = m.def_submodule("sparse_matrix");
   pybindCRS<double>(m_sparse_matrix);
   
   
};
