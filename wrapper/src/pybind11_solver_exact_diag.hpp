//
//  pybind11_solver_exact_diag.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/10.
//

#ifndef COMPNAL_PYBIND11_SOLVER_EXACT_DIAG_HPP_
#define COMPNAL_PYBIND11_SOLVER_EXACT_DIAG_HPP_

#include "../../src/solver/all.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

template<class ModelClass>
void pybind11SolverExactDiag(py::module &m, const std::string &mtype_str) {
   
   using ED = compnal::solver::ExactDiag<ModelClass>;
   auto str = std::string("ExactDiag")+mtype_str;
   
   py::class_<ED>(m, str.c_str(), py::module_local())
   .def(py::init<const ModelClass&>(), "model"_a)
   .def(py::init<const ModelClass&, compnal::sparse_matrix::ParametersAll&>(), "model"_a, "params"_a)
   .def_readonly("model", &ED::model)
   .def_readonly("params", &ED::params)
   .def("get_eigenvectors", &ED::GetEigenvectors)
   .def("get_eigenvalues", &ED::GetEigenvalues)
   .def("calculate_expectation_value", &ED::CalculateExpectationValue, "operator"_a, "site"_a, "level"_a = 0)
   .def("calculate_correlation_function", &ED::CalculateCorrelationFunction, "m_1"_a, "site_1"_a, "m_2"_a, "site_2"_a, "level"_a = 0)
   .def("calculate_ground_state", [](ED &self, const std::string &diag_method) {
      py::scoped_ostream_redirect stream(
          std::cout,                                // std::ostream&
          py::module_::import("sys").attr("stdout") // Python output
      );
      self.CalculateGroundState(diag_method);
   }, "diag_method"_a = "Lanczos");
   
   auto mkci_str = std::string("make_exact_diag_system");
   m.def(mkci_str.c_str(), [](const ModelClass &model) {
      return compnal::solver::ExactDiag<ModelClass>(model);
   }, "model"_a);
   

}




#endif /* COMPNAL_PYBIND11_SOLVER_EXACT_DIAG_HPP_ */
