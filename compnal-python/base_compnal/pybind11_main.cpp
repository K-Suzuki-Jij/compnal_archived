//  Copyright 2022 Kohei Suzuki
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Created by Kohei Suzuki on 2022/07/13.
//

#include "src/pybind11_model.cpp"
#include "src/pybind11_solver.cpp"


PYBIND11_MODULE(base_compnal, m) {
   namespace py = pybind11;
   
   using RealType = double;
   
   compnal::wrapper::pybind11ModelLattice(m);
   compnal::wrapper::pybind11SolverUpdater(m);

   
   py::module_ m_model = m.def_submodule("model");
   compnal::wrapper::pybind11ModelPolynomialIsing<RealType>(m_model);
   
   py::module_ m_solver = m.def_submodule("solver");
   compnal::wrapper::pybind11SolverClassicalMonteCarlo<compnal::model::PolynomialIsing<RealType>>(m_solver);

   

};
