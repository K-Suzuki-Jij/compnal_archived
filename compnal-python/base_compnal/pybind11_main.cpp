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

#include "src/pybind11_lattice.hpp"
#include "src/pybind11_model.hpp"
//#include "src/pybind11_solver.hpp"


PYBIND11_MODULE(base_compnal, m) {
   namespace py = pybind11;
   
   using RealType = double;
      
   py::module_ m_lattice = m.def_submodule("lattice");
   compnal::wrapper::pybind11LatticeChain(m_lattice);
   compnal::wrapper::pybind11LatticeSquare(m_lattice);
   compnal::wrapper::pybind11LatticeTriangle(m_lattice);
   compnal::wrapper::pybind11LatticeHoneycomb(m_lattice);
   compnal::wrapper::pybind11LatticeCubic(m_lattice);
   compnal::wrapper::pybind11LatticeInfiniteRange(m_lattice);

   py::module_ m_model = m.def_submodule("model");
   compnal::wrapper::pybind11ModelPolynomialIsing<compnal::lattice::Chain, RealType>(m_model);
   //compnal::wrapper::pybind11ModelPolynomialIsing<compnal::lattice::Square, RealType>(m_model);
   //
   //py::module_ m_solver = m.def_submodule("solver");
   //compnal::wrapper::pybind11SolverUpdater(m_solver);
   //compnal::wrapper::pybind11SolverClassicalMonteCarlo<compnal::model::PolynomialIsing<RealType>>(m_solver);



};
