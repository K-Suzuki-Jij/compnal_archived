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

#ifndef COMPNAL_PYBIND11_MODEL_HPP_
#define COMPNAL_PYBIND11_MODEL_HPP_

#include "../../../cpp/src/model/all.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>


namespace compnal {
namespace wrapper {

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;

void pybind11ModelLattice(py::module &m) {
   
   py::enum_<model::Lattice>(m, "Lattice")
      .value("NONE"         , model::Lattice::NONE         )
      .value("CHAIN"        , model::Lattice::CHAIN        )
      .value("SQUARE"       , model::Lattice::SQUARE       )
      .value("TRIANGLE"     , model::Lattice::TRIANGLE     )
      .value("HONEYCOMB"    , model::Lattice::HONEYCOMB    )
      .value("CUBIC"        , model::Lattice::CUBIC        )
      .value("INFINIT_RANGE", model::Lattice::INFINIT_RANGE)
      .value("ANY_TYPE"     , model::Lattice::ANY_TYPE     );
   
}


template<typename RealType>
void pybind11ModelPolynomialIsing(py::module &m) {
   
   using PolyIsing = model::PolynomialIsing<RealType>;
   
   auto py_class = py::class_<PolyIsing>(m, "PolynomialIsing", py::module_local());
   
   //Constructors
   py_class.def(py::init<const int, const std::vector<RealType>, const model::Lattice>(), "system_size"_a, "interaction"_a, "lattice"_a);
   
   //Public Member Functions
   py_class.def("set_system_size", &PolyIsing::SetSystemSize , "system_size"_a);
   py_class.def("set_interaction", &PolyIsing::SetInteraction, "interaction"_a);
   py_class.def("set_lattice"    , &PolyIsing::SetLattice    , "lattice"_a    );
   py_class.def("get_system_size", &PolyIsing::GetSystemSize );
   py_class.def("get_interaction", &PolyIsing::GetInteraction);
   py_class.def("get_lattice"    , &PolyIsing::GetLattice    );
   
}





} // namespace wrapper
} // namespace compnal



#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
