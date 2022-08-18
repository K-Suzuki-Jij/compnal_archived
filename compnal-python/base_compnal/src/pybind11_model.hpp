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

#include "../../../compnal-cpp/src/model/all.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>


namespace compnal {
namespace wrapper {

namespace py = pybind11;

//The following does not bring in anything else from the pybind11 namespace except for literals.
using namespace pybind11::literals;


template<class LatticeType, typename RealType>
void pybind11ModelPolynomialIsing(py::module &m, const std::string &post_name = "") {
   
   using PolyIsing = model::PolynomialIsing<LatticeType, RealType>;
   std::string name = std::string("PolynomialIsing") + post_name;

   auto py_class = py::class_<PolyIsing>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const LatticeType&, const std::unordered_map<std::int32_t, RealType>&>(), "lattice"_a, "interaction"_a);
   
   //Public Member Functions
   py_class.def("set_interaction", &PolyIsing::SetInteraction, "degree"_a, "value"_a);
   py_class.def("get_interaction", &PolyIsing::GetInteraction);
   py_class.def("get_system_size", &PolyIsing::GetSystemSize);
   py_class.def("get_degree", &PolyIsing::GetDegree);
   
}

template<typename RealType>
void pybind11ModelPolynomialIsingAnyLattice(py::module &m, const std::string &post_name = "") {
   
   using PolyIsing = model::PolynomialIsing<lattice::AnyLattice, RealType>;
   std::string name = std::string("PolynomialIsing") + post_name;

   auto py_class = py::class_<PolyIsing>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const lattice::AnyLattice&>(), "lattice"_a);
   
   //Public Member Functions
   py_class.def("add_interaction", &PolyIsing::AddInteraction, "index_list"_a, "value"_a);
   py_class.def("get_interaction", &PolyIsing::GetInteraction);
   py_class.def("get_index_set", &PolyIsing::GetIndexSet);
   py_class.def("generate_index_list", &PolyIsing::GenerateIndexList);
   py_class.def("get_system_size", &PolyIsing::GetSystemSize);
   
}





} // namespace wrapper
} // namespace compnal



#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
