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
void pybind11ModelIsing(py::module &m, const std::string &post_name = "") {
   
   using Ising = model::Ising<LatticeType, RealType>;
   std::string name = std::string("Ising") + post_name;
   
   auto py_class = py::class_<Ising>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const LatticeType&, const RealType, const RealType>(), "lattice"_a, "interaction_deg_1"_a, "interaction_deg_2"_a);
   
   //Public Member Functions
   py_class.def("set_constant", &Ising::SetConstant, "constant"_a);
   py_class.def("get_interaction", &Ising::GetInteraction);
   py_class.def("get_system_size", &Ising::GetSystemSize);
   py_class.def("get_boundary_condition", &Ising::GetBoundaryCondition);
   py_class.def("get_degree", &Ising::GetDegree);
   py_class.def("calculate_energy", py::overload_cast<const std::vector<typename Ising::OPType>&>(&Ising::CalculateEnergy, py::const_), "sample"_a);
 
   m.def("make_ising", [](const LatticeType &lattice, const RealType interaction_deg_1, const RealType interaction_deg_2) {
      return model::make_ising<LatticeType, RealType>(lattice, interaction_deg_1, interaction_deg_2);
   }, "lattice"_a, "interaction_deg_1"_a, "interaction_deg_2"_a);
}


template<typename RealType>
void pybind11ModelIsingAnyLattice(py::module &m, const std::string &post_name = "") {
   
   using Ising = model::Ising<lattice::AnyLattice, RealType>;
   std::string name = std::string("Ising") + post_name;
   
   using IndexType = typename Ising::IndexType;
   using IndexHash = typename Ising::IndexHash;
   using PairHash  = typename Ising::PairHash;
   
   auto py_class = py::class_<Ising>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const lattice::AnyLattice&, const std::unordered_map<IndexType, RealType, IndexHash>&, const std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash>&>(), "lattice"_a, "linear"_a, "quadratic"_a);
   
   //Public Member Functions
   py_class.def("set_constant", &Ising::SetConstant, "constant"_a);
   py_class.def("generate_index_list", &Ising::GenerateIndexList);
   py_class.def("get_constant", &Ising::GetConstant);
   py_class.def("generate_linear_interaction_as_pair", &Ising::GenerateLinearInteractionAsPair);
   py_class.def("generate_quadratic_interaction_as_pair", &Ising::GenerateQuadraticInteractionAsPair);
   py_class.def("generate_index_list", &Ising::GenerateIndexList);
   py_class.def("get_system_size", &Ising::GetSystemSize);
   py_class.def("get_degree", &Ising::GetDegree);
   py_class.def("get_boundary_condition", &Ising::GetBoundaryCondition);
   py_class.def("calculate_energy", &Ising::CalculateEnergy, "sample"_a);

   m.def("make_ising", [](const lattice::AnyLattice &lattice,
                          const std::unordered_map<IndexType, RealType, IndexHash> &linear,
                          const std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash> &quadratic) {
      return model::make_ising<RealType>(lattice, linear, quadratic);
   }, "lattice"_a, "linear"_a, "quadratic"_a);
}

template<class LatticeType, typename RealType>
void pybind11ModelPolynomialIsing(py::module &m, const std::string &post_name = "") {
   
   using PolyIsing = model::PolynomialIsing<LatticeType, RealType>;
   std::string name = std::string("PolynomialIsing") + post_name;

   auto py_class = py::class_<PolyIsing>(m, name.c_str(), py::module_local());
   
   //Constructors
   py_class.def(py::init<const LatticeType&, const typename PolyIsing::InteractionType&>(), "lattice"_a, "interaction"_a);
   
   //Public Member Functions
   py_class.def("get_system_size", &PolyIsing::GetSystemSize);
   py_class.def("get_boundary_condition", &PolyIsing::GetBoundaryCondition);
   py_class.def("get_degree", &PolyIsing::GetDegree);
   py_class.def("calculate_energy", py::overload_cast<const std::vector<typename PolyIsing::OPType>&>(&PolyIsing::CalculateEnergy, py::const_), "sample"_a);
   
   m.def("make_polynomial_ising", [](const LatticeType &lattice, const typename PolyIsing::InteractionType &interaction) {
      return model::make_polynomial_ising<LatticeType, RealType>(lattice, interaction);
   }, "lattice"_a, "interaction"_a);
}




} // namespace wrapper
} // namespace compnal



#endif /* COMPNAL_PYBIND11_MODEL_HPP_ */
