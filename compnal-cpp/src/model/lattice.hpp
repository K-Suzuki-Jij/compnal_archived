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
//  Created by Kohei Suzuki on 2022/06/11.
//

#ifndef COMPNAL_MODEL_LATTICE_HPP_
#define COMPNAL_MODEL_LATTICE_HPP_

#include <stdexcept>

namespace compnal {
namespace model {

//! @brief Enumerated type for lattice.
enum class Lattice {
   
   NONE         , // None Type
   CHAIN        , // One-dimensional chain
   SQUARE       , // Two-dimensional square lattice
   TRIANGLE     , // Two-dimensional triangle lattice
   HONEYCOMB    , // Two-dimensional honeycomb lattice
   CUBIC        , // Three-dimensional cubic lattice
   INFINIT_RANGE, // Infinite range model
   ANY_TYPE       // Users can create any types of lattice.
   
};

enum class BoundaryCondition {
  
   NONE, // None Type
   OBC,  // Open boundary condition
   PBC   // Periodic boundary condition
   
};

BoundaryCondition GetDefaultBoundaryCondition(const Lattice lattice) {
   if (lattice == Lattice::NONE) {
      throw std::runtime_error("Please set Lattice");
   }
   else if (lattice == Lattice::ANY_TYPE) {
      return BoundaryCondition::NONE;
   }
   else {
      return BoundaryCondition::OBC;
   }
}

} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_LATTICE_HPP_ */
