//
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
//  infinite_range.hpp
//  compnal
//
//  Created by kohei on 2022/08/10.
//  
//

#ifndef COMPNAL_LATTICE_INFINITE_RANGE_HPP_
#define COMPNAL_LATTICE_INFINITE_RANGE_HPP_

#include "boundary_condition.hpp"
#include <stdexcept>

namespace compnal {
namespace lattice {

//! @brief Class to represent the infinite-dimensional lattice.
class InfiniteRange {
   
public:
   
   //! @brief Constructor of InfiniteRange class.
   //! @param system_size System size.
   InfiniteRange(const int system_size) {
      SetSystemSize(system_size);
   }
   
   //! @brief Set system size.
   //! @param system_size System size.
   void SetSystemSize(const int system_size) {
      if (system_size < 0) {
         throw std::runtime_error("system_size must be larger than or equal to 0.");
      }
      system_size_ = system_size;
   }
   
   //! @brief Get system size.
   //! @return System size.
   int GetSystemSize() const {
      return system_size_;
   }
   
   //! @brief Get boundary condition.
   //! @return Boundary condition, BoundaryCondition::NONE.
   BoundaryCondition GetBoundaryCondition() const {
      return BoundaryCondition::NONE;
   }
   
private:
   //! @brief System size.
   int system_size_ = 0;
   
};

} // namespace lattice
} // namespace compnal


#endif /* COMPNAL_LATTICE_INFINITE_RANGE_HPP_ */
