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
//  square.hpp
//  compnal
//
//  Created by kohei on 2022/08/10.
//  
//

#ifndef COMPNAL_LATTICE_SQUARE_HPP_
#define COMPNAL_LATTICE_SQUARE_HPP_

#include "boundary_condition.hpp"
#include <stdexcept>

namespace compnal {
namespace lattice {

class Square {
   
public:
   
   Square(const int x_size, const int y_size) {
      SetSystemSize(x_size, y_size);
   }
   
   Square(const int x_size, const int y_size, const BoundaryCondition boundary_condition) {
      SetSystemSize(x_size, y_size);
      SetBoundaryCondition(boundary_condition);
   }
   
   void SetXSize(const int x_size) {
      if (x_size < 0) {
         throw std::runtime_error("x_size must be larger than or equal to 0.");
      }
      x_size_ = x_size;
   }
   
   void SetYSize(const int y_size) {
      if (y_size < 0) {
         throw std::runtime_error("y_size must be larger than or equal to 0.");
      }
      y_size_ = y_size;
   }
   
   void SetSystemSize(const int x_size, const int y_size) {
      SetXSize(x_size);
      SetYSize(y_size);
   }
   
   //! @brief Set boundary condition.
   //! @param boundary_condition Boundary condtion.
   void SetBoundaryCondition(const BoundaryCondition boundary_condition) {
      if (boundary_condition == BoundaryCondition::NONE) {
         throw std::runtime_error("BoundaryCondition::NONE cannot be set.");
      }
      bc_ = boundary_condition;
   }
   
   int GetXSize() const {
      return x_size_;
   }
   
   int GetYSize() const {
      return y_size_;
   }
   
   //! @brief Get system size.
   //! @return System size.
   int GetSystemSize() const {
      return x_size_*y_size_;
   }
   
   //! @brief Get boundary condition.
   //! @return Boundary condition.
   BoundaryCondition GetBoundaryCondition() const {
      return bc_;
   }
   
private:
   int x_size_ = 0;
   int y_size_ = 0;
   
   //! @brief Boundary condition.
   BoundaryCondition bc_ = BoundaryCondition::OBC;
};


} // namespace lattice
} // namespace compnal



#endif /* COMPNAL_LATTICE_SQUARE_HPP_ */
