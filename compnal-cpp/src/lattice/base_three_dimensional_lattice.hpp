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
//  base_three_dimensional_lattice.hpp
//  compnal
//
//  Created by kohei on 2022/08/11.
//  
//

#ifndef COMPNAL_LATTICE_BASE_THREE_DIMENSIONAL_LATTICE_HPP_
#define COMPNAL_LATTICE_BASE_THREE_DIMENSIONAL_LATTICE_HPP_

#include "boundary_condition.hpp"
#include <stdexcept>

namespace compnal {
namespace lattice {

//! @brief Base class to represent the three-dimensional lattice.
class BaseThreeDimensionalLattice {
  
public:
   
   //! @brief Constructor of BaseThreeDimensionalLattice.
   //! @param x_size The size of the x-direction.
   //! @param y_size The size of the y-direction.
   //! @param z_size The size of the z-direction.
   BaseThreeDimensionalLattice(const int x_size, const int y_size, const int z_size) {
      SetSystemSize(x_size, y_size, z_size);
   }
   
   //! @brief Constructor of BaseThreeDimensionalLattice.
   //! @param x_size The size of the x-direction.
   //! @param y_size The size of the y-direction.
   //! @param z_size The size of the z-direction.
   //! @param boundary_condition Boundary condtion. BoundaryCondition::NONE cannot be used here.
   BaseThreeDimensionalLattice(const int x_size, const int y_size, const int z_size, const BoundaryCondition boundary_condition) {
      SetSystemSize(x_size, y_size, z_size);
      SetBoundaryCondition(boundary_condition);
   }
   
   //! @brief Set size of the x-direction.
   //! @param x_size Size of the x-direction.
   void SetXSize(const int x_size) {
      if (x_size < 0) {
         throw std::runtime_error("x_size must be larger than or equal to 0.");
      }
      x_size_ = x_size;
   }
   
   //! @brief Set size of the y-direction.
   //! @param y_size Size of the y-direction.
   void SetYSize(const int y_size) {
      if (y_size < 0) {
         throw std::runtime_error("y_size must be larger than or equal to 0.");
      }
      y_size_ = y_size;
   }
   
   //! @brief Set size of the z-direction.
   //! @param z_size Size of the z-direction.
   void SetZSize(const int z_size) {
      if (z_size < 0) {
         throw std::runtime_error("z_size must be larger than or equal to 0.");
      }
      z_size_ = z_size;
   }
   
   //! @brief Set sizes of the x-, y- and z-direction.
   //! @param x_size Size of the x-direction.
   //! @param y_size Size of the y-direction.
   //! @param z_size Size of the z-direction.
   void SetSystemSize(const int x_size, const int y_size, const int z_size) {
      SetXSize(x_size);
      SetYSize(y_size);
      SetZSize(z_size);
   }
   
   //! @brief Set boundary condition.
   //! @param boundary_condition Boundary condtion. BoundaryCondition::NONE cannot be used here.
   void SetBoundaryCondition(const BoundaryCondition boundary_condition) {
      if (boundary_condition == BoundaryCondition::NONE) {
         throw std::runtime_error("BoundaryCondition::NONE cannot be set.");
      }
      bc_ = boundary_condition;
   }
   
   //! @brief Get size of the x-direction.
   //! @return Size of the x-direction.
   int GetXSize() const {
      return x_size_;
   }
   
   //! @brief Get size of the y-direction.
   //! @return Size of the y-direction.
   int GetYSize() const {
      return y_size_;
   }
   
   //! @brief Get size of the z-direction.
   //! @return Size of the z-direction.
   int GetZSize() const {
      return z_size_;
   }
   
   //! @brief Get system size.
   //! @return System size.
   int GetSystemSize() const {
      return x_size_*y_size_*z_size_;
   }
   
   //! @brief Get boundary condition.
   //! @return Boundary condition.
   BoundaryCondition GetBoundaryCondition() const {
      return bc_;
   }
   
   
private:
   //! @brief Size of the x-direction.
   int x_size_ = 0;
   
   //! @brief Size of the y-direction.
   int y_size_ = 0;
   
   //! @brief Size of the z-direction.
   int z_size_ = 0;
   
   //! @brief Boundary condition.
   BoundaryCondition bc_ = BoundaryCondition::OBC;
   
};



} // namespace lattice
} // namespace compnal


#endif /* COMPNAL_LATTICE_BASE_THREE_DIMENSIONAL_LATTICE_HPP_ */
