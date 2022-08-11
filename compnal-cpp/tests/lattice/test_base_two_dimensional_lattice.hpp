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
//  test_two_dimensional_lattice.hpp
//  compnal
//
//  Created by kohei on 2022/08/11.
//  
//

#ifndef COMPNAL_TEST_LATTICE_TWO_DIMENSIONAL_LATTICE_HPP_
#define COMPNAL_TEST_LATTICE_TWO_DIMENSIONAL_LATTICE_HPP_

#include "../../src/lattice/base_two_dimensional_lattice.hpp"
#include "../../src/lattice/square.hpp"
#include "../../src/lattice/triangle.hpp"
#include "../../src/lattice/honeycomb.hpp"
#include <gtest/gtest.h>

namespace compnal {
namespace test {

TEST(LatticeBaseTwoDimensionalLattice, Constructor) {
   
   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{8, 7}.GetSystemSize()), 56);
   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{8, 7}.GetBoundaryCondition()), lattice::BoundaryCondition::OBC);

   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{
      8, 7, lattice::BoundaryCondition::PBC
   }.GetSystemSize()), 56);
   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{
      8, 7, lattice::BoundaryCondition::PBC
   }.GetXSize()), 8);
   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{
      8, 7, lattice::BoundaryCondition::PBC
   }.GetYSize()), 7);
   EXPECT_EQ((lattice::BaseTwoDimensionalLattice{
      8, 7, lattice::BoundaryCondition::PBC
   }.GetBoundaryCondition()), lattice::BoundaryCondition::PBC);

   EXPECT_THROW((lattice::BaseTwoDimensionalLattice{-1, +9}), std::runtime_error);
   EXPECT_THROW((lattice::BaseTwoDimensionalLattice{+3, -1}), std::runtime_error);
   EXPECT_THROW((lattice::BaseTwoDimensionalLattice{-2, -1}), std::runtime_error);

   EXPECT_THROW((lattice::BaseTwoDimensionalLattice{+2, +2, lattice::BoundaryCondition::NONE}), std::runtime_error);
   EXPECT_THROW((lattice::BaseTwoDimensionalLattice{-2, -2, lattice::BoundaryCondition::NONE}), std::runtime_error);

}

TEST(LatticeBaseTwoDimensionalLattice, SetSystemSize) {
   
   lattice::BaseTwoDimensionalLattice square{0, 0};
   square.SetXSize(3);
   square.SetYSize(1);

   EXPECT_EQ(square.GetXSize(), 3);
   EXPECT_EQ(square.GetYSize(), 1);
   
   square.SetSystemSize(2, 4);
   EXPECT_EQ(square.GetXSize(), 2);
   EXPECT_EQ(square.GetYSize(), 4);
   
   EXPECT_THROW(square.SetXSize(-1), std::runtime_error);
   EXPECT_THROW(square.SetYSize(-1), std::runtime_error);
   EXPECT_THROW(square.SetSystemSize(-1, +3), std::runtime_error);
   EXPECT_THROW(square.SetSystemSize(+1, -3), std::runtime_error);
   EXPECT_THROW(square.SetSystemSize(-1, -3), std::runtime_error);

}

TEST(LatticeBaseTwoDimensionalLattice, SetBoundaryCondition) {
   
   lattice::BaseTwoDimensionalLattice square{0, 0, lattice::BoundaryCondition::PBC};
   square.SetBoundaryCondition(lattice::BoundaryCondition::OBC);
   EXPECT_EQ(square.GetBoundaryCondition(), lattice::BoundaryCondition::OBC);
   
   EXPECT_THROW(square.SetBoundaryCondition(lattice::BoundaryCondition::NONE), std::runtime_error);

}

TEST(LatticeSquare, Basic) {
   EXPECT_NE(typeid(lattice::Square), typeid(lattice::BaseTwoDimensionalLattice));
}

TEST(LatticeTriangle, Basic) {
   EXPECT_NE(typeid(lattice::Triangle), typeid(lattice::BaseTwoDimensionalLattice));
}

TEST(LatticeHoneycomb, Basic) {
   EXPECT_NE(typeid(lattice::Honeycomb), typeid(lattice::BaseTwoDimensionalLattice));
}


} // namespace test
} // namespace compnal


#endif /* COMPNAL_TEST_LATTICE_TWO_DIMENSIONAL_LATTICE_HPP_ */
