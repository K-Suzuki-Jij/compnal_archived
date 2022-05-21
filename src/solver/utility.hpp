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
//  Created by Kohei Suzuki on 2022/01/29.
//

#ifndef COMPNAL_SOLVER_UTILITY_HPP_
#define COMPNAL_SOLVER_UTILITY_HPP_

namespace compnal {
namespace solver {

enum DiagMethod {

   LANCZOS = 0,
   LOBPCG = 1

};

//! @brief Information for calculating the matrix elements of the Hamiltonian.
template <typename RealType>
struct ExactDiagMatrixComponents {
   //! @brief Values of the matrix elements.
   std::vector<RealType> val;

   //! @brief Column number of the matrix elements.
   std::vector<std::int64_t> basis_affected;

   //! @brief The onsite basis.
   std::vector<int> basis_onsite;

   //! @brief Constants for calculating matrix elements.
   std::vector<std::int64_t> site_constant;

   //! @brief Inverse basis.
   std::unordered_map<std::int64_t, std::int64_t> inv_basis_affected;

   //! @brief Calculation accuracy of the matrix elements.
   double zero_precision = std::pow(10, -15);
};

}  // namespace solver
}  // namespace compnal

#endif /*  COMPNAL_SOLVER_UTILITY_HPP_ */
