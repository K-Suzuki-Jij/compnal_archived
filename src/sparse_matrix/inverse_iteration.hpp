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
//  Created by Kohei Suzuki on 2021/11/12.
//

#ifndef COMPNAL_SPARSE_MATRIX_INVERSE_ITERATION_HPP_
#define COMPNAL_SPARSE_MATRIX_INVERSE_ITERATION_HPP_

#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
std::pair<int, double> InverseIteration(type::CRS<RealType>          *matrix_in,
                                        type::BraketVector<RealType> *eigenvector,
                                        const RealType         eigenvalue,
                                        const ParametersII     &params = ParametersII(),
                                        const std::vector<type::BraketVector<RealType>> &subspace_vectors = {}
                                        ) {
   
   if (matrix_in->row_dim != matrix_in->col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in->row_dim << ", col=" << matrix_in->col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   const auto start = std::chrono::system_clock::now();
   
   type::BraketVector<RealType> improved_eigenvector;
   type::BraketVector<RealType> vectors_work(matrix_in->row_dim);
   
   if (params.cg.flag_use_initial_vec) {
      if (static_cast<std::int64_t>(eigenvector->val.size()) != matrix_in->row_dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The dimension of the initial vector is not equal to that of the input matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
      improved_eigenvector = *eigenvector;
   }

   matrix_in->DiagonalScaling(params.diag_add - eigenvalue);
   
   for (int step = 0; step < params.max_step; ++step) {
      if (params.cg.flag_symmetric_crs) {
         CalculateSymmetricMatrixVectorProduct(&vectors_work, 1.0, *matrix_in, *eigenvector);
      }
      else {
         CalculateMatrixVectorProduct(&vectors_work, 1.0, *matrix_in, *eigenvector);
      }
      const RealType residual_error = CalculateL1Norm(params.diag_add, *eigenvector, 1.0, vectors_work);
      
      if (params.flag_output_info) {
         std::cout << "\rII_Step[" << step << "]=" << std::scientific << std::setprecision(1);
         std::cout << residual_error << std::flush;
      }
      
      if (residual_error < params.acc) {
         matrix_in->DiagonalScaling(-(params.diag_add - eigenvalue));
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
         std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of inverse iteration:" << time_sec << "[sec]";
         std::cout << " (" << residual_error << ")" << std::flush;
         std::cout << std::endl;
         return {step, time_sec};
      }

      ConjugateGradient(&improved_eigenvector, *matrix_in, *eigenvector, params.cg, subspace_vectors);
      improved_eigenvector.Normalize();
      eigenvector->Assign(improved_eigenvector);
   }
   
   std::stringstream ss;
   ss << "Error in " << __func__ << std::endl;
   ss << "Does not converge" << std::endl;
   throw std::runtime_error(ss.str());
   
}

}
}


#endif /* COMPNAL_SPARSE_MATRIX_INVERSE_ITERATION_HPP_ */
