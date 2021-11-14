//
//  inverse_iteration.hpp
//  compnal
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
std::pair<int, double> InverseIteration(CRS<RealType>          *matrix_in,
                                        BraketVector<RealType> *eigenvector,
                                        const RealType         eigenvalue,
                                        const ParametersII     &params = ParametersII()
                                        ) {
   
   if (matrix_in->row_dim != matrix_in->col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in->row_dim << ", col=" << matrix_in->col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   const auto start = std::chrono::system_clock::now();
   
   BraketVector<RealType> improved_eigenvector;
   BraketVector<RealType> vectors_work(matrix_in->row_dim);
   
   if (params.cg.flag_use_initial_vec) {
      if (eigenvector->val.size() != matrix_in->row_dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The dimension of the initial vector is not equal to that of the input matrix." << std::endl;
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
         std::cout << "\rII_Step[" << step + 1 << "]=" << std::scientific << std::setprecision(1);
         std::cout << residual_error << std::flush;
      }
      
      if (residual_error < params.acc) {
         matrix_in->DiagonalScaling(-(params.diag_add - eigenvalue));
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/TIME_UNIT_CONSTANT;
         std::cout << std::defaultfloat << std::setprecision(8) << "\rElapsed time of inver iteration:" << time_sec << "[sec]" << std::flush;
         std::cout << std::endl;
         return {step, time_sec};
      }

      ConjugateGradient(&improved_eigenvector, *matrix_in, *eigenvector, params.cg);
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
