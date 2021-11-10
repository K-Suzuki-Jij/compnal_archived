//
//  matrix_vector_operation.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/05.
//

#ifndef COMPNAL_SPARSE_MATRIX_MATRIX_VECTOR_OPERATION_HPP_
#define COMPNAL_SPARSE_MATRIX_MATRIX_VECTOR_OPERATION_HPP_

#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
void CalculateMatrixVectorProduct(BraketVector<RealType> *vector_out,
                                  const RealType coeef,
                                  const CRS<RealType> &matrix_in,
                                  const BraketVector<RealType> &vector_in) {
   
   if (matrix_in.col_dim != vector_in.val.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim  << std::endl;
      ss << "The dimension of the input vector is " << vector_in.val.size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(matrix_in.row_dim);
#pragma omp parallel for
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      RealType temp = 0.0;
      for (std::size_t j = matrix_in.row[i]; j < matrix_in.row[i+1]; ++j) {
         temp += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
      }
      vector_out->val[i] = temp*coeef;
   }
}

template<typename RealType>
void CalculateSymmetricMatrixVectorProduct(BraketVector<RealType> *vector_out,
                                           const RealType coeef,
                                           const CRS<RealType> &matrix_in,
                                           const BraketVector<RealType> &vector_in,
                                           std::vector<std::vector<RealType>> *vectors_work) {
   
   if (matrix_in.row_dim != matrix_in.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not symmetric" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (matrix_in.col_dim != vector_in.val.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is " << matrix_in.col_dim  << std::endl;
      ss << "The dimension of the input vector is " << vector_in.val.size() << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(matrix_in.row_dim);
   
#ifdef _OPENMP
   
   const int num_threads = omp_get_max_threads();
   if (static_cast<int>(vectors_work->size()) != num_threads) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Working vector (vectors_work) must be arrays of the number of parallel threads";
      throw std::runtime_error(ss.str());
   }
   
#pragma omp parallel for schedule (guided)
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      const int      thread_num  = omp_get_thread_num();
      const RealType temp_vec_in = vector_in.val[i];
      RealType       temp_val    = matrix_in.val[matrix_in.row[i + 1] - 1]*temp_vec_in;
      for (std::size_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
         (*vectors_work)[thread_num][matrix_in.col[j]] += matrix_in.val[j]*temp_vec_in;
      }
      (*vectors_work)[thread_num][i] += temp_val;
   }
   
#pragma omp parallel for
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      RealType temp_val = 0.0;
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         temp_val += (*vectors_work)[thread_num][i];
         (*vectors_work)[thread_num][i] = 0.0;
      }
      vector_out->val[i] = temp_val*coeef;
   }
   
#else
   vector_out->Fill(0.0);
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      const RealType temp_vec_in = vector_in.val[i];
      RealType       temp_val    = matrix_in.val[matrix_in.row[i + 1] - 1]*temp_vec_in;
      for (std::size_t j = matrix_in.row[i]; j < matrix_in.row[i + 1] - 1; ++j) {
         temp_val += matrix_in.val[j]*vector_in.val[matrix_in.col[j]];
         vector_out->val[matrix_in.col[j]] += matrix_in.val[j]*temp_vec_in;
      }
      vector_out->val[i] += temp_val*coeef;
   }
   
#endif
}

template<typename RealType>
void CalculateSymmetricMatrixVectorProduct(BraketVector<RealType> *vector_out,
                                           const RealType coeef,
                                           const CRS<RealType> &matrix_in,
                                           const BraketVector<RealType> &vector_in) {
#ifdef _OPENMP
   std::vector<std::vector<RealType>> vector_work(omp_get_max_threads(), std::vector<RealType>(matrix_in.row_dim, 0.0));
   CalculateSymmetricMatrixVectorProduct(vector_out, coeef, matrix_in, vector_in, &vector_work);
#else
   CalculateSymmetricMatrixVectorProduct(vector_out, coeef, matrix_in, vector_in, nullptr);
#endif
}


template<typename RealType>
BraketVector<RealType> CalculateMatrixVectorProduct(const RealType coeef,
                                                    const CRS<RealType> &matrix_in,
                                                    const BraketVector<RealType> &vector_in) {
   BraketVector<RealType> vector_out;
   CalculateMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in);
   return vector_out;
}

template<typename RealType>
BraketVector<RealType> CalculateSymmetricMatrixVectorProduct(const RealType coeef,
                                                             const CRS<RealType> &matrix_in,
                                                             const BraketVector<RealType> &vector_in,
                                                             std::vector<std::vector<RealType>> *vectors_work) {
   BraketVector<RealType> vector_out;
   CalculateSymmetricMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in, vectors_work);
   return vector_out;
}

template<typename RealType>
BraketVector<RealType> CalculateSymmetricMatrixVectorProduct(const RealType coeef,
                                                             const CRS<RealType> &matrix_in,
                                                             const BraketVector<RealType> &vector_in) {
   BraketVector<RealType> vector_out;
   CalculateSymmetricMatrixVectorProduct(&vector_out, coeef, matrix_in, vector_in);
   return vector_out;
}


}
}



#endif /* COMPNAL_SPARSE_MATRIX_MATRIX_VECTOR_OPERATION_HPP_ */
