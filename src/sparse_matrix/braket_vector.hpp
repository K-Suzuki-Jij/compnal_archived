//
//  braket_vector.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_
#define COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_

#include "compressed_row_storage.hpp"
#include <sstream>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
struct BraketVector {
   
   std::vector<RealType> val;
   
   explicit BraketVector(const std::size_t dim = 0) {
      this->val.resize(dim);
#pragma omp parallel for
      for (std::size_t i = 0; i < dim; ++i) {
         this->val[i] = 0.0;
      }
   }
   
   explicit BraketVector(const std::vector<RealType> &vector) {
      Assign(vector);
   }
   
   explicit BraketVector(const BraketVector &vector) {
      Assign(vector);
   }
   
   BraketVector &operator=(const BraketVector &vector) & {
      Assign(vector);
      return *this;
   }
   
   void Fill(const RealType val) {
#pragma omp parallel for
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         this->val[i] = val;
      }
   }
   
   void Free() {
      std::vector<RealType>().swap(val);
   }
   
   void Clear() {
      val.clear();
   }
   
   void Assign(const BraketVector &vector) {
      this->val.resize(vector.val.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.val.size(); ++i) {
         this->val[i] = vector.val[i];
      }
   }
   
   void Assign(const std::vector<RealType> &vector) {
      this->val.resize(vector.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.size(); ++i) {
         this->val[i] = vector[i];
      }
   }
   
   void Normalize(const RealType normalization_factor = 1.0) {
      MultiplyByScalar(normalization_factor/L2Norm());
   }
   
   void MultiplyByScalar(const RealType coeef) {
   #pragma omp parallel for
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         this->val[i] *= coeef;
      }
   }
      
   RealType L2Norm() const {
      RealType inner_product = 0.0;
#pragma omp parallel for reduction (+:inner_product)
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         inner_product += val[i]*val[i];
      }
      return std::sqrt(inner_product);
   }
   
   void Print(const std::string display_name = "BraketVector") const {
      for (std::size_t i = 0; i < this->val.size(); i++) {
         std::cout << display_name << "[" << i << "]=" << this->val[i] << std::endl;
      }
   }
   
};

template<typename RealType>
void CalculateVectorSum(BraketVector<RealType> *vector_out,
                        const RealType coeef_1,
                        const BraketVector<RealType> &braket_vector_1,
                        const RealType coeef_2,
                        const BraketVector<RealType> &braket_vector_2) {
   if (braket_vector_1.val.size() != braket_vector_2.val.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.val.size() << ", dim_2 = " << braket_vector_2.val.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   vector_out->val.resize(braket_vector_1.val.size());
#pragma omp parallel for
   for (std::size_t i = 0; i < braket_vector_1.val.size(); ++i) {
      vector_out->val[i] = coeef_1*braket_vector_1.val[i] + coeef_2*braket_vector_2.val[i];
   }
}

template<typename RealType>
RealType CalculateInnerProduct(const BraketVector<RealType> &braket_vector_1,
                               const BraketVector<RealType> &braket_vector_2) {
   if (braket_vector_1.val.size() != braket_vector_2.val.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.val.size() << ", dim_2 = " << braket_vector_2.val.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   RealType val_out = 0.0;
#pragma omp parallel for reduction (+: val_out)
   for (std::size_t i = 0; i < braket_vector_1.val.size(); ++i) {
      val_out += braket_vector_1.val[i]*braket_vector_2.val[i];
   }
   return val_out;
}

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
                                           std::vector<std::vector<RealType>> *vectors_work = nullptr) {
   
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
   if (vectors_work->size() != num_threads) {
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



} // namespace sparse_matrix
} // namespace compnal


#endif /* COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_ */
