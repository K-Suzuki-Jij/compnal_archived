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

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
struct BraketVector {
   
   std::vector<RealType> val;
   
   explicit BraketVector(const std::size_t dim = 0) {
      this->val.resize(dim);
   }
   
   BraketVector(const BraketVector &vector) {
      this->val.resize(vector.val.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.val.size(); ++i) {
         this->val[i] = vector.val[i];
      }
   }
   
   BraketVector &operator=(const BraketVector &vector) & {
      this->val.resize(vector.val.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.val.size(); ++i) {
         this->val[i] = vector.val[i];
      }
      return *this;
   }
   
   void Free() {
      std::vector<RealType>().swap(val);
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
      ss << "dim_1 = " << braket_vector_1.val.size() << ", dim_2 = " << braket_vector_2.size() << std::endl;
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
      ss << "dim_1 = " << braket_vector_1.val.size() << ", dim_2 = " << braket_vector_2.size() << std::endl;
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
void CreateMatrixVectorProduct(BraketVector<RealType> *vector_out,
                               const CRS<RealType> &matrix_in,
                               const BraketVector<RealType> &vector_in,
                               const RealType coeef = 1.0
                               ) {
   
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
         temp += matrix_in.val[j]*vector_in.val(matrix_in.col[j]);
      }
      vector_out->val[i] = temp*coeef;
   }
}



} // namespace sparse_matrix
} // namespace compnal


#endif /* COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_ */
