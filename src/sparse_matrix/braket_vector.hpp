//
//  braket_vector.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_
#define COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_

#include <iostream>
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
   
   explicit BraketVector(const std::int64_t dim = 0) {
      this->val.resize(dim);
#pragma omp parallel for
      for (std::int64_t i = 0; i < dim; ++i) {
         this->val[i] = 0.0;
      }
   }
   
   explicit BraketVector(const std::vector<RealType> &vector) {
      Assign(vector);
   }
   
   BraketVector(const BraketVector &vector) {
      Assign(vector);
   }
   
   BraketVector &operator=(const BraketVector &vector) & {
      Assign(vector);
      return *this;
   }
   
   void Fill(const RealType val_in) {
#pragma omp parallel for
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         this->val[i] = val_in;
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
      const RealType norm = L2Norm();
      if (norm == 0.0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "All the elements are zero" << std::endl;
         throw std::runtime_error(ss.str());
      }
      MultiplyByScalar(normalization_factor/norm);
   }
   
   void MultiplyByScalar(const RealType coeef) {
#pragma omp parallel for
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         this->val[i] *= coeef;
      }
   }
   
   RealType L1Norm() const {
      RealType norm = 0.0;
#pragma omp parallel for reduction (+:norm)
      for (std::int64_t i = 0; i < this->val.size(); ++i) {
         norm += std::abs(val[i]);
      }
      return norm;
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
RealType CalculateL1Norm(const RealType coeef_1,
                         const BraketVector<RealType> &braket_vector_1,
                         const RealType coeef_2,
                         const BraketVector<RealType> &braket_vector_2
                         ) {
   
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
      val_out += std::abs(coeef_1*braket_vector_1.val[i] - coeef_2*braket_vector_2.val[i]);
   }
   return val_out;
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
BraketVector<RealType> CalculateVectorSum(const RealType coeef_1,
                                          const BraketVector<RealType> &braket_vector_1,
                                          const RealType coeef_2,
                                          const BraketVector<RealType> &braket_vector_2) {
   BraketVector<RealType> vector_out;
   CalculateVectorSum(&vector_out, coeef_1, braket_vector_1, coeef_2, braket_vector_2);
   return vector_out;
}

} // namespace sparse_matrix
} // namespace compnal


#endif /* COMPNAL_SPARSE_MATRIX_BRAKET_VECTOR_HPP_ */
