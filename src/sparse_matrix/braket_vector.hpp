//
//  braket_vector.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef braket_vector_hpp
#define braket_vector_hpp

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
class BraketVector {
   
public:
   
   explicit BraketVector(const int64_t dim = 0) {
      if (dim < 0) {
         throw std::runtime_error("The size of BraketVector must be larger than or equal to zero");
      }
      val_.resize(dim);
   }
   
   BraketVector(const int64_t dim, const RealType val) {
      if (dim < 0) {
         throw std::runtime_error("The size of BraketVector must be larger than or equal to zero");
      }
      val_.resize(dim);
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         val_[i] = val;
      }
   }
   
   explicit BraketVector(const std::vector<RealType> &val): val_(val) {}
      
   inline double &Val(const int64_t index)       { return val_[index]; }
   inline double  Val(const int64_t index) const { return val_[index]; }
   
   inline void PushVal(const RealType val) { val_.push_back(val); }
   
   inline int64_t GetDim() const { return val_.size(); }
   
   inline void Resize(const int64_t dim) {
      if (dim < 0) {
         throw std::runtime_error("The size of BraketVector must be larger than or equal to zero");
      }
      val_.resize(dim);
   }
   
   void ResetVal(const RealType fill_in = 0.0) {
      const int64_t dim = val_.size();
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         val_[i] = fill_in;
      }
   }
      
   RealType L2Norm() const {
      RealType inner_product = 0.0;
      const int64_t dim = val_.size();
#pragma omp parallel for reduction (+:inner_product)
      for (int64_t i = 0; i < dim; ++i) {
         inner_product += val_[i]*val_[i];
      }
      return std::sqrt(inner_product);
   }
   
   void Normalize(const RealType normalization_factor = 1.0) {
      const RealType coeef = normalization_factor/L2Norm();
      const int64_t dim = val_.size();
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         val_[i] *= coeef;
      }
   }
   
   BraketVector CreateCopy(const RealType coeef = 1.0) const {
      const int64_t dim = val_.size();
      BraketVector vec_out(dim);
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         vec_out.Val(i) = coeef*val_[i];
      }
      return vec_out;
   }
   
   void Print(const std::string display_name = "BraketVector") const {
      for (int64_t i = 0; i < val_.size(); i++) {
         std::cout << display_name << "[" << i << "]=" << val_[i] << std::endl;
      }
   }
   
   BraketVector operator+() const { return *this; }
   BraketVector operator-() const { return CreateCopy(-1.0); }
   BraketVector& operator+=(const BraketVector &vector_rhs) { return *this = CreateVectorSum(*this, vector_rhs); }
   BraketVector& operator-=(const BraketVector &vector_rhs) { return *this = CreateVectorSum(*this, vector_rhs, 1.0, -1.0); }
   
   BraketVector& operator*=(const RealType coeef) {
      MultiplyVectorByScalar(this, coeef);
      return *this;
   }
   
   template<typename CrsRealType>
   BraketVector& operator*=(const CRS<CrsRealType> &matrix_lhs) {
      const auto temp = CreateCopy();
      CreateMatrixVectorProduct(this, matrix_lhs, temp);
      return *this;
   }
   
private:
   std::vector<RealType> val_;
   
};

template<typename RealType>
BraketVector<RealType> operator+(const BraketVector<RealType> &vector_lhs, const BraketVector<RealType> &vector_rhs) {
   return BraketVector<RealType>(vector_lhs) += vector_rhs;
}

template<typename RealType>
BraketVector<RealType> operator-(const BraketVector<RealType> &vector_lhs, const BraketVector<RealType> &vector_rhs) {
   return BraketVector<RealType>(vector_lhs) -= vector_rhs;
}

template<typename RealType>
RealType operator*(const BraketVector<RealType> &vector_lhs, const BraketVector<RealType> &vector_rhs) {
   return CalculateInnerProduct(vector_lhs, vector_rhs);
}

template<typename RealType>
BraketVector<RealType> operator*(const RealType coeef_lhs, const BraketVector<RealType> &vector_rhs) {
   return BraketVector<RealType>(vector_rhs) *= coeef_lhs;
}

template<typename RealType>
BraketVector<RealType> operator*(const BraketVector<RealType> &vector_lhs, const RealType coeef_rhs) {
   return BraketVector<RealType>(vector_lhs) *= coeef_rhs;
}

template<typename RealType>
BraketVector<RealType> operator*(const BraketVector<RealType> &vector_lhs, const CRS<RealType> matrix_rhs) {
   return BraketVector<RealType>(vector_lhs) *= matrix_rhs;
}

template<typename RealType>
BraketVector<RealType> operator*(const CRS<RealType> matrix_lhs, const BraketVector<RealType> &vector_rhs) {
   return BraketVector<RealType>(vector_rhs) *= matrix_lhs;
}

template<typename RealType>
void MultiplyVectorByScalar(BraketVector<RealType> *vector, const RealType coeef) {
   const int64_t dim = vector->GetDim();
#pragma omp parallel for
   for (int64_t i = 0; i < dim; ++i) {
      vector->Val(i) *= coeef;
   }
}

template<typename RealType>
BraketVector<RealType> CreateVectorSum(const BraketVector<RealType> &vector_1, const BraketVector<RealType> &vector_2, const RealType coeef_1 = 1.0, const RealType coeef_2 = 1.0) {
   if (vector_1.GetDim() != vector_2.GetDim()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << vector_1.GetDim() << ", dim_2 = " << vector_2.GetDim() << std::endl;
      throw std::runtime_error(ss.str());
   }
   const int64_t dim = vector_1.GetDim();
   BraketVector<RealType> vector_out(dim);
#pragma omp parallel for
   for (int64_t i = 0; i < dim; ++i) {
      vector_out.Val(i) = coeef_1*vector_1.Val(i) + coeef_2*vector_2.Val(i);
   }
   return vector_out;
}

template<typename RealType>
RealType CalculateInnerProduct(const BraketVector<RealType> &vector_1, const BraketVector<RealType> &vector_2) {
   if (vector_1.GetDim() != vector_2.GetDim()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << vector_1.GetDim() << ", dim_2 = " << vector_2.GetDim() << std::endl;
      throw std::runtime_error(ss.str());
   }
   const int64_t dim = vector_1.GetDim();
   RealType val_out = 0.0;
#pragma omp parallel for reduction (+: val_out)
   for (int64_t i = 0; i < dim; ++i) {
      val_out += vector_1.Val(i)*vector_2.Val(i);
   }
   return val_out;
}

template<typename RealType>
void CreateMatrixVectorProduct(BraketVector<RealType> *vector_out,
                               const CRS<RealType> &matrix_in,
                               const BraketVector<RealType> &vector_in,
                               const RealType coeef = 1.0
                               ) {
   
   if (matrix_in.GetColDim() != vector_in.GetDim()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The column of the input matrix is "    << matrix_in.GetColDim() << std::endl;
      ss << "The dimension of the input vector is " << vector_in.GetDim()    << std::endl;
      ss << "Both must be equal" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   const int64_t row_dim = matrix_in.GetRowDim();
   vector_out->Resize(row_dim);
   
#pragma omp parallel for
   for (int64_t i = 0; i < row_dim; ++i) {
      RealType temp = 0.0;
      const int64_t begin = matrix_in.Row(i);
      const int64_t end   = matrix_in.Row(i+1);
      for (int64_t j = begin; j < end; ++j) {
         temp += matrix_in.Val(j)*vector_in.Val(matrix_in.Col(j));
      }
      vector_out->Val(i) = temp*coeef;
   }
}

template<typename RealType>
BraketVector<RealType> CreateMatrixVectorProduct(const CRS<RealType> &matrix_in,
                                                 const BraketVector<RealType> &vector_in,
                                                 const RealType coeef = 1.0
                                                 ) {
   
   BraketVector<RealType> vector_out;
   CreateMatrixVectorProduct(&vector_out, matrix_in, vector_in, coeef);
   return vector_out;
}

} // namespace sparse_matrix
} // namespace compnal


#endif /* braket_vector_hpp */
