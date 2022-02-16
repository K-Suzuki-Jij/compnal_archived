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
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef COMPNAL_TYPE_BRAKET_VECTOR_HPP_
#define COMPNAL_TYPE_BRAKET_VECTOR_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace type {

template<typename ElementType>
class BraketVector {
   
public:
   using ValueType = ElementType;
   
   
   BraketVector() {};
   
   explicit BraketVector(const std::vector<ElementType> &vector) {
      Assign(vector);
   }
   
   BraketVector(const BraketVector &vector) {
      Assign(vector);
   }
   
   BraketVector &operator=(const BraketVector &vector) & {
      Assign(vector);
      return *this;
   }

   template<typename T>
   void Fill(const T value) {
#pragma omp parallel for
      for (std::size_t i = 0; i < value_.size(); ++i) {
         value_[i] = value;
      }
   }
   
   void Free() {
      std::vector<ElementType>().swap(value_);
   }
   
   void Clear() {
      value_.clear();
   }
   
   void Assign(const BraketVector &vector) {
      value_.resize(vector.Size());
#pragma omp parallel for
      for (std::size_t i = 0; i < value_.size(); ++i) {
         value_[i] = vector.Val(i);
      }
   }
   
   template<typename T>
   void Assign(const std::vector<T> &vector) {
      value_.resize(vector.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.size(); ++i) {
         value_[i] = vector[i];
      }
   }
   
   template<typename T>
   void Normalize(const T normalization_factor = 1) {
      const auto norm = L2Norm();
      if (norm == 0.0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "All the elements are zero" << std::endl;
         throw std::runtime_error(ss.str());
      }
      MultiplyByScalar(normalization_factor/norm);
   }
   
   template<typename T>
   void MultiplyByScalar(const T coeef) {
#pragma omp parallel for
      for (std::size_t i = 0; i < value_.size(); ++i) {
         value_[i] *= coeef;
      }
   }
   
   auto L1Norm() const {
      if constexpr (std::is_floating_point<ElementType>::value) {
         return L1Norm<ElementType>(value_);
      }
      else {
         return L1Norm<double>(value_);
      }
   }
   
   auto L2Norm() const {
      if constexpr (std::is_floating_point<ElementType>::value) {
         return L2Norm<ElementType>(value_);
      }
      else {
         return L2Norm<double>(value_);
      }
   }
   
   void Print(const std::string display_name = "BraketVector") const {
      for (std::size_t i = 0; i < value_.size(); i++) {
         std::cout << display_name << "[" << i << "]=" << value_[i] << std::endl;
      }
   }
   
   inline std::int64_t Size() const { return static_cast<std::int64_t>(value_.size()); }
   
   template<typename T>
   inline const ElementType &Val(const T i) const { return value_[i]; }
   
   template<typename T>
   inline ElementType &Val(const T i) { return value_[i]; }
   
   template<typename T>
   inline void Resize(const T n) {
      value_.resize(n);
   }
   
private:
   std::vector<ElementType> value_;
   
   
   template<typename ReturnType, typename T>
   ReturnType L1Norm(const std::vector<T> &vec) const {
      ReturnType norm = 0.0;
#pragma omp parallel for reduction (+:norm)
      for (std::size_t i = 0; i < vec.size(); ++i) {
         norm += std::abs(vec[i]);
      }
      return norm;
   }
   
   template<typename ReturnType, typename T>
   ReturnType L2Norm(const std::vector<T> &vec) const {
      ReturnType inner_product = 0.0;
#pragma omp parallel for reduction (+:inner_product)
      for (std::size_t i = 0; i < vec.size(); ++i) {
         inner_product += vec[i]*vec[i];
      }
      return std::sqrt(inner_product);
   }
   
};

//------------------------------------------------------------------
//----------------------Operator overloading------------------------
//------------------------------------------------------------------
//! @brief Operator overloading: addition operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
BraketVector<decltype(T1{0}+T2{0})> operator+(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   return CalculateVectorVectorSum(T1{1}, lhs, T2{1}, rhs);
}

//! @brief Operator overloading: subtraction operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
BraketVector<decltype(T1{0}-T2{0})> operator-(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   return CalculateVectorVectorSum(T1{1}, lhs, T2{-1}, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
decltype(T1{0}*T2{0}) operator*(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   return CalculateVectorVectorProduct(lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
BraketVector<decltype(T1{0}*T2{0})> operator*(const T1 lhs, const BraketVector<T2> &rhs) {
   return CalculateScalarVectorProduct(lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
BraketVector<decltype(T1{0}*T2{0})> operator*(const BraketVector<T1> &lhs, const T2 rhs) {
   return CalculateScalarVectorProduct(rhs, lhs);
}

//! @brief Operator overloading: equality operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
bool operator==(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   if (lhs.Size() != rhs.Size()) {
      return false;
   }
   for (std::int64_t i = 0; i < lhs.Size(); ++i) {
      if (lhs.Val(i) != rhs.Val(i)) {
         return false;
      }
   }
   return true;
}

//! @brief Operator overloading: output operator.
//! @tparam T Type of BraketVector.
//! @param os Ostream object.
//! @param v BraketVector.
template<typename T>
std::ostream& operator<<(std::ostream &os, const BraketVector<T> &v) {
   os << std::fixed;
   os << std::setprecision(std::numeric_limits<T>::max_digits10);
   for (std::int64_t i = 0; i < v.Size(); ++i) {
      os << "v[";
      os << std::noshowpos << std::left << std::setw(3) << i << "]=";
      os << std::showpos << v.Val(i) << std::endl;
   }
   return os;
}


template<typename T1, typename T2, typename T3, typename T4>
BraketVector<decltype(T1{0}*T2{0}+T3{0}*T4{0})> CalculateVectorVectorSum(const T1 coeef_1,
                                                                         const BraketVector<T2> &braket_vector_1,
                                                                         const T3 coeef_2,
                                                                         const BraketVector<T4> &braket_vector_2) {
   if (braket_vector_1.Size() != braket_vector_2.Size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.Size() << ", dim_2 = " << braket_vector_2.Size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   BraketVector<decltype(T1{0}*T2{0}+T3{0}*T4{0})> vector_out;
   vector_out.Resize(braket_vector_1.Size());
#pragma omp parallel for
   for (std::int64_t i = 0; i < braket_vector_1.Size(); ++i) {
      vector_out.Val(i) = coeef_1*braket_vector_1.Val(i) + coeef_2*braket_vector_2.Val(i);
   }
   return vector_out;
}

template<typename T1, typename T2>
decltype(T1{0}*T2{0}) CalculateVectorVectorProduct(const BraketVector<T1> &braket_vector_1,
                                                   const BraketVector<T2> &braket_vector_2) {
   if (braket_vector_1.Size() != braket_vector_2.Size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.Size() << ", dim_2 = " << braket_vector_2.Size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   decltype(T1{0}*T2{0}) val_out = 0.0;
#pragma omp parallel for reduction (+: val_out)
   for (std::int64_t i = 0; i < braket_vector_1.Size(); ++i) {
      val_out += braket_vector_1.Val(i)*braket_vector_2.Val(i);
   }
   return val_out;
}

template<typename T1, typename T2>
BraketVector<decltype(T1{0}*T2{0})> CalculateScalarVectorProduct(const T1 value,
                                                                 const BraketVector<T2> &braket_vector) {
   
   BraketVector<decltype(T1{0}*T2{0})> out;
   out.Resize(braket_vector.Size());
#pragma omp parallel for
   for (std::int64_t i = 0; i < braket_vector.Size(); ++i) {
      out.Val(i) = value*braket_vector.Val(i);
   }
   return out;
}

template<typename T1, typename T2, typename T3, typename T4>
decltype(T1{0}*T2{0}-T3{0}*T4{0}) CalculateL1Norm(const T1 coeef_1,
                                                  const BraketVector<T2> &braket_vector_1,
                                                  const T3 coeef_2,
                                                  const BraketVector<T4> &braket_vector_2
                                                  ) {
   
   if (braket_vector_1.Size() != braket_vector_2.Size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.Size() << ", dim_2 = " << braket_vector_2.Size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   decltype(T1{0}*T2{0}-T3{0}*T4{0}) val_out = 0.0;
#pragma omp parallel for reduction (+: val_out)
   for (std::size_t i = 0; i < braket_vector_1.Size(); ++i) {
      val_out += std::abs(coeef_1*braket_vector_1.Val(i) - coeef_2*braket_vector_2.Val(i));
   }
   return val_out;
}

} // namespace type
} // namespace compnal


#endif /* COMPNAL_TYPE_BRAKET_VECTOR_HPP_ */
