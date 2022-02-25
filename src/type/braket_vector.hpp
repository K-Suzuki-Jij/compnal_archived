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

namespace compnal {
namespace type {

//! @brief The wrapper class of std::vector to represent braket.
//! Note that there is no difference of <bra| and |ket>,
//! since the template parameter ElementType cannot accept Complex type.
//! @tparam ElementType The value type of the vector elements.
template<typename ElementType>
class BraketVector {
   
public:
   //------------------------------------------------------------------
   //------------------------Public Type Alias-------------------------
   //------------------------------------------------------------------
   //! @brief Alias of ElementType.
   using ValueType = ElementType;
   
   
   //------------------------------------------------------------------
   //----------------------Public Member Variables---------------------
   //------------------------------------------------------------------
   //! @brief The elements of BraketVector
   std::vector<ElementType> value_list;
   
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of BraketVector class.
   BraketVector() {};
   
   //! @brief Constructor of BraketVector class.
   //! @tparam T Value type of std::vector.
   //! @param vector The vector elements as std::vector.
   template<typename T=ElementType>
   explicit BraketVector(const std::vector<T> &vector) {
      Assign(vector);
   }
   
   //! @brief Copy constructor of BraketVector class.
   //! @tparam T Value type of the BraketVector.
   //! @param vector The vector elements as BraketVector.
   template<typename T>
   BraketVector(const BraketVector<T> &vector) {
      Assign(vector);
   }
   
   
   //------------------------------------------------------------------
   //----------------------Public Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Fill the value to BraketVector.
   //! @tparam T Value type.
   //! @param value To be filled in.
   template<typename T>
   void Fill(const T value) {
#pragma omp parallel for
      for (std::size_t i = 0; i < this->value_list.size(); ++i) {
         this->value_list[i] = value;
      }
   }
   
   //! @brief Clear the elements and free memory.
   void Free() {
      std::vector<ElementType>().swap(this->value_list);
   }
   
   
   //! @brief Assign BraketVector
   //! @tparam T Value type of the BraketVector.
   //! @param vector The vector elements as BraketVector.
   template<typename T>
   void Assign(const BraketVector<T> &vector) {
      this->value_list.resize(vector.value_list.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < this->value_list.size(); ++i) {
         this->value_list[i] = vector.value_list[i];
      }
   }
   
   //! @brief Assign BraketVector
   //! @tparam T Value type of std::vector.
   //! @param vector The vector elements as std::vector.
   template<typename T>
   void Assign(const std::vector<T> &vector) {
      this->value_list.resize(vector.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < vector.size(); ++i) {
         this->value_list[i] = vector[i];
      }
   }
   
   //! @brief Noemalize BraketVector.
   //! @tparam T Value type.
   //! @param normalization_factor Normalized as \f$\|\boldsymbol{v}\|={\rm normalization\_factor} \f$.
   template<typename T>
   void Normalize(const T normalization_factor = T{1.0}) {
      const auto norm = CalculateL2Norm();
      if (norm == 0.0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "All the elements are zero" << std::endl;
         throw std::runtime_error(ss.str());
      }
      MultiplyByScalar(normalization_factor/norm);
   }
   
   //! @brief Multiply by scalar to the elements.
   //! @tparam T Value type.
   //! @param coeff The value to be multiplied.
   template<typename T>
   void MultiplyByScalar(const T coeff) {
#pragma omp parallel for
      for (std::size_t i = 0; i < this->value_list.size(); ++i) {
         this->value_list[i] *= coeff;
      }
   }
   
   //! @brief Calculate L1 norm.
   //! @return L1 norm.
   auto CalculateL1Norm() const {
      if constexpr (std::is_floating_point<ElementType>::value) {
         return CalculateL1Norm<ElementType>();
      }
      else {
         return CalculateL1Norm<double>();
      }
   }
   
   //! @brief Calculate L2 norm.
   //! @return L2 norm.
   auto CalculateL2Norm() const {
      if constexpr (std::is_floating_point<ElementType>::value) {
         return CalculateL2Norm<ElementType>();
      }
      else {
         return CalculateL2Norm<double>();
      }
   }
   
   //! @brief Print BraketVector elements.
   void Print(const std::string display_name = "BraketVector") const {
      for (std::size_t i = 0; i < this->value_list.size(); i++) {
         std::cout << display_name << "[" << i << "]=" << this->value_list[i] << std::endl;
      }
   }
   
   //! @brief Resize BraketVector
   //! @tparam T Value type.
   //! @param n BraketVector will be resized to n.
   template<typename T>
   void Resize(const T n) { value_list.resize(n); }
   
   //! @brief Return the BraketVector size as std::int64_t.
   //! @return The size.
   std::int64_t Size() const { return static_cast<std::int64_t>(value_list.size()); }
   
   //------------------------------------------------------------------
   //-----------------------Operator Overloading-----------------------
   //------------------------------------------------------------------
   //! @brief Operator overloading: unary plus operator.
   //! @return BraketVector object, \f$ +\boldsymbol{v} \f$.
   BraketVector operator+() const {
      return *this;
   }
   
   //! @brief Operator overloading: unary negation operator.
   //! @return BraketVector object, \f$ -\boldsymbol{v} \f$.
   BraketVector operator-() const {
      return CalculateScalarVectorProduct(-1, *this);
   }
   
   //! @brief Operator overloading: compound assignment plus operator.
   //! @tparam T Value type of the right-hand side BraketVector.
   //! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
   //! @return BraketVector object, \f$ \boldsymbol{v} + \boldsymbol{v}_{\rm rhs}\f$.
   template<typename T>
   BraketVector& operator+=(const BraketVector<T> &rhs) {
      return *this = *this + rhs;
   }
   
   //! @brief Operator overloading: compound assignment plus operator.
   //! @tparam T Value type of the right-hand side BraketVector.
   //! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
   //! @return BraketVector object, \f$ \boldsymbol{v} - \boldsymbol{v}_{\rm rhs}\f$.
   template<typename T>
   BraketVector& operator-=(const BraketVector<T> &rhs) {
      return *this = *this - rhs;
   }
   
   //! @brief Operator overloading: casting operator to std::vector<T>.
   //! @tparam T Value type of std::vector.
   template<typename T>
   operator std::vector<T>() const noexcept {
      std::vector<T> out(this->value_list.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < out.size(); ++i) {
         out[i] = static_cast<T>(this->value_list[i]);
      }
      return out;
   }
   
   //! @brief Operator overloading: assignment operator.
   //! @param vector The BraketVector to be assigned.
   //! @return BraketVector object.
   BraketVector &operator=(const BraketVector &vector) & {
      Assign(vector);
      return *this;
   }
   
private:
   //------------------------------------------------------------------
   //----------------------Private Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Calculate L1 norm.
   //! @tparam ReturnType Type of the return value.
   //! @return L1 norm.
   template<typename ReturnType>
   ReturnType CalculateL1Norm() const {
      ReturnType norm = 0.0;
#pragma omp parallel for reduction (+:norm)
      for (std::size_t i = 0; i < this->value_list.size(); ++i) {
         norm += std::abs(this->value_list[i]);
      }
      return norm;
   }
   
   //! @brief Calculate L2 norm.
   //! @tparam ReturnType Type of the return value.
   //! @return L2 norm.
   template<typename ReturnType>
   ReturnType CalculateL2Norm() const {
      ReturnType inner_product = 0.0;
#pragma omp parallel for reduction (+:inner_product)
      for (std::size_t i = 0; i < this->value_list.size(); ++i) {
         inner_product += this->value_list[i]*this->value_list[i];
      }
      return std::sqrt(inner_product);
   }
   
};

//------------------------------------------------------------------
//----------------------Operator overloading------------------------
//------------------------------------------------------------------
//! @brief Operator overloading: addition operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return BraketVector object, \f$ \boldsymbol{v}_{\rm lhs} + \boldsymbol{v}_{\rm rhs}\f$.
template<typename T1, typename T2>
auto operator+(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) ->
BraketVector<decltype(std::declval<T1>() + std::declval<T2>())> {
   return CalculateVectorVectorSum(T1{1}, lhs, T2{1}, rhs);
}

//! @brief Operator overloading: subtraction operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return BraketVector object, \f$ \boldsymbol{v}_{\rm lhs} - \boldsymbol{v}_{\rm rhs}\f$.
template<typename T1, typename T2>
auto operator-(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) ->
BraketVector<decltype(std::declval<T1>() - std::declval<T2>())> {
   return CalculateVectorVectorSum(T1{1}, lhs, T2{-1}, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return Inner product, \f$ \boldsymbol{v}_{\rm lhs}\cdot\boldsymbol{v}_{\rm rhs}\f$.
template<typename T1, typename T2>
auto operator*(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) ->
decltype(std::declval<T1>()*std::declval<T2>()) {
   return CalculateVectorVectorProduct(lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The value of the left-hand side, \f$ c_{\rm lhs}\f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return BraketVector object, \f$ c_{\rm lhs}\boldsymbol{v}_{\rm rhs}\f$.
template<typename T1, typename T2>
auto operator*(const T1 lhs, const BraketVector<T2> &rhs) ->
BraketVector<decltype(std::declval<T1>()*std::declval<T2>())> {
   return CalculateScalarVectorProduct(lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The value of the right-hand side, \f$ c_{\rm rhs}\f$.
//! @return BraketVector object, \f$ \boldsymbol{v}_{\rm lhs}c_{\rm rhs} = c_{\rm rhs}\boldsymbol{v}_{\rm lhs}\f$.
template<typename T1, typename T2>
auto operator*(const BraketVector<T1> &lhs, const T2 rhs) ->
BraketVector<decltype(std::declval<T1>()*std::declval<T2>())> {
   return CalculateScalarVectorProduct(rhs, lhs);
}

//! @brief Operator overloading: equality operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return Return true if \f$ \boldsymbol{v}_{\rm lhs}=\boldsymbol{v}_{\rm rhs}\f$, otherwise false.
template<typename T1, typename T2>
bool operator==(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   if (lhs.Size() != rhs.Size()) {
      return false;
   }
   const std::int64_t size = lhs.Size();
   for (std::int64_t i = 0; i < size; ++i) {
      if (lhs.value_list[i] != rhs.value_list[i]) {
         return false;
      }
   }
   return true;
}

//! @brief Operator overloading: equality operator.
//! @tparam T1 Value type of the left-hand side BraketVector.
//! @tparam T2 Value type of the right-hand side BraketVector.
//! @param lhs The BraketVector of the left-hand side, \f$ \boldsymbol{v}_{\rm lhs} \f$.
//! @param rhs The BraketVector of the right-hand side, \f$ \boldsymbol{v}_{\rm rhs} \f$.
//! @return Return true if \f$ \boldsymbol{v}_{\rm lhs}\neq\boldsymbol{v}_{\rm rhs}\f$, otherwise false.
template<typename T1, typename T2>
bool operator!=(const BraketVector<T1> &lhs, const BraketVector<T2> &rhs) {
   return !(lhs == rhs);
}

//! @brief Operator overloading: output operator.
//! @tparam T Type of BraketVector.
//! @param os Ostream object.
//! @param v BraketVector.
template<typename T>
std::ostream& operator<<(std::ostream &os, const BraketVector<T> &v) {
   os << std::fixed;
   os << std::setprecision(std::numeric_limits<T>::max_digits10);
   for (std::size_t i = 0; i < v.value_list.size(); ++i) {
      os << "v[";
      os << std::noshowpos << std::left << std::setw(3) << i << "]=";
      os << std::showpos << v.value_list[i] << std::endl;
   }
   return os;
}

//------------------------------------------------------------------
//------------------------Global Functions--------------------------
//------------------------------------------------------------------
//! @brief Calculate BraketVector summation, \f$ c_1\boldsymbol{v}_1 + c_2\boldsymbol{v}_2 \f$.
//! @tparam T1 Value type of coeff_1.
//! @tparam T2 Value type of braket_vector_1.
//! @tparam T3 Value type of coeff_2.
//! @tparam T4 Value type of braket_vector_2.
//! @param coeff_1 The coefficient \f$ c_1 \f$.
//! @param braket_vector_1 BraketVector \f$ \boldsymbol{v}_1\f$.
//! @param coeff_2 The coefficient \f$ c_2 \f$.
//! @param braket_vector_2 BraketVector \f$ \boldsymbol{v}_2\f$.
//! @return BraketVector object, \f$ c_1\boldsymbol{v}_{1} + c_2\boldsymbol{v}_{2}\f$.
template<typename T1, typename T2, typename T3, typename T4>
auto CalculateVectorVectorSum(const T1 coeff_1,
                              const BraketVector<T2> &braket_vector_1,
                              const T3 coeff_2,
                              const BraketVector<T4> &braket_vector_2) ->
BraketVector<decltype(std::declval<T1>()*std::declval<T2>() + std::declval<T3>()*std::declval<T4>())> {
   if (braket_vector_1.value_list.size() != braket_vector_2.value_list.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.value_list.size() << ", dim_2 = " << braket_vector_2.value_list.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   using T1T2T3T4 = decltype(std::declval<T1>()*std::declval<T2>() + std::declval<T3>()*std::declval<T4>());
   BraketVector<T1T2T3T4> vector_out;
   
   vector_out.Resize(braket_vector_1.value_list.size());
   const std::int64_t size = braket_vector_1.value_list.size();
#pragma omp parallel for
   for (std::int64_t i = 0; i < size; ++i) {
      vector_out.value_list[i] = coeff_1*braket_vector_1.value_list[i] + coeff_2*braket_vector_2.value_list[i];
   }
   return vector_out;
}

//! @brief Calculate BraketVector innner product, \f$ \boldsymbol{v}_1\cdot\boldsymbol{v}_2 \f$.
//! @tparam T1 Value type of braket_vector_1.
//! @tparam T2 Value type of braket_vector_2.
//! @param braket_vector_1 BraketVector \f$ \boldsymbol{v}_1\f$.
//! @param braket_vector_2 BraketVector \f$ \boldsymbol{v}_2\f$.
//! @return Inner product \f$ \boldsymbol{v}_1\cdot\boldsymbol{v}_2 \f$.
template<typename T1, typename T2>
auto CalculateVectorVectorProduct(const BraketVector<T1> &braket_vector_1,
                                  const BraketVector<T2> &braket_vector_2) ->
decltype(std::declval<T1>()*std::declval<T2>()) {
   if (braket_vector_1.value_list.size() != braket_vector_2.value_list.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.value_list.size() << ", dim_2 = " << braket_vector_2.value_list.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   using T1T2 = decltype(std::declval<T1>()*std::declval<T2>());
   T1T2 val_out = 0.0;
   const std::int64_t size = braket_vector_1.value_list.size();
#pragma omp parallel for reduction (+: val_out)
   for (std::int64_t i = 0; i < size; ++i) {
      val_out += braket_vector_1.value_list[i]*braket_vector_2.value_list[i];
   }
   return val_out;
}

//! @brief Calculate scalar BraketVector product, \f$ c\boldsymbol{v} \f$.
//! @tparam T1 Value type of value.
//! @tparam T2 Value type of braket_vector_1.
//! @param value The coefficient \f$ c\f$
//! @param braket_vector BraketVector \f$ \boldsymbol{v}\f$.
//! @return BraketVector object \f$ c\boldsymbol{v} \f$.
template<typename T1, typename T2>
auto CalculateScalarVectorProduct(const T1 value,
                                  const BraketVector<T2> &braket_vector) ->
BraketVector<decltype(std::declval<T1>()*std::declval<T2>())> {
   
   using T1T2 = decltype(std::declval<T1>()*std::declval<T2>());
   BraketVector<T1T2> out;
   out.Resize(braket_vector.value_list.size());
   const std::int64_t size = braket_vector.value_list.size();
#pragma omp parallel for
   for (std::int64_t i = 0; i < size; ++i) {
      out.value_list[i] = value*braket_vector.value_list[i];
   }
   return out;
}

//! @brief Calculate L1Norm between two BraketVector, aka Manhattan distance,
//! \f$ \|c_1\boldsymbol{v}_1 - c_2\boldsymbol{v}_2\|_1=\sum_{i}\|c_1v_{1,i} - c_2v_{2,i}\| \f$.
//! @tparam T1 Value type of coeff_1.
//! @tparam T2 Value type of braket_vector_1.
//! @tparam T3 Value type of coeff_2.
//! @tparam T4 Value type of braket_vector_2.
//! @param coeff_1 The coefficient \f$ c_1 \f$.
//! @param braket_vector_1 BraketVector \f$ \boldsymbol{v}_1\f$.
//! @param coeff_2 The coefficient \f$ c_2 \f$.
//! @param braket_vector_2 BraketVector \f$ \boldsymbol{v}_2\f$.
//! @return L1Norm \f$ \|c_1\boldsymbol{v}_1 - c_2\boldsymbol{v}_2\|_1 \f$.
template<typename T1, typename T2, typename T3, typename T4>
auto CalculateL1Norm(const T1 coeff_1,
                     const BraketVector<T2> &braket_vector_1,
                     const T3 coeff_2,
                     const BraketVector<T4> &braket_vector_2
                     ) ->
decltype(std::declval<T1>()*std::declval<T2>() - std::declval<T3>()*std::declval<T4>()) {
   
   if (braket_vector_1.value_list.size() != braket_vector_2.value_list.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.value_list.size() << ", dim_2 = " << braket_vector_2.value_list.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   using T1T2T3T4 = decltype(std::declval<T1>()*std::declval<T2>() - std::declval<T3>()*std::declval<T4>());
   T1T2T3T4 val_out = 0.0;

   const std::int64_t size = braket_vector_1.value_list.size();
#pragma omp parallel for reduction (+: val_out)
   for (std::int64_t i = 0; i < size; ++i) {
      val_out += std::abs(coeff_1*braket_vector_1.value_list[i] - coeff_2*braket_vector_2.value_list[i]);
   }
   return val_out;
}

//! @brief Calculate L2Norm, aka Euclidean distance, \f$ \|c_1\boldsymbol{v}_1 - c_2\boldsymbol{v}_2\| \f$.
//! @tparam T1 Value type of coeff_1.
//! @tparam T2 Value type of braket_vector_1.
//! @tparam T3 Value type of coeff_2.
//! @tparam T4 Value type of braket_vector_2.
//! @param coeff_1 The coefficient \f$ c_1 \f$.
//! @param braket_vector_1 BraketVector \f$ \boldsymbol{v}_1\f$.
//! @param coeff_2 The coefficient \f$ c_2 \f$.
//! @param braket_vector_2 BraketVector \f$ \boldsymbol{v}_2\f$.
//! @return L2Norm \f$ \|c_1\boldsymbol{v}_1 - c_2\boldsymbol{v}_2\| \f$.
template<typename T1, typename T2, typename T3, typename T4>
auto CalculateL2Norm(const T1 coeff_1,
                     const BraketVector<T2> &braket_vector_1,
                     const T3 coeff_2,
                     const BraketVector<T4> &braket_vector_2
                     ) ->
decltype(std::declval<T1>()*std::declval<T2>() - std::declval<T3>()*std::declval<T4>()) {
   
   if (braket_vector_1.value_list.size() != braket_vector_2.value_list.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "BraketVector types do not match each other" << std::endl;
      ss << "dim_1 = " << braket_vector_1.value_list.size() << ", dim_2 = " << braket_vector_2.value_list.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   using T1T2T3T4 = decltype(std::declval<T1>()*std::declval<T2>() - std::declval<T3>()*std::declval<T4>());
   T1T2T3T4 val_out = 0.0;
   
   const std::int64_t size = braket_vector_1.value_list.size();
#pragma omp parallel for reduction (+: val_out)
   for (std::int64_t i = 0; i < size; ++i) {
      const T1T2T3T4 v = coeff_1*braket_vector_1.value_list[i] - coeff_2*braket_vector_2.value_list[i];
      val_out += v*v;
   }
   return std::sqrt(val_out);
}

} // namespace type
} // namespace compnal


#endif /* COMPNAL_TYPE_BRAKET_VECTOR_HPP_ */
