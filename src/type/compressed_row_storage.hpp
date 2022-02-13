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
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef COMPNAL_TYPE_COMPRESSED_ROW_STORAGE_HPP_
#define COMPNAL_TYPE_COMPRESSED_ROW_STORAGE_HPP_

#include "../utility/sort.hpp"

#include <iostream>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <sstream>

namespace compnal {
namespace type {

//! @brief Enumerated type to represent Fermion or not.
enum CRSTag {
   
   //! @brief No taged.
   NONE = 0,
   
   //! @brief Fermionic operator.
   FERMION = 1,
   
   //! @brief Bosonic operator.
   BOSON = 2,
   
   //! @brief Sum of Fermion and Boson.
   MIX = 3
};

//! @brief Sparse matrix classs (Compressed Row Strage: CRS).
//! @tparam ValueType Value type.
template<typename ValueType>
struct CRS {
   //------------------------------------------------------------------
   //----------------------Public Member Variables---------------------
   //------------------------------------------------------------------
   std::int64_t row_dim = 0;
   std::int64_t col_dim = 0;
   std::vector<std::int64_t> row;
   std::vector<std::int64_t> col;
   std::vector<ValueType> val;
   CRSTag tag = CRSTag::NONE;

   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   CRS(const std::int64_t row_dim_in = 0, const std::int64_t col_dim_in = 0, const CRSTag tag_in = CRSTag::NONE) {
      
      this->row_dim = row_dim_in;
      this->col_dim = col_dim_in;
      this->row.resize(row_dim_in + 1);
#pragma omp parallel for
      for (std::int64_t i = 0; i <= row_dim_in; ++i) {
         this->row[i] = 0;
      }
      this->tag = tag_in;
   }
   
   explicit CRS(const std::vector<std::vector<ValueType>> &mat_vec, const CRSTag tag_in = CRSTag::NONE) {
      this->row_dim = mat_vec.size();
      this->col_dim = 0;
      this->row.resize(this->row_dim + 1);
      this->row[0] = 0;
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         const std::int64_t size = static_cast<std::int64_t>(mat_vec[i].size());
         for (std::int64_t j = 0; j < size;++j) {
            if (mat_vec[i][j] != 0.0) {
               this->col.push_back(j);
               this->val.push_back(mat_vec[i][j]);
            }
            if (this->col_dim < j + 1) {
               this->col_dim = j + 1;
            }
         }
         this->row[i + 1] = this->col.size();
      }
      this->tag = tag_in;
   }
   
   CRS(const CRS &matrix) {
      Assign(matrix);
   }
   
   CRS &operator=(const CRS &matrix) & {
      Assign(matrix);
      return *this;
   }
      
   //------------------------------------------------------------------
   //----------------------Public Member Functions---------------------
   //------------------------------------------------------------------
   void Assign(const CRS &matrix) {
      this->row_dim = matrix.row_dim;
      this->col_dim = matrix.col_dim;
      this->row.resize(matrix.row_dim + 1);
      
      this->col.resize(matrix.col.size());
      this->val.resize(matrix.col.size());
      
#pragma omp parallel for
      for (std::int64_t i = 0; i <= matrix.row_dim; ++i) {
         this->row[i] = matrix.row[i];
      }
      
#pragma omp parallel for
      for (std::size_t i = 0; i < matrix.col.size(); ++i) {
         this->col[i] = matrix.col[i];
         this->val[i] = matrix.val[i];
      }
      this->tag = matrix.tag;
   }
   
   void MultiplyByScalar(const ValueType coeef) {
      if (coeef == 0.0) {
#pragma omp parallel for
         for (std::int64_t i = 0; i < this->row_dim; ++i) {
            this->row[i + 1] = 0;
         }
         this->col.clear();
         this->val.clear();
      }
      else {
#pragma omp parallel for
         for (std::size_t i = 0; i < this->col.size(); ++i) {
            this->val[i] *= coeef;
         }
      }
   }
      
   void DiagonalScaling(const ValueType diag_add) {
      if (this->row_dim != this->col_dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
#pragma omp parallel for
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         bool flag = true;
         for (std::int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
            if (i == this->col[j]) {
               this->val[j] += diag_add;
               flag = false;
               break;
            }
         }
         if (flag) {
            std::stringstream ss;
            ss << "Error in " << __func__ << std::endl;
            ss << "Some of the diagonal components are not registered." << std::endl;
            throw std::runtime_error(ss.str());
         }
      }
   }
   
   void Free() {
      this->row_dim = 0;
      this->col_dim = 0;
      std::vector<std::int64_t>().swap(this->row);
      std::vector<std::int64_t>().swap(this->col);
      std::vector<ValueType>().swap(this->val);
      this->row.push_back(0);
      this->tag = CRSTag::NONE;
   }
   
   void Clear() {
      this->row_dim = 0;
      this->col_dim = 0;
      this->row.clear();
      this->col.clear();
      this->val.clear();
      this->row.push_back(0);
      this->tag = CRSTag::NONE;
   }
      
   void SortCol() {
#pragma omp parallel for schedule (guided)
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         utility::QuickSort<std::int64_t, ValueType>(&this->col, &this->val, this->row[i], this->row[i + 1]);
      }
   }
   
   void Print(const std::string display_name = "Matrix") const {
      std::cout << std::fixed;
      std::cout << std::setprecision(std::numeric_limits<ValueType>::max_digits10);
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         for (std::int64_t j = this->row.at(i); j < this->row.at(i+1); ++j) {
            std::cout << display_name << "[";
            std::cout << std::noshowpos << std::left << std::setw(3) << i << "][";
            std::cout << std::left << std::setw(3) << this->col[j] << "]=";
            std::cout << std::showpos << this->val[j] << std::endl;
         }
      }
      std::cout << std::noshowpos;
   }
   
   void PrintInfo(const std::string display_name = "Matrix") const {
      std::cout << std::fixed;
      std::cout << std::setprecision(std::numeric_limits<ValueType>::max_digits10);
      std::cout << "Print information about CRS: " << display_name << std::endl;
      std::cout << "row_dim = " << this->row_dim << std::endl;
      std::cout << "col_dim = " << this->col_dim << std::endl;
      for (std::size_t i = 0; i < this->row.size(); ++i) {
         std::cout << "row[" << i << "] = " << this->row.at(i) << std::endl;
      }
      for (std::size_t i = 0; i < this->col.size(); ++i) {
         std::cout << "col[" << i << "] = " << this->col.at(i) << std::endl;
      }
      for (std::size_t i = 0; i < this->val.size(); ++i) {
         std::cout << "val[" << i << "] = " << this->val.at(i) << std::endl;
      }
   }
   
   bool isSymmetric(const ValueType threshold = 0.000000000000001/*pow(10,-15)*/) const {
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         for (std::int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
            const auto iter_begin = this->col.begin() + this->row[col[j]];
            const auto iter_end   = this->col.begin() + this->row[col[j] + 1];
            const auto iter_find  = std::lower_bound(iter_begin, iter_end, i);
            if (iter_find == iter_end || *iter_find != i) {
               std::cout << "The input matrix is not symmetric." << std::endl;
               std::cout << "Corresponding element does not exist." << std::endl;
               std::cout << "row=" << i << ", col=" << col[j] << ", val=" << val[j] << std::endl;
               return false;
            }
            const auto inv = std::distance(iter_begin, iter_find);
            if (std::abs(this->val[j] - this->val[inv]) > threshold) {
               std::cout << "The input matrix is not symmetric." << std::endl;
               std::cout << "M[" << i << "][" << this->col[j] << "]=" << this->val[j] << ", " << this->val[inv] << "=M[" << this->col[j] << "][" << i << "]" << std::endl;
               return false;
            }
         }
      }
      return true;
   }
   
   //------------------------------------------------------------------
   //-----------------------Operator Overloading-----------------------
   //------------------------------------------------------------------
   //! @brief Operator overloading: unary plus operator.
   CRS operator+() const {
      return *this;
   }
   
   //! @brief Operator overloading: unary negation operator.
   CRS operator-() const {
      MultiplyByScalar(ValueType{-1.0});
      return *this;
   }
   
   //! @brief Operator overloading: compound assignment plus operator.
   //! @tparam T Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template<typename T>
   CRS& operator+=(const CRS<T> rhs) {
      return *this = *this + rhs;
   }
   
   //! @brief Operator overloading: compound assignment subtraction operator.
   //! @tparam T Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template<typename T>
   CRS& operator-=(const CRS<T> rhs) {
      return *this = *this - rhs;
   }
   
   //! @brief Operator overloading: compound assignment multiplication operator.
   //! @tparam T Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template<typename T>
   CRS& operator*=(const CRS<T> rhs) {
      return *this = *this * rhs;
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
CRS<decltype(T1{0}+T2{0})> operator+(const CRS<T1> &lhs, const CRS<T2> &rhs) {
   return CalculateMatrixMatrixSum(T1{1}, lhs, T2{1}, rhs);
}

//! @brief Operator overloading: subtraction operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
CRS<decltype(T1{0}-T2{0})> operator-(const CRS<T1> &lhs, const CRS<T2> &rhs) {
   return CalculateMatrixMatrixSum(T1{1}, lhs, T2{-1}, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
CRS<decltype(T1{0}*T2{0})> operator*(const CRS<T1> &lhs, const CRS<T2> &rhs) {
   return CalculateMatrixMatrixProduct(T1{1}, lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
CRS<decltype(T1{0}*T2{0})> operator*(const T1 lhs, const CRS<T2> &rhs) {
   return CalculateScalarMatrixProduct(lhs, rhs);
}

//! @brief Operator overloading: multiplication operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
CRS<decltype(T1{0}*T2{0})> operator*(const CRS<T1> &lhs, const T2 rhs) {
   return CalculateScalarMatrixProduct(rhs, lhs);
}

//! @brief Operator overloading: equality operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
bool operator==(const CRS<T1> &lhs, const CRS<T2> &rhs) {
   if (lhs.row_dim != rhs.row_dim) {
      return false;
   }
   if (lhs.col_dim != rhs.col_dim) {
      return false;
   }
   if (lhs.row.size() != rhs.row.size()) {
      return false;
   }
   if (lhs.tag != rhs.tag) {
      return false;
   }
   if (lhs.col.size() != rhs.col.size()) {
      return false;
   }
   if (lhs.val.size() != rhs.val.size()) {
      return false;
   }
   for (std::size_t i = 0; i < lhs.row.size(); ++i) {
      if (lhs.row[i] != rhs.row[i]) {
         return false;
      }
   }
   for (std::size_t i = 0; i < lhs.col.size(); ++i) {
      if (lhs.col[i] != rhs.col[i]) {
         return false;
      }
   }
   for (std::size_t i = 0; i < lhs.val.size(); ++i) {
      if (lhs.val[i] != rhs.val[i]) {
         return false;
      }
   }
   return true;
}

//! @brief Operator overloading: inequality operator.
//! @tparam T1 Value type of the left-hand side.
//! @tparam T2 Value type of the right-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename T1, typename T2>
bool operator!=(const CRS<T1> &lhs, const CRS<T2> &rhs) {
   return !(lhs == rhs);
}

//------------------------------------------------------------------
//----------------Operator overloading: I/O Stream------------------
//------------------------------------------------------------------
//! @brief Operator overloading: output operator.
//! @tparam ValueType Value type CRS matrix.
//! @param os Ostream object.
//! @param m The matrix as CRS form.
template<typename ValueType>
std::ostream& operator<<(std::ostream &os, const CRS<ValueType> &m) {
   os << std::fixed;
   os << std::setprecision(std::numeric_limits<ValueType>::max_digits10);
   for (std::int64_t i = 0; i < m.row_dim; ++i) {
      for (std::int64_t j = m.row.at(i); j < m.row.at(i+1); ++j) {
         os << "M[";
         os << std::noshowpos << std::left << std::setw(3) << i << "][";
         os << std::left << std::setw(3) << m.col[j] << "]=";
         os << std::showpos << m.val[j] << std::endl;
      }
   }
   return os;
}

template<typename T1, typename T2, typename T3>
CRS<decltype(T1{0}*T2{0}*T3{0})> CalculateMatrixMatrixProduct(const T1 coeef_1,
                                                              const CRS<T2> &matrix_1,
                                                              const CRS<T3> &matrix_2) {
   
   if (matrix_1.col_dim != matrix_2.row_dim) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "Matrix product cannot be defined" << std::endl;
      ss << "matrix_1.col_dim = " << matrix_1.col_dim << ", matrix_2.row_dim = " << matrix_2.row_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   CRS<decltype(T1{0}*T2{0}*T3{0})> matrix_out(matrix_1.row_dim, matrix_2.col_dim);
   
   std::vector<decltype(T1{0}*T2{0})>       temp_v1(matrix_1.col_dim, 0.0);
   std::vector<decltype(T1{0}*T2{0}*T3{0})> temp_v2(matrix_2.col_dim, 0.0);
   
   for (std::int64_t i = 0; i < matrix_1.row_dim; ++i) {
      for (std::int64_t j = matrix_1.row[i]; j < matrix_1.row[i + 1]; ++j) {
         temp_v1[matrix_1.col[j]] = coeef_1*matrix_1.val[j];
      }
      for (std::int64_t j = 0; j < matrix_1.col_dim; ++j) {
         for (std::int64_t k = matrix_2.row[j]; k < matrix_2.row[j + 1]; ++k) {
            temp_v2[matrix_2.col[k]] += temp_v1[j]*matrix_2.val[k];
         }
      }
      for (std::int64_t j = 0; j < matrix_2.col_dim; ++j) {
         if (std::abs(temp_v2[j]) > 0.0) {
            matrix_out.val.push_back(temp_v2[j]);
            matrix_out.col.push_back(j);
         }
      }
      
      matrix_out.row[i + 1] = matrix_out.col.size();
      
      for (std::int64_t j = matrix_1.row[i]; j < matrix_1.row[i + 1]; ++j) {
         temp_v1[matrix_1.col[j]] = 0.0;
      }
      for (std::int64_t j = matrix_out.row[i]; j < matrix_out.row[i + 1]; ++j) {
         temp_v2[matrix_out.col[j]] = 0.0;
      }
   }
   
   if (matrix_1.tag == CRSTag::NONE) {
      matrix_out.tag = matrix_2.tag;
   }
   else if (matrix_1.tag == CRSTag::FERMION) {
      if (matrix_2.tag == CRSTag::NONE || matrix_2.tag == CRSTag::BOSON) {
         matrix_out.tag = CRSTag::FERMION;
      }
      else if (matrix_2.tag == CRSTag::FERMION) {
         matrix_out.tag = CRSTag::BOSON;
      }
      else if (matrix_2.tag == CRSTag::MIX) {
         matrix_out.tag = CRSTag::MIX;
      }
      else {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Unknown CRSTag detected.";
         throw std::runtime_error(ss.str());
      }
   }
   else if (matrix_1.tag == CRSTag::BOSON) {
      if (matrix_2.tag == CRSTag::FERMION || matrix_2.tag == CRSTag::BOSON || matrix_2.tag == CRSTag::MIX) {
         matrix_out.tag = matrix_2.tag;
      }
      else if (matrix_2.tag == CRSTag::NONE) {
         matrix_out.tag = CRSTag::BOSON;
      }
      else {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Unknown CRSTag detected.";
         throw std::runtime_error(ss.str());
      }
   }
   else if (matrix_1.tag == CRSTag::MIX) {
      matrix_out.tag = CRSTag::MIX;
   }
   else {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "Unknown CRSTag detected.";
      throw std::runtime_error(ss.str());
   }
   
   return matrix_out;
}

template<typename T1, typename T2, typename T3, typename T4>
CRS<decltype(T1{0}+T2{0}+T3{0}+T4{0})> CalculateMatrixMatrixSum(const T1 coeef_1,
                                                                const CRS<T2> &matrix_1,
                                                                const T3 coeef_2,
                                                                const CRS<T4> &matrix_2) {
   
   if (matrix_1.row_dim != matrix_2.row_dim || matrix_1.col_dim != matrix_2.col_dim) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "The summation of the matrices cannot be defined." << std::endl;
      ss << "matrix_1.row_dim = " << matrix_1.row_dim << ", matrix_1.col_dim = " << matrix_1.col_dim << std::endl;
      ss << "matrix_2.row_dim = " << matrix_2.row_dim << ", matrix_2.col_dim = " << matrix_2.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
      
   CRS<decltype(T1{0}+T2{0}+T3{0}+T4{0})> matrix_out(matrix_1.row_dim, matrix_1.col_dim);
   
   for (std::int64_t i = 0; i < matrix_1.row_dim; ++i) {
      
      int check = 0;
      std::int64_t count_1 = 0;
      std::int64_t count_2 = 0;
      
      const std::int64_t row_lower_1 = matrix_1.row[  i  ];
      const std::int64_t row_upper_1 = matrix_1.row[i + 1];
      const std::int64_t row_lower_2 = matrix_2.row[  i  ];
      const std::int64_t row_upper_2 = matrix_2.row[i + 1];
      
      const std::int64_t m1_count = row_upper_1 - row_lower_1;
      const std::int64_t m2_count = row_upper_2 - row_lower_2;
      
      if (m1_count != 0 && m2_count == 0) {
         for (std::int64_t j = row_lower_1; j < row_upper_1; ++j) {
            matrix_out.val.push_back(coeef_1*matrix_1.val[j]);
            matrix_out.col.push_back(matrix_1.col[j]);
         }
      }
      else if (m1_count == 0 && m2_count != 0) {
         for (std::int64_t j = row_lower_2; j < row_upper_2; ++j) {
            matrix_out.val.push_back(coeef_2*matrix_2.val[j]);
            matrix_out.col.push_back(matrix_2.col[j]);
         }
      }
      else if (m1_count != 0 && m2_count != 0) {
         for (std::int64_t j = 0; j < m1_count + m2_count; ++j) {
            if (matrix_1.col[row_lower_1 + count_1] < matrix_2.col[row_lower_2 + count_2]) {
               matrix_out.val.push_back(coeef_1*matrix_1.val[row_lower_1 + count_1]);
               matrix_out.col.push_back(matrix_1.col[row_lower_1 + count_1]);
               count_1++;
               if (row_lower_1 + count_1 == row_upper_1) {
                  check = 1;
                  break;
               }
            }
            else if (matrix_1.col[row_lower_1 + count_1] == matrix_2.col[row_lower_2 + count_2]) {
               const auto val = coeef_1*matrix_1.val[row_lower_1 + count_1] + coeef_2*matrix_2.val[row_lower_2 + count_2];
               if (std::abs(val) > 0.0) {
                  matrix_out.val.push_back(val);
                  matrix_out.col.push_back(matrix_1.col[row_lower_1 + count_1]);
                  count_1++;
                  count_2++;
                  const std::int64_t temp_count_1 = row_lower_1 + count_1;
                  const std::int64_t temp_count_2 = row_lower_2 + count_2;
                  if (temp_count_1 == row_upper_1 && temp_count_2 < row_upper_2) {
                     check = 1;
                     break;
                  }
                  else if (temp_count_1 < row_upper_1 && temp_count_2 == row_upper_2) {
                     check = 2;
                     break;
                  }
                  else if (temp_count_1 == row_upper_1 && temp_count_2 == row_upper_2) {
                     check = 0;
                     break;
                  }
               }
               else {
                  count_1++;
                  count_2++;
                  const std::int64_t temp_count_1 = row_lower_1 + count_1;
                  const std::int64_t temp_count_2 = row_lower_2 + count_2;
                  if (temp_count_1 == row_upper_1 && temp_count_2 < row_upper_2) {
                     check = 1;
                     break;
                  }
                  else if (temp_count_1 < row_upper_1 && temp_count_2 == row_upper_2) {
                     check = 2;
                     break;
                  }
                  else if (temp_count_1 == row_upper_1 && temp_count_2 == row_upper_2) {
                     check = 0;
                     break;
                  }
               }
            }
            else {
               matrix_out.val.push_back(coeef_2*matrix_2.val[row_lower_2 + count_2]);
               matrix_out.col.push_back(matrix_2.col[row_lower_2 + count_2]);
               count_2++;
               if (row_lower_2 + count_2 == row_upper_2) {
                  check = 2;
                  break;
               }
            }
         }
         if (check == 1) {
            for (std::int64_t j = row_lower_2 + count_2; j < row_upper_2; ++j) {
               matrix_out.val.push_back(coeef_2*matrix_2.val[j]);
               matrix_out.col.push_back(matrix_2.col[j]);
            }
         }
         else if (check == 2) {
            for (std::int64_t j = row_lower_1 + count_1; j < row_upper_1; ++j) {
               matrix_out.val.push_back(coeef_1*matrix_1.val[j]);
               matrix_out.col.push_back(matrix_1.col[j]);
            }
         }
      }
      matrix_out.row[i + 1] = matrix_out.col.size();
   }
   
   matrix_out.tag = matrix_1.tag;
   return matrix_out;
}

template<typename T1, typename T2>
CRS<decltype(T1{0}*T2{0})> CalculateScalarMatrixProduct(const T1 lhs, const CRS<T2> &rhs) {
   
   CRS<decltype(T1{0}*T2{0})> out(rhs.row_dim, rhs.col_dim, rhs.tag);
   out.col.resize(rhs.col.size());
   out.val.resize(rhs.val.size());
   
#pragma omp parallel for
   for (std::size_t i = 0; i < rhs.col.size(); ++i) {
      out.col[i] = rhs.col[i];
      out.val[i] = lhs*rhs.val[i];
   }
   
#pragma omp parallel for
   for (std::int64_t i = 1; i <= rhs.row_dim; ++i) {
      out.row[i] = rhs.row[i];
   }
   return out;
}

template<typename ValueType>
CRS<ValueType> CalculateTransposedMatrix(const CRS<ValueType> &matrix_in) {

   CRS<ValueType> matrix_out(matrix_in.col_dim, matrix_in.row_dim, matrix_in.tag);
   
   std::vector<std::int64_t> row_count(matrix_in.row_dim);
   for (std::int64_t i = 0; i < matrix_in.col_dim; ++i) {
      for (std::int64_t j = 0; j < matrix_in.row_dim; ++j) {
         const std::int64_t row = matrix_in.row[j] + row_count[j];
         if (row < matrix_in.row[j + 1] && matrix_in.col[row] == i) {
            matrix_out.val.push_back(matrix_in.val[row]);
            matrix_out.col.push_back(j);
            row_count[j]++;
         }
      }
      matrix_out.row[i + 1] = matrix_out.col.size();
   }
   
   return matrix_out;
}

} // namespace type
} // namespace compnal


#endif /* COMPNAL_TYPE_COMPRESSED_ROW_STORAGE_HPP_ */
