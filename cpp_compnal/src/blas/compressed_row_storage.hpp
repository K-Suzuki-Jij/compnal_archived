//
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
//  compressed_row_storage.hpp
//  compnal
//
//  Created by kohei on 2022/10/06.
//  
//

#ifndef COMPNAL_BLAS_COMPRESSED_ROW_STORAGE_HPP_
#define COMPNAL_BLAS_COMPRESSED_ROW_STORAGE_HPP_

#include "../utility/sort.hpp"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace compnal {
namespace blas {

//! @brief Enumerated type to represent Fermion, Boson, etc..
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


template<typename RealType>
struct CRS {
   
   //! @brief The dimension of the row.
   std::int64_t row_dim = 0;
   
   //! @brief The dimension of the colmun.
   std::int64_t col_dim = 0;
   
   //! @brief The locations in the val vector that start a row.
   std::vector<std::int64_t> row = {0};
   
   //! @brief The column indexes of the elements in the val vector.
   std::vector<std::int64_t> col;
   
   //! @brief The values of the nonzero elements of the matrix.
   std::vector<RealType> val;
   
   //! @brief Type of the the matrix.
   CRSTag tag = CRSTag::NONE;
   
   //! @brief Matrix name.
   std::string name = "";
   
   template <typename T>
   void MultiplyByScalar(const T coeff) {
      if (std::abs(coeff) < std::numeric_limits<RealType>::epsilon() < T{0}) {
#pragma omp parallel for
         for (std::int64_t i = 0; i < this->row_dim; ++i) {
            this->row[i + 1] = 0;
         }
         this->col.clear();
         this->val.clear();
      } else {
#pragma omp parallel for
         for (std::size_t i = 0; i < this->col.size(); ++i) {
            this->val[i] *= coeff;
         }
      }
   }
   
   //! @brief Add a value to the diagonal elements.
   //! @tparam T Value type.
   //! @param diag_add The value to be added to the diagonal elements.
   template <typename T>
   void AddDiagonalElements(const T diag_add) {
      if (this->row_dim != this->col_dim) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "The matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }

      bool flag = false;
#pragma omp parallel for
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         bool temp_flag = true;
         for (std::int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
            if (i == this->col[j]) {
               this->val[j] += diag_add;
               temp_flag = false;
               break;
            }
         }
         if (temp_flag) {
#pragma omp critical
            { flag = temp_flag; }
         }
      }
      if (flag) {
         // Restore the diagonal elements.
#pragma omp parallel for
         for (std::int64_t i = 0; i < this->row_dim; ++i) {
            bool temp_flag = true;
            for (std::int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
               if (i == this->col[j]) {
                  this->val[j] -= diag_add;
                  temp_flag = false;
                  break;
               }
            }
         }
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "Could not add the value to the diagonal elements." << std::endl;
         ss << "The matrix need to have all the diagonal elements even if "
               "they "
               "are zero."
            << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
   void SortCol() {
#pragma omp parallel for schedule(guided)
      for (std::int64_t i = 0; i < this->row_dim; i++) {
         utility::QuickSortVector(&this->col, &this->val, this->row[i], this->row[i + 1]);
      }
   }
   
   //! @brief Check if the matrix is symmetric or not within the threshold.
   //! @param threshold The threshold. Defaults to 10^-15.
   //! @param flag_display_info Display information if the matrix is not
   //! symmetric.
   //! @return Return true if the matrix is symmetric, otherwise false.
   bool CheckSymmetric(const RealType threshold = std::numeric_limits<RealType>::epsilon(),
                       const bool flag_display_info = true) const {
      if (this->row_dim != this->col_dim) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "The matrix is not square one." << std::endl;
         throw std::runtime_error(ss.str());
      }
      // Check sorted
      bool flag_sorted = false;
#pragma omp parallel for
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         bool temp_flag = false;
         for (std::int64_t j = this->row[i]; j < this->row[i + 1] - 1; ++j) {
            if (this->col[j] >= this->col[j + 1]) {
               temp_flag = true;
            }
         }
         if (temp_flag) {
#pragma omp critical
            { flag_sorted = temp_flag; }
         }
      }
      if (flag_sorted) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "The colmun indexes of the matrix is not sorted." << std::endl;
         throw std::runtime_error(ss.str());
      }
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         for (std::int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
            const auto iter_begin = this->col.begin() + this->row[this->col[j]];
            const auto iter_end = this->col.begin() + this->row[this->col[j] + 1];
            const auto iter_find = std::lower_bound(iter_begin, iter_end, i);
            if (iter_find == iter_end || *iter_find != i) {
               if (flag_display_info) {
                  std::cout << "The input matrix is not symmetric." << std::endl;
                  std::cout << "Corresponding element does not exist." << std::endl;
                  std::cout << "row=" << i << ", col=" << col[j] << ", val=" << val[j] << std::endl;
               }
               return false;
            }
            const auto inv = std::distance(this->col.begin(), iter_find);
            if (std::abs(this->val[j] - this->val[inv]) > threshold) {
               if (flag_display_info) {
                  std::cout << "The input matrix is not symmetric." << std::endl;
                  std::cout << "M[" << i << "][" << this->col[j] << "]=" << this->val[j] << ", " << this->val[inv]
                            << "=M[" << this->col[j] << "][" << i << "]" << std::endl;
               }
               return false;
            }
         }
      }
      return true;
   }
   
   //! @brief Print the matrix.
   void Print(const std::string display_name = "") const {
      if (this->name != "") {
         std::cout << "Name: " << this->name << std::endl;
      }
      std::cout << this->tag << std::endl;
      std::cout << std::fixed;
      std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10);
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         for (std::int64_t j = this->row.at(i); j < this->row.at(i + 1); ++j) {
            std::cout << display_name << "[";
            std::cout << std::noshowpos << std::left << std::setw(3) << i << "][";
            std::cout << std::left << std::setw(3) << this->col[j] << "]=";
            std::cout << std::showpos << this->val[j] << std::endl;
         }
      }
      std::cout << std::noshowpos;
   }
   
   
   //! @brief Print the information about the matrix.
   void PrintInfo(const std::string display_name = "Matrix") const {
      std::cout << std::fixed;
      std::cout << std::setprecision(std::numeric_limits<RealType>::max_digits10);
      std::cout << "Print information about CRS: " << display_name << std::endl;
      std::cout << this->tag << std::endl;
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
   
};



} // namespace blas
} // namespace compnal

#endif /* COMPNAL_BLAS_COMPRESSED_ROW_STORAGE_HPP_ */
