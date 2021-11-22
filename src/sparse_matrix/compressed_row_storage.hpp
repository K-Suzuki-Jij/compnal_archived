//
//  compressed_row_storage.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef COMPNAL_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_
#define COMPNAL_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_

#include "../utility/all.hpp"
#include <iostream>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <sstream>

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
struct CRS {
   
   std::int64_t row_dim = 0;
   std::int64_t col_dim = 0;
   std::vector<std::int64_t> row;
   std::vector<std::int64_t> col;
   std::vector<RealType>     val;
   
   CRS(const std::int64_t row_dim_in = 0, const std::int64_t col_dim_in = 0) {
      this->row_dim = row_dim_in;
      this->col_dim = col_dim_in;
      this->row.resize(row_dim_in + 1);
#pragma omp parallel for
      for (std::int64_t i = 0; i <= row_dim_in; ++i) {
         this->row[i] = 0;
      }
   }
   
   explicit CRS(const std::vector<std::vector<RealType>> &mat_vec) {
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
   }
   
   CRS(const CRS &matrix) {
      Assign(matrix);
   }
   
   CRS &operator=(const CRS &matrix) & {
      Assign(matrix);
      return *this;
   }
   
   void UpdateRow() {
      this->row.push_back(this->col.size());
   }
   
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
   }
   
   void MultiplyByScalar(const RealType coeef) {
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
   
   CRS MultiplyByScalar(const RealType coeef) const {
      if (coeef == 0.0) {
         return CRS(this->row_dim, this->col_dim);
      }
      else {
         CRS m = *this;
#pragma omp parallel for
         for (std::size_t i = 0; i < m.col.size(); ++i) {
            m.val[i] *= coeef;
         }
         return m;
      }
   }
   
   void DiagonalScaling(const RealType diag_add) {
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
#pragma omp critical
            {
               std::stringstream ss;
               ss << "Error in " << __func__ << std::endl;
               ss << "Some of the diagonal components are not registered." << std::endl;
               throw std::runtime_error(ss.str());
            }
         }
      }
   }
   
   void Free() {
      this->row_dim = 0;
      this->col_dim = 0;
      std::vector<std::int64_t>().swap(this->row);
      std::vector<std::int64_t>().swap(this->col);
      std::vector<RealType>().swap(this->val);
      this->row.push_back(0);
   }
   
   void Clear() {
      this->row_dim = 0;
      this->col_dim = 0;
      this->row.clear();
      this->col.clear();
      this->val.clear();
      this->row.push_back(0);
   }
   
   void SortCol() {
#pragma omp parallel for schedule (guided)
      for (std::int64_t i = 0; i < this->row_dim; ++i) {
         utility::QuickSort<std::int64_t, RealType>(&this->col, &this->val, this->row[i], this->row[i + 1]);
      }
   }
   
   void PrintMatrix(const std::string display_name = "Matrix") const {
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
   
   bool isSymmetric(const RealType threshold = 0.000000000000001/*pow(10,-15)*/) const {
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
   
};

template<typename RealType>
void CalculateTransposedMatrix(CRS<RealType> *matrix_out,
                               const CRS<RealType> &matrix_in) {
   
   std::vector<std::int64_t> row_count(matrix_in.row_dim);
   matrix_out->Clear();
   for (std::int64_t i = 0; i < matrix_in.col_dim; ++i) {
      for (std::int64_t j = 0; j < matrix_in.row_dim; ++j) {
         const std::int64_t row = matrix_in.row[j] + row_count[j];
         if (row < matrix_in.row[j + 1] && matrix_in.col[row] == i) {
            matrix_out->val.push_back(matrix_in.val[row]);
            matrix_out->col.push_back(j);
            row_count[j]++;
         }
      }
      matrix_out->row.push_back(matrix_out->col.size());
   }
   matrix_out->row_dim = matrix_in.col_dim;
   matrix_out->col_dim = matrix_in.row_dim;
}


template<typename RealType>
CRS<RealType> CalculateTransposedMatrix(const CRS<RealType> &matrix_in) {
   CRS<RealType> matrix_out;
   CalculateTransposedMatrix(&matrix_out, matrix_in);
   return matrix_out;
}

template<typename RealType>
CRS<RealType> CalculateMatrixMatrixProduct(const RealType coeef_1,
                                           const CRS<RealType> &matrix_1,
                                           const RealType coeef_2,
                                           const CRS<RealType> &matrix_2) {
   CRS<RealType> m;
   CalculateMatrixMatrixProduct(&m, coeef_1, matrix_1, coeef_2, matrix_2);
   return m;
}

template<typename RealType>
void CalculateMatrixMatrixProduct(CRS<RealType> *matrix_out,
                                  const RealType coeef_1,
                                  const CRS<RealType> &matrix_1,
                                  const RealType coeef_2,
                                  const CRS<RealType> &matrix_2) {
   
   if (matrix_1.col_dim != matrix_2.row_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Matrix product cannot be defined" << std::endl;
      ss << "matrix_1.col_dim = " << matrix_1.col_dim << ", matrix_2.row_dim = " << matrix_2.row_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   matrix_out->Clear();
   matrix_out->row.resize(matrix_1.row_dim + 1);
   
   std::vector<RealType> temp_v1(matrix_1.col_dim, 0.0);
   std::vector<RealType> temp_v2(matrix_2.col_dim, 0.0);
   
   for (std::int64_t i = 0; i < matrix_1.row_dim; ++i) {
      for (std::int64_t j = matrix_1.row[i]; j < matrix_1.row[i + 1]; ++j) {
         temp_v1[matrix_1.col[j]] = coeef_1*matrix_1.val[j];
      }
      for (std::int64_t j = 0; j < matrix_1.col_dim; ++j) {
         for (std::int64_t k = matrix_2.row[j]; k < matrix_2.row[j + 1]; ++k) {
            temp_v2[matrix_2.col[k]] += temp_v1[j]*coeef_2*matrix_2.val[k];
         }
      }
      for (std::int64_t j = 0; j < matrix_2.col_dim; ++j) {
         if (std::abs(temp_v2[j]) > 0.0) {
            matrix_out->val.push_back(temp_v2[j]);
            matrix_out->col.push_back(j);
         }
      }
      
      matrix_out->row[i + 1] = matrix_out->col.size();
      
      for (std::int64_t j = matrix_1.row[i]; j < matrix_1.row[i + 1]; ++j) {
         temp_v1[matrix_1.col[j]] = 0.0;
      }
      for (std::int64_t j = matrix_out->row[i]; j < matrix_out->row[i + 1]; ++j) {
         temp_v2[matrix_out->col[j]] = 0.0;
      }
   }
   
   matrix_out->row_dim = matrix_1.row_dim;
   matrix_out->col_dim = matrix_2.col_dim;
   
}

template<typename RealType>
CRS<RealType> CreateMatrixSum(const RealType coeef_1,
                              const CRS<RealType> &matrix_1,
                              const RealType coeef_2,
                              const CRS<RealType> &matrix_2) {
   
   CRS<RealType> m;
   CreateMatrixSum(&m, coeef_1, matrix_1, coeef_2, matrix_2);
   return m;
   
}

template<typename RealType>
void CreateMatrixSum(CRS<RealType> *matrix_out,
                     const RealType coeef_1,
                     const CRS<RealType> &matrix_1,
                     const RealType coeef_2,
                     const CRS<RealType> &matrix_2) {
   
   if (matrix_1.row_dim != matrix_2.row_dim || matrix_1.col_dim != matrix_2.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The summation of the matrices cannot be defined." << std::endl;
      ss << "matrix_1.row_dim = " << matrix_1.row_dim << ", matrix_1.col_dim = " << matrix_1.col_dim << std::endl;
      ss << "matrix_2.row_dim = " << matrix_2.row_dim << ", matrix_2.col_dim = " << matrix_2.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   *matrix_out = CRS<RealType>(matrix_1.row_dim, matrix_1.col_dim);
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
            matrix_out->val.push_back(coeef_1*matrix_1.val[j]);
            matrix_out->col.push_back(matrix_1.col[j]);
         }
      }
      else if (m1_count == 0 && m2_count != 0) {
         for (std::int64_t j = row_lower_2; j < row_upper_2; ++j) {
            matrix_out->val.push_back(coeef_2*matrix_2.val[j]);
            matrix_out->col.push_back(matrix_2.col[j]);
         }
      }
      else if (m1_count != 0 && m2_count != 0) {
         for (std::int64_t j = 0; j < m1_count + m2_count; ++j) {
            if (matrix_1.col[row_lower_1 + count_1] < matrix_2.col[row_lower_2 + count_2]) {
               matrix_out->val.push_back(coeef_1*matrix_1.val[row_lower_1 + count_1]);
               matrix_out->col.push_back(matrix_1.col[row_lower_1 + count_1]);
               count_1++;
               if (row_lower_1 + count_1 == row_upper_1) {
                  check = 1;
                  break;
               }
            }
            else if (matrix_1.col[row_lower_1 + count_1] == matrix_2.col[row_lower_2 + count_2]) {
               const RealType val = coeef_1*matrix_1.val[row_lower_1 + count_1] + coeef_2*matrix_2.val[row_lower_2 + count_2];
               if (std::abs(val) > 0.0) {
                  matrix_out->val.push_back(val);
                  matrix_out->col.push_back(matrix_1.col[row_lower_1 + count_1]);
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
               matrix_out->val.push_back(coeef_2*matrix_2.val[row_lower_2 + count_2]);
               matrix_out->col.push_back(matrix_2.col[row_lower_2 + count_2]);
               count_2++;
               if (row_lower_2 + count_2 == row_upper_2) {
                  check = 2;
                  break;
               }
            }
         }
         if (check == 1) {
            for (std::int64_t j = row_lower_2 + count_2; j < row_upper_2; ++j) {
               matrix_out->val.push_back(coeef_2*matrix_2.val[j]);
               matrix_out->col.push_back(matrix_2.col[j]);
            }
         }
         else if (check == 2) {
            for (std::int64_t j = row_lower_1 + count_1; j < row_upper_1; ++j) {
               matrix_out->val.push_back(coeef_1*matrix_1.val[j]);
               matrix_out->col.push_back(matrix_1.col[j]);
            }
         }
      }
      matrix_out->row[i + 1] = matrix_out->col.size();
   }
   
}




} // namespace sparse_matrix
} // namespace compnal


#endif /* COMPNAL_SPARSE_MATRIX_COMPRESSED_ROW_STORAGE_HPP_ */
