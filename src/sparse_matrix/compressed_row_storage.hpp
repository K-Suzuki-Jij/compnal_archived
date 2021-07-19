//
//  compressed_row_storage.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef compressed_row_storage_hpp
#define compressed_row_storage_hpp

#include "../utility/utility.hpp"
#include <iostream>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <sstream>

namespace compnal {
namespace sparse_matrix {

template<typename RealType>
struct CRS {
  
   //-----------Member Variables-----------
   int64_t row_dim = 0;
   int64_t col_dim = 0;
   std::vector<int64_t>  row;
   std::vector<int64_t>  col;
   std::vector<RealType> val;
   //--------------------------------------
   
   //-------------Constructor--------------
   CRS() {}
   CRS(const int64_t row_dim_in, const int64_t col_dim_in) {
      if (row_dim_in < 0 || col_dim_in < 0) {
         throw std::runtime_error("row_dim and col_dim of CRS must be larger than or equal to zero");
      }
      this->row_dim = row_dim_in;
      this->col_dim = col_dim_in;
      this->row.resize(row_dim_in + 1);
#pragma omp parallel for
      for (int64_t i = 0; i <= row_dim_in; ++i) {
         this->row[i] = 0;
      }
   }
   CRS(const std::vector<std::vector<RealType>> &mat_vec) {
      this->row_dim = mat_vec.size();
      this->col_dim = 0;
      this->row.resize(this->row_dim + 1);
      this->row[0] = 0;
      for (int64_t i = 0; i < this->row_dim; ++i) {
         for (int64_t j = 0; j < mat_vec[i].size();++j) {
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
   //--------------------------------------
   
   //----------Copy Constructor------------
   CRS(const CRS &matrix) { Assign(matrix); }
   //--------------------------------------
   
   //-----------Member Functions-----------
   void Assign(const CRS &matrix) {
      this->row_dim = matrix.row_dim;
      this->col_dim = matrix.col_dim;
      this->row.resize(matrix.row_dim + 1);
      
      const int64_t num_interaction = static_cast<int64_t>(matrix.col.size());
      this->col.resize(num_interaction);
      this->val.resize(num_interaction);
      
#pragma omp parallel for
      for (int64_t i = 0; i <= matrix.row_dim; ++i) {
         this->row[i] = matrix.row[i];
      }
      
#pragma omp parallel for
      for (int64_t i = 0; i < num_interaction; ++i) {
         this->col[i] = matrix.col[i];
         this->val[i] = matrix.val[i];
      }
   }
   
   void Free() {
      this->row_dim = 0; this->col_dim = 0;
      std::vector<int64_t> ().swap(this->row);
      std::vector<int64_t> ().swap(this->col);
      std::vector<RealType>().swap(this->val);
   }
   
   void Clear() {
      this->row_dim = 0; this->col_dim = 0;
      this->row.clear(); this->col.clear(); this->val.clear();
   }
   
   void SortCol() {
#pragma omp parallel for schedule(guided)
      for (int64_t i = 0; i < this->row_dim; ++i) {
         utility::QuickSort<int64_t, RealType>(&this->col, &this->val, this->row[i], this->row[i + 1]);
      }
   }
   
   void PrintMatrix(const std::string display_name = "Matrix") const {
      for (int64_t i = 0; i < this->row_dim; ++i) {
         for (int64_t j = this->row.at(i); j < this->row.at(i+1); ++j) {
            std::cout << display_name << "[";
            std::cout << std::noshowpos << std::left << std::setw(3) << i << "][";
            std::cout << std::left << std::setw(3) << this->col[j] << "]=";
            std::cout << std::showpos << this->val[j] << std::endl;
         }
      }
      std::cout << std::noshowpos;
   }
   
   void PrintInfo() const {
      std::cout << "Print information about CRS" << std::endl;
      std::cout << "row_dim = " << this->row_dim << std::endl;
      std::cout << "col_dim = " << this->col_dim << std::endl;
      for (int64_t i = 0; i < this->row.size(); ++i) {
         std::cout << "row[" << i << "] = " << this->row.at(i) << std::endl;
      }
      for (int64_t i = 0; i < this->col.size(); ++i) {
         std::cout << "col[" << i << "] = " << this->col.at(i) << std::endl;
      }
      for (int64_t i = 0; i < this->val.size(); ++i) {
         std::cout << "val[" << i << "] = " << this->val.at(i) << std::endl;
      }
   }
   
   bool isSymmetric(const RealType threshold = 0.000000000000001/*pow(10,-15)*/) const {
      for (int64_t i = 0; i < this->row_dim; ++i) {
         for (int64_t j = this->row[i]; j < this->row[i + 1]; ++j) {
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
   //--------------------------------------
   
   //---------Operator Overloading---------
   CRS operator+() const { return *this; }
   CRS operator-() const {
      CRS matrix(this->row_dim, this->col_dim);
#pragma omp parallel for
      for (int64_t i = 1; i <= this->row_dim; ++i) {
         matrix.row[i] = this->row[i];
      }
      int64_t num_interaction = static_cast<int64_t>(this->col.size());
      matrix.col.resize(num_interaction);
      matrix.val.resize(num_interaction);
#pragma omp parallel for
      for (int64_t i = 0; i < num_interaction; ++i) {
         matrix.col[i] = this->col[i];
         matrix.val[i] = -1.0*this->val[i];
      }
      return matrix;
   }
   CRS& operator =(const CRS &matrix) {Assign(matrix); return *this; }
   CRS& operator+=(const CRS &matrix) { return *this = CreateMatrixSum(*this, matrix); }
   CRS& operator-=(const CRS &matrix) { return *this = CreateMatrixSum(*this, matrix, +1.0, -1.0); }
   CRS& operator*=(const CRS &matrix) { return *this = CreateMatrixProduct(*this, matrix); }
   CRS& operator*=(const RealType coeef) { MultiplyMatrixByScalar(this, coeef); return *this; }
   
};

template<typename RealType> CRS<RealType> operator+(const CRS<RealType> &matrix_1, const CRS<RealType> &matrix_2) { return CreateMatrixSum(matrix_1, matrix_2); }
template<typename RealType> CRS<RealType> operator-(const CRS<RealType> &matrix_1, const CRS<RealType> &matrix_2) { return CreateMatrixSum(matrix_1, matrix_2, 1.0, -1.0); }
template<typename RealType> CRS<RealType> operator*(const CRS<RealType> &matrix_1, const CRS<RealType> &matrix_2) { return CreateMatrixProduct(matrix_1, matrix_2); }
template<typename RealType>
CRS<RealType> operator*(const CRS<RealType> &matrix, const RealType coeef) {
   CRS<RealType> matrix_out(matrix);
   MultiplyMatrixByScalar(&matrix_out, coeef);
   return matrix_out;
}
template<typename RealType>
CRS<RealType> operator*(const RealType coeef, const CRS<RealType> &matrix) { return matrix*coeef; }

template <typename RealType>
void MultiplyMatrixByScalar(CRS<RealType> *matrix, const RealType coeef) {
   if (coeef == 0.0) {
      for (int64_t i = 0; i < matrix->row_dim; ++i) {
         matrix->row[i + 1] = 0;
      }
      matrix->col.clear();
      matrix->val.clear();
   }
   else {
      int64_t num_interaction = static_cast<int64_t>(matrix->val.size());
#pragma omp parallel for
      for (int64_t i = 0; i < num_interaction; ++i) {
         matrix->val[i] *= coeef;
      }
   }
}

template<typename RealType>
CRS<RealType> CreateTransposedMatrix(const CRS<RealType> &matrix_in) {
   
   CRS<RealType> matrix_out;
   
   std::vector<int64_t> row_count(matrix_in.row_dim);
   
   matrix_out.row.push_back(0);
   for (int64_t i = 0; i < matrix_in.col_dim; ++i) {
      for (int64_t j = 0; j < matrix_in.row_dim; ++j) {
         const int64_t row = matrix_in.row[j] + row_count[j];
         if (row < matrix_in.row[j + 1] && matrix_in.col[row] == i) {
            matrix_out.val.push_back(matrix_in.Val(row));
            matrix_out.col.push_back(j);
            row_count[j]++;
         }
      }
      matrix_out.row.push_back(matrix_out.col.size());
   }
   matrix_out.row_dim = matrix_in.col_dim;
   matrix_out.col_dim = matrix_in.row_dim;
   
   return matrix_out;
}

template<typename RealType>
CRS<RealType> CreateMatrixProduct(const CRS<RealType> &matrix_lhs,
                                  const CRS<RealType> &matrix_rhs,
                                  const RealType coeef_lhs = 1.0,
                                  const RealType coeef_rhs = 1.0) {
      
   if (matrix_lhs.col_dim != matrix_rhs.row_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Matrix product cannot be defined" << std::endl;
      ss << "matrix_lhs.col_dim = " << matrix_lhs.col_dim << ", matrix_rhs.row_dim = " << matrix_rhs.row_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   CRS<RealType> matrix_out;
   matrix_out.row.resize(matrix_lhs.row_dim + 1);
   
   std::vector<RealType> temp_v1(matrix_lhs.col_dim, 0.0);
   std::vector<RealType> temp_v2(matrix_rhs.col_dim, 0.0);
   
   matrix_out.row[0] = 0;
   for (int64_t i = 0; i < matrix_lhs.row_dim; ++i) {
      for (int64_t j = matrix_lhs.row[i]; j < matrix_lhs.row[i + 1]; ++j) {
         temp_v1[matrix_lhs.Col(j)] = coeef_lhs*matrix_lhs.Val(j);
      }
      for (int64_t j = 0; j < matrix_lhs.col_dim; ++j) {
         for (int64_t k = matrix_rhs.row[j]; k < matrix_rhs.row[j + 1]; ++k) {
            temp_v2[matrix_rhs.col[k]] += temp_v1[j]*coeef_rhs*matrix_rhs.val[k];
         }
      }
      for (int64_t j = 0; j < matrix_rhs.col_dim; ++j) {
         if (std::abs(temp_v2[j]) > 0.0) {
            matrix_out.val.push_back(temp_v2[j]);
            matrix_out.col.push_back(j);
         }
      }
      
      matrix_out.row[i + 1] = matrix_out.col.size();
      
      for (int64_t j = matrix_lhs.row[i]; j < matrix_lhs.row[i + 1]; ++j) {
         temp_v1[matrix_lhs.col[j]] = 0.0;
      }
      for (int64_t j = matrix_out.row[i]; j < matrix_out.row[i + 1]; ++j) {
         temp_v2[matrix_out.col[j]] = 0.0;
      }
   }
   
   return matrix_out;
}

template<typename RealType>
CRS<RealType> CreateMatrixSum(const CRS<RealType> &matrix_1,
                              const CRS<RealType> &matrix_2,
                              const RealType coeef_1 = 1.0,
                              const RealType coeef_2 = 1.0) {
   
   if (matrix_1.row_dim != matrix_2.row_dim || matrix_1.col_dim != matrix_2.col_dim) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The summation of the matrices cannot be defined." << std::endl;
      ss << "matrix_1.row_dim = " << matrix_1.row_dim << ", matrix_1.col_dim = " << matrix_1.col_dim << std::endl;
      ss << "matrix_2.row_dim = " << matrix_2.row_dim << ", matrix_2.col_dim = " << matrix_2.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   CRS<RealType> matrix_out(matrix_1.row_dim, matrix_1.col_dim);
   int64_t total_count;
   
   total_count = 0;
   for (int64_t i = 0; i < matrix_1.row_dim; ++i) {
      
      int check = 0;
      int64_t count_1 = 0;
      int64_t count_2 = 0;

      const int64_t row_lower_1 = matrix_1.row[  i  ];
      const int64_t row_upper_1 = matrix_1.row[i + 1];
      const int64_t row_lower_2 = matrix_2.row[  i  ];
      const int64_t row_upper_2 = matrix_2.row[i + 1];

      const int64_t m1_count = row_upper_1 - row_lower_1;
      const int64_t m2_count = row_upper_2 - row_lower_2;
      
      if (m1_count != 0 && m2_count == 0) {
         for (int64_t j = row_lower_1; j < row_upper_1; ++j) {
            matrix_out.val.push_back(coeef_1*matrix_1.val[j]);
            matrix_out.col.push_back(matrix_1.col[j]);
         }
      }
      else if (m1_count == 0 && m2_count != 0) {
         for (int64_t j = row_lower_2; j < row_upper_2; ++j) {
            matrix_out.val.push_back(coeef_2*matrix_2.val[j]);
            matrix_out.col.push_back(matrix_2.col[j]);
         }
      }
      else if (m1_count != 0 && m2_count != 0) {
         for (int64_t j = 0; j < m1_count + m2_count; ++j) {
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
               const RealType val = coeef_1*matrix_1.val[row_lower_1 + count_1] + coeef_2*matrix_2.val[row_lower_2 + count_2];
               if (std::abs(val) > 0.0) {
                  matrix_out.val.push_back[val];
                  matrix_out.col.push_back[matrix_1.col[row_lower_1 + count_1]];
                  count_1++;
                  count_2++;
                  const int64_t temp_count_1 = row_lower_1 + count_1;
                  const int64_t temp_count_2 = row_lower_2 + count_2;
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
                  const int64_t temp_count_1 = row_lower_1 + count_1;
                  const int64_t temp_count_2 = row_lower_2 + count_2;
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
               matrix_out.val.push_back[coeef_2*matrix_2.val[row_lower_2 + count_2]];
               matrix_out.col.push_back[matrix_2.col[row_lower_2 + count_2]];
               count_2++;
               if (row_lower_2 + count_2 == row_upper_2) {
                  check = 2;
                  break;
               }
            }
         }
         if (check == 1) {
            for (int64_t j = row_lower_2 + count_2; j < row_upper_2; ++j) {
               matrix_out.val.push_back[coeef_2*matrix_2.val[j]];
               matrix_out.col.push_back[matrix_2.col[j]];
            }
         }
         else if (check == 2) {
            for (int64_t j = row_lower_1 + count_1; j < row_upper_1; ++j) {
               matrix_out.val.push_back[coeef_1*matrix_1.val[j]];
               matrix_out.col.push_back[matrix_1.col[j]];
            }
         }
      }
      matrix_out.row[i + 1] = matrix_out.col.size();
   }
   return matrix_out;
}




} // namespace sparse_matrix
} // namespace compnal


#endif /* compressed_row_storage_hpp */
