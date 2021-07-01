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
class CRS {
   
public:
   CRS(): row_dim_(0), col_dim_(0), row_(1) {}
   CRS(const int64_t row_dim, const int64_t col_dim) {
      if (row_dim < 0 || col_dim < 0) {
         throw std::runtime_error("row_dim or col_dim of CRS must be larger than or equal to zero");
      }
      row_dim_ = row_dim;
      col_dim_ = col_dim;
      row_.resize(row_dim + 1);
#pragma omp parallel for
      for (int64_t i = 0; i < row_dim + 1; ++i) {
         row_[i] = 0;
      }
   }
   CRS(const std::vector<std::vector<RealType>> &mat_vec): row_dim_(mat_vec.size()) {
      col_dim_ = 0;
      row_.resize(row_dim_ + 1);
      row_[0] = 0;
      for (int64_t i = 0; i < row_dim_; ++i) {
         for (int64_t j = 0; j < mat_vec[i].size();++j) {
            if (mat_vec[i][j] != 0.0) {
               col_.push_back(j);
               val_.push_back(mat_vec[i][j]);
            }
            if (col_dim_ < j + 1) {
               col_dim_ = j + 1;
            }
         }
         row_[i + 1] = col_.size();
      }
   }
   
   void ResizeRow(const int64_t row_dim) {
      row_.resize(row_dim + 1);
      row_[0] = 0;
   }
   void ResizeColVal(const int64_t num_elements) {
      col_.resize(num_elements);
      val_.resize(num_elements);
   }
   
   void SetRowDim(const int64_t dim) { row_dim_ = dim; }
   void SetColDim(const int64_t dim) { col_dim_ = dim; }

   inline int64_t GetRowDim() const { return row_dim_; }
   inline int64_t GetColDim() const { return col_dim_; }
   
   inline int64_t   Row(const int64_t index) const { return row_[index]; }
   inline int64_t   Col(const int64_t index) const { return col_[index]; }
   inline RealType  Val(const int64_t index) const { return val_[index]; }
   
   inline int64_t  &Row(const int64_t index)       { return row_[index]; }
   inline int64_t  &Col(const int64_t index)       { return col_[index]; }
   inline RealType &Val(const int64_t index)       { return val_[index]; }
   
   inline int64_t GetNumElements() const { return val_.size(); }
   
   inline int64_t GetSizeVal() const { return val_.size(); }
   inline int64_t GetSizeCol() const { return col_.size(); }
   inline int64_t GetSizeRow() const { return row_.size(); }
         
   inline void PushRow(const int64_t  size) { row_.push_back(size); }
   inline void PushCol(const int64_t  col ) { col_.push_back(col) ; }
   inline void PushVal(const RealType val ) { val_.push_back(val) ; }
   
   const std::vector<int64_t> &GetRow() const {
      return row_;
   }
   
   const std::vector<int64_t> &GetCol() const {
      return col_;
   }
   
   const std::vector<RealType> &GetVal() const {
      return val_;
   }
      
   void Clear(const int64_t row_dim = 0, const int64_t col_dim = 0) {
      row_dim_ = row_dim;
      col_dim_ = col_dim;
      col_.clear();
      val_.clear();
      row_.resize(row_dim + 1);
      for (int64_t i = 0; i < row_dim + 1; ++i) {
         row_[i] = 0;
      }
   }
   
   void Free() {
      row_dim_ = 0;
      col_dim_ = 0;
      std::vector<int64_t>().swap(row_);
      std::vector<int64_t>().swap(col_);
      std::vector<RealType>().swap(val_);
      row_ = std::vector<int64_t>{0};
   }
   
   CRS CreateCopy(const RealType coeef = 1.0) const {
      const int64_t num_elements = GetNumElements();
      CRS matrix_out(row_dim_, col_dim_);
      matrix_out.ResizeColVal(num_elements);
#pragma omp parallel for
      for (int64_t i = 0; i < num_elements; ++i) {
         matrix_out.Col(i) = col_[i];
         matrix_out.Val(i) = coeef*val_[i];
      }
      
#pragma omp parallel for
      for (int64_t i = 0; i < row_dim_ + 1; ++i) {
         matrix_out.Row(i) = row_[i];
      }
      matrix_out.SetRowDim(row_dim_);
      matrix_out.SetColDim(col_dim_);
      return matrix_out;
   }
   
   void SortCol() {
#pragma omp parallel for schedule(guided)
      for (int64_t i = 0; i < row_dim_; ++i) {
         utility::QuickSort<int64_t, RealType>(&col_, &val_, row_[i], row_[i + 1]);
      }
   }
   
   void Print(const std::string display_name = "Matrix") const {
      for (int64_t i = 0; i < row_dim_; ++i) {
         for (int64_t j = row_.at(i); j < row_.at(i+1); ++j) {
            std::cout << display_name << "[";
            std::cout << std::noshowpos << std::left << std::setw(3) << i << "][";
            std::cout << std::left << std::setw(3) << col_[j] << "]=";
            std::cout << std::showpos << val_[j] << std::endl;
         }
      }
      std::cout << std::noshowpos;
   }
   
   void PrintInfo() const {
      std::cout << "Print information about CRS" << std::endl;
      std::cout << "row_dim = " << row_dim_ << std::endl;
      std::cout << "col_dim = " << col_dim_ << std::endl;
      for (int64_t i = 0; i < row_.size(); ++i) {
         std::cout << "row_[" << i << "] = " << row_.at(i) << std::endl;
      }
      for (int64_t i = 0; i < col_.size(); ++i) {
         std::cout << "col_[" << i << "] = " << col_.at(i) << std::endl;
      }
      for (int64_t i = 0; i < val_.size(); ++i) {
         std::cout << "val_[" << i << "] = " << val_.at(i) << std::endl;
      }
   }
   
   bool CheckSymmetric(const RealType threshold = 0.000000000000001/*pow(10,-15)*/) const {
      for (int64_t i = 0; i < row_dim_; ++i) {
         for (int64_t j = row_[i]; j < row_[i + 1]; ++j) {
            
            auto iter_begin = col_.begin() + row_[col_[j]];
            auto iter_end   = col_.begin() + row_[col_[j] + 1];
            auto iter_find  = std::lower_bound(iter_begin, iter_end, i);
            if (iter_find == iter_end || *iter_find != i) {
               std::cout << "The input matrix is not symmetric." << std::endl;
               std::cout << "Corresponding element does not exist." << std::endl;
               std::cout << "row=" << i << ", col=" << col_[j] << ", val=" << val_[j] << std::endl;
               return false;
            }
            auto inv = std::distance(iter_begin, iter_find);
            if (std::abs(val_[j] - val_[inv]) > threshold) {
               std::cout << "The input matrix is not symmetric." << std::endl;
               std::cout << "M[" << i << "][" << col_[j] << "]=" << val_[j] << ", " << val_[inv] << "=M[" << col_[j] << "][" << i << "]" << std::endl;
               return false;
            }
         }
      }
      return true;
   }
   
   
   
   CRS operator+() const { return *this; }
   CRS operator-() const { return CreateCopy(-1.0); }
   
   CRS& operator+=(const CRS &matrix_rhs) { return *this = CreateMatrixSum(*this, matrix_rhs); }
   CRS& operator-=(const CRS &matrix_rhs) { return *this = CreateMatrixSum(*this, matrix_rhs, 1.0, -1.0); }
   CRS& operator*=(const CRS &matrix_rhs) { return *this = CreateMatrixProduct(*this, matrix_rhs); }
   
   CRS& operator*=(const RealType coeef) {
      MultiplyMatrixByScalar(this, coeef);
      return *this;
   }
   
private:
   int64_t row_dim_;
   int64_t col_dim_;
   std::vector<int64_t>  row_;
   std::vector<int64_t>  col_;
   std::vector<RealType> val_;
};

template<typename RealType>
bool operator==(const CRS<RealType> &matrix_lhs, const CRS<RealType> &matrix_rhs) {
   if(matrix_lhs.GetRowDim() != matrix_rhs.GetRowDim()) {
      return false;
   }
   if (matrix_lhs.GetColDim() != matrix_rhs.GetColDim()) {
      return false;
   }
   if (matrix_lhs.GetSizeCol() != matrix_rhs.GetSizeCol()) {
      return false;
   }
   const int64_t row_dim = matrix_lhs.GetRowDim();
   for (int64_t i = 0; i < row_dim + 1; ++i) {
      if (matrix_lhs.Row(i) != matrix_rhs.Row(i)) {
         return false;
      }
   }
   const int64_t num_elements = matrix_lhs.GetSizeCol();
   for (int64_t i = 0; i < num_elements; ++i) {
      if (matrix_lhs.Col(i) != matrix_rhs.Col(i)) {
         return false;
      }
      if (matrix_lhs.Val(i) != matrix_rhs.Val(i)) {
         return false;
      }
   }
   return true;
}

template<typename RealType>
CRS<RealType> operator+(const CRS<RealType> &matrix_lhs, const CRS<RealType> &matrix_rhs) {
   return CRS<RealType>(matrix_lhs) += matrix_rhs;
}

template<typename RealType>
CRS<RealType> operator-(const CRS<RealType> &matrix_lhs, const CRS<RealType> &matrix_rhs) {
   return CRS<RealType>(matrix_lhs) -= matrix_rhs;
}

template<typename RealType>
CRS<RealType> operator*(const CRS<RealType> &matrix_lhs, const CRS<RealType> &matrix_rhs) {
   return CRS<RealType>(matrix_lhs) *= matrix_rhs;
}

template<typename RealType>
CRS<RealType> operator*(const RealType coeef_lhs, const CRS<RealType> &matrix_rhs) {
   return CRS<RealType>(matrix_rhs) *= coeef_lhs;
}

template<typename RealType>
CRS<RealType> operator*(const CRS<RealType> &matrix_lhs, const RealType coeef_rhs) {
   return CRS<RealType>(matrix_lhs) *= coeef_rhs;
}

template<typename RealType>
CRS<RealType> CreateTransposedMatrix(const CRS<RealType> &matrix_in) {
   
   CRS<RealType> matrix_out;
   
   const int64_t row_dim = matrix_in.GetRowDim();
   const int64_t col_dim = matrix_in.GetColDim();
   
   std::vector<int64_t> row_count(row_dim);
   
   for (int64_t i = 0; i < col_dim; ++i) {
      for (int64_t j = 0; j < row_dim; ++j) {
         int64_t row = matrix_in.Row(j) + row_count[j];
         if (row < matrix_in.Row(j + 1) && matrix_in.Col(row) == i) {
            matrix_out.PushVal(matrix_in.Val(row));
            matrix_out.PushCol(j);
            row_count[j]++;
         }
      }
      matrix_out.UpdateRow();
   }
   matrix_out.SetRowDim(col_dim);
   matrix_out.SetColDim(row_dim);
   
   return matrix_out;
}

template<typename RealType>
CRS<RealType> CreateMatrixProduct(const CRS<RealType> &matrix_lhs, const CRS<RealType> &matrix_rhs, const RealType coeef_lhs = 1.0, const RealType coeef_rhs = 1.0) {
   
   const int64_t row_dim_1 = matrix_lhs.GetRowDim();
   const int64_t col_dim_1 = matrix_lhs.GetColDim();
   const int64_t row_dim_2 = matrix_rhs.GetRowDim();
   const int64_t col_dim_2 = matrix_rhs.GetColDim();
   
   if (col_dim_1 != row_dim_2) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Matrix product cannot be defined" << std::endl;
      ss << "col_dim_1 = " << col_dim_1 << ", row_dim_2 = " << row_dim_2 << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   CRS<RealType> matrix_out(row_dim_1, col_dim_2);
   
   std::vector<RealType> temp_v1(col_dim_1, 0.0);
   std::vector<RealType> temp_v2(col_dim_2, 0.0);
   
   for (int64_t i = 0; i < row_dim_1; ++i) {
      const int64_t begin_m1 = matrix_lhs.Row(i);
      const int64_t end_m1   = matrix_lhs.Row(i + 1);
      for (int64_t j = begin_m1; j < end_m1; ++j) {
         temp_v1[matrix_lhs.Col(j)] = coeef_lhs*matrix_lhs.Val(j);
      }
      
      for (int64_t j = 0; j < col_dim_1; ++j) {
         const int64_t begin = matrix_rhs.Row(j);
         const int64_t end   = matrix_rhs.Row(j + 1);
         for (int64_t k = begin; k < end; ++k) {
            temp_v2[matrix_rhs.Col(k)] += temp_v1[j]*coeef_rhs*matrix_rhs.Val(k);
         }
      }
      
      for (int64_t j = 0; j < col_dim_2; ++j) {
         if (std::abs(temp_v2[j]) > 0.0) {
            matrix_out.PushVal(temp_v2[j]);
            matrix_out.PushCol(j);
         }
      }
      
      matrix_out.Row(i+1) = matrix_out.GetSizeCol();
      
      for (int64_t j = begin_m1; j < end_m1; ++j) {
         temp_v1[matrix_lhs.Col(j)] = 0.0;
      }

      const int64_t begin_out = matrix_out.Row(i);
      const int64_t end_out   = matrix_out.Row(i + 1);
      for (int64_t j = begin_out; j < end_out; ++j) {
         temp_v2[matrix_out.Col(j)] = 0.0;
      }
      
   }
   
   return matrix_out;
}

template <typename RealType>
void MultiplyMatrixByScalar(CRS<RealType> *matrix, const RealType coeef) {
   if (coeef == 0.0) {
      matrix->Clear(matrix->GetRowDim(), matrix->GetColDim());
   }
   else {
      int64_t dim = matrix->GetRowDim();
#pragma omp parallel for
      for (int64_t i = 0; i < dim; ++i) {
         const int64_t begin = matrix->Row(i);
         const int64_t end   = matrix->Row(i + 1);
         for (int64_t j = begin; j < end; ++j) {
            matrix->Val(j) *= coeef;
         }
      }
   }
}

template<typename RealType>
CRS<RealType> CreateMatrixSum(const CRS<RealType> &matrix_1, const CRS<RealType> &matrix_2, const RealType coeef_1 = 1.0, const RealType coeef_2 = 1.0) {
   
   const int64_t row_dim_1 = matrix_1.GetRowDim();
   const int64_t col_dim_1 = matrix_1.GetColDim();
   
   const int64_t row_dim_2 = matrix_2.GetRowDim();
   const int64_t col_dim_2 = matrix_2.GetColDim();
   
   if (row_dim_1 != row_dim_2 || col_dim_1 != col_dim_2) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "Matrix types do not match each other" << std::endl;
      ss << "row_dim_1 = " << row_dim_1 << ", col_dim_1 = " << col_dim_1 << std::endl;
      ss << "row_dim_2 = " << row_dim_2 << ", col_dim_2 = " << col_dim_2 << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   CRS<RealType> matrix_out(row_dim_1, col_dim_1);
   int64_t total_count;
   
   total_count = 0;
   for (int64_t i = 0; i < row_dim_1; ++i) {
      
      int check = 0;
      int64_t count_1 = 0;
      int64_t count_2 = 0;

      const int64_t row_lower_1 = matrix_1.Row(  i  );
      const int64_t row_upper_1 = matrix_1.Row(i + 1);
      const int64_t row_lower_2 = matrix_2.Row(  i  );
      const int64_t row_upper_2 = matrix_2.Row(i + 1);

      const int64_t m1_count = row_upper_1 - row_lower_1;
      const int64_t m2_count = row_upper_2 - row_lower_2;
      
      if (m1_count != 0 && m2_count == 0) {
         for (int64_t j = row_lower_1; j < row_upper_1; ++j) {
            matrix_out.PushVal(coeef_1*matrix_1.Val(j));
            matrix_out.PushCol(matrix_1.Col(j));
         }
      }
      else if (m1_count == 0 && m2_count != 0) {
         for (int64_t j = row_lower_2; j < row_upper_2; ++j) {
            matrix_out.PushVal(coeef_2*matrix_2.Val(j));
            matrix_out.PushCol(matrix_2.Col(j));
         }
      }
      else if (m1_count != 0 && m2_count != 0) {
         for (int64_t j = 0; j < m1_count + m2_count; ++j) {
            if (matrix_1.Col(row_lower_1 + count_1) < matrix_2.Col(row_lower_2 + count_2)) {
               matrix_out.PushVal(coeef_1*matrix_1.Val(row_lower_1 + count_1));
               matrix_out.PushCol(matrix_1.Col(row_lower_1 + count_1));
               count_1++;
               if (row_lower_1 + count_1 == row_upper_1) {
                  check = 1;
                  break;
               }
            }
            else if (matrix_1.Col(row_lower_1 + count_1) == matrix_2.Col(row_lower_2 + count_2)) {
               const RealType val = coeef_1*matrix_1.Val(row_lower_1 + count_1) + coeef_2*matrix_2.Val(row_lower_2 + count_2);
               if (std::abs(val) > 0.0) {
                  matrix_out.PushVal(val);
                  matrix_out.PushCol(matrix_1.Col(row_lower_1 + count_1));
                  count_1++;
                  count_2++;
                  int64_t temp_count_1 = row_lower_1 + count_1;
                  int64_t temp_count_2 = row_lower_2 + count_2;
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
                  int64_t temp_count_1 = row_lower_1 + count_1;
                  int64_t temp_count_2 = row_lower_2 + count_2;
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
               matrix_out.PushVal(coeef_2*matrix_2.Val(row_lower_2 + count_2));
               matrix_out.PushCol(matrix_2.Col(row_lower_2 + count_2));
               count_2++;
               if (row_lower_2 + count_2 == row_upper_2) {
                  check = 2;
                  break;
               }
            }
         }
         if (check == 1) {
            for (int64_t j = row_lower_2 + count_2; j < row_upper_2; ++j) {
               matrix_out.PushVal(coeef_2*matrix_2.Val(j));
               matrix_out.PushCol(matrix_2.Col(j));
            }
         }
         else if (check == 2) {
            for (int64_t j = row_lower_1 + count_1; j < row_upper_1; ++j) {
               matrix_out.PushVal(coeef_1*matrix_1.Val(j));
               matrix_out.PushCol(matrix_1.Col(j));
            }
         }
      }
      matrix_out.Row(i + 1) = matrix_out.GetSizeCol();
   }
   return matrix_out;
}

} // namespace sparse_matrix
} // namespace compnal


#endif /* compressed_row_storage_hpp */
