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

#ifndef COMPNAL_BLAS_LAPACK_HPP_
#define COMPNAL_BLAS_LAPACK_HPP_

#include "compressed_row_storage.hpp"
#include <sstream>

namespace compnal {
namespace blas {

extern "C" {
void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double **A, const int &LDA, double *W, double *work, const int &Lwork, int &INFO);
};

extern "C" {
void dstev_(const char &JOBZ, const int &N, double *D, double *E, double **Z,  const int &LDZ, double *WORK, int &INFO);
};

extern "C" {
void dspgv_(const int &ITYPE,const char &JOBZ, const char &UPLO, const int &N, double *AP, double *BP, double *W, double **Z, const int &LDZ, double *WORK, int &INFO);
};

template <typename RealType>
void LapackDsyev(RealType *gs_value,
                 std::vector<RealType> *gs_vector,
                 const CRS<RealType> &matrix_in,
                 const int target_level = 0
                 ) {
   
   if (matrix_in.row_dim != matrix_in.col_dim || matrix_in.row_dim < 1 || matrix_in.col_dim < 1) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in.row_dim << ", col=" << matrix_in.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (target_level < 0) {
      std::stringstream ss;
      ss << "Invalid target_level" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   double matrix_array[matrix_in.row_dim][matrix_in.col_dim];
   
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      for (std::int64_t j = 0; j < matrix_in.col_dim; ++j) {
         matrix_array[i][j] = 0.0;
      }
   }
   
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
      for (std::int64_t j = matrix_in.row[i]; j < matrix_in.row[i + 1]; ++j) {
         matrix_array[i][matrix_in.col[j]] = static_cast<double>(matrix_in.val[j]);
         matrix_array[matrix_in.col[j]][i] = static_cast<double>(matrix_in.val[j]);
      }
   }
   
   int info;
   double val_array[matrix_in.row_dim];
   double work[3 * matrix_in.row_dim];
   
   dsyev_('V', 'L', static_cast<int>(matrix_in.row_dim),
          (double**)matrix_array,
          static_cast<int>(matrix_in.row_dim), val_array, work,
          static_cast<int>(3 * matrix_in.row_dim), info);
   
   gs_vector->resize(matrix_in.row_dim);
   
   for (std::int64_t i = 0; i < matrix_in.row_dim; ++i) {
       (*gs_vector)[i] = static_cast<RealType>(matrix_array[target_level][i]);
   }
   
   *gs_value = static_cast<RealType>(val_array[target_level]);
}

template <typename RealType>
void LapackDstev(RealType *gs_value,
                 std::vector<RealType> *gs_vector,
                 const std::vector<RealType> &diag,
                 const std::vector<RealType> &off_diag) {
   
   if (off_diag.size() + 1 != diag.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "diag size=" << diag.size() << ", off_diag size=" << off_diag.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   int dim = static_cast<int>(diag.size());

   int info;
   double Lap_D[dim];
   double Lap_E[dim - 1];
   double Lap_Vec[dim][dim];
   double Lap_Work[2*dim];
   
   for (int i = 0; i < dim; i++) {
      Lap_D[i] = static_cast<double>(diag[i]);
   }
   
   for (int i = 0; i < dim - 1; i++) {
      Lap_E[i] = static_cast<double>(off_diag[i]);
   }
   
   dstev_('V', dim, Lap_D, Lap_E, (double**)Lap_Vec, dim, Lap_Work, info);
   
   gs_vector->resize(dim);
   
   for (int i = 0; i < dim; ++i) {
      (*gs_vector)[i] = static_cast<RealType>(Lap_Vec[0][i]);
   }
   
   *gs_value = static_cast<RealType>(Lap_D[0]);
}

template <typename RealType>
void LapackDspgv(std::vector<RealType> *eigenvalues,
                 std::vector<RealType> *eigenvectors,
                 const int i_type,
                 const int dim,
                 const std::vector<RealType> &mat_a,
                 const std::vector<RealType> &mat_b) {
   
   const int size = dim*(dim + 1)/2;
   
   if (static_cast<int>(mat_a.size()) < size || static_cast<int>(mat_b.size()) < size) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   double lap_ap[size];
   double lap_bp[size];
   double lap_w[dim];
   double lap_z[dim][dim];
   double lap_work[3*dim];
   int lap_ldz = dim;
   int lap_info;
   
   for (int i = 0; i < size; ++i) {
      lap_ap[i] = static_cast<double>(mat_a[i]);
      lap_bp[i] = static_cast<double>(mat_b[i]);
   }
   
   dspgv_(i_type, 'V', 'U', dim, lap_ap, lap_bp, lap_w, (double**)lap_z, lap_ldz, lap_work, lap_info);
   
   eigenvalues->resize(dim);
   eigenvectors->resize(dim*dim);
   
   for (int i = 0; i < dim; ++i) {
      (*eigenvalues)[i] = static_cast<RealType>(lap_w[i]);
      for (int j = 0; j < dim; ++j) {
         (*eigenvectors)[i*dim + j] = static_cast<RealType>(lap_z[i][j]);
      }
   }
   
}

} // namespace blas
} // namespace compnel

#endif /* COMPNAL_BLAS_LAPACK_HPP_ */
