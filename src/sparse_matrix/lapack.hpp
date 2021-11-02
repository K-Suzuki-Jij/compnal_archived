//
//  lapack.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef lapack_hpp
#define lapack_hpp

#include "compressed_row_storage.hpp"
#include "braket_vector.hpp"
#include <sstream>

namespace compnal {
namespace sparse_matrix {

extern "C" {
void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double **A,
            const int &LDA, double *W, double *work, const int &Lwork,
            int &INFO);
};

extern "C" {
void dstev_(const char &JOBZ, const int &N, double *D, double *E, double **Z,  const int &LDZ, double* WORK, int& INFO);
};


template <typename RealType>
void LapackDsyev(RealType *gs_value,
                 std::vector<RealType> *gs_vector,
                 const CRS<RealType> &matrix_in) {
   
   if (matrix_in.row_dim != matrix_in.col_dim || matrix_in.row_dim < 1 || matrix_in.col_dim < 1) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The input matrix is not a square one" << std::endl;
      ss << "row=" << matrix_in.row_dim << ", col=" << matrix_in.col_dim << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   double matrix_array[matrix_in.row_dim][matrix_in.col_dim];
   
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      for (std::size_t j = 0; j < matrix_in.col_dim; ++j) {
         matrix_array[i][j] = 0.0;
      }
   }
   
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
      for (std::size_t j = matrix_in.row[i]; j < matrix_in.row[i + 1]; ++j) {
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
   
   for (std::size_t i = 0; i < matrix_in.row_dim; ++i) {
       (*gs_vector)[i] = static_cast<RealType>(matrix_array[0][i]);
   }
   
   *gs_value = static_cast<RealType>(val_array[0]);
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
   
   for (std::size_t i = 0; i < dim; ++i) {
      (*gs_vector)[i] = static_cast<RealType>(Lap_Vec[0][i]);
   }
   
   *gs_value = static_cast<RealType>(Lap_D[0]);
}

}  // namespace sparse_matrix
}  // namespace compnal


#endif /* lapack_hpp */
