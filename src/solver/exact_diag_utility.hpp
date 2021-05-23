//
//  model_utility.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/24.
//

#ifndef model_utility_hpp
#define model_utility_hpp

#include "sparse_matrix.hpp"
#include <sstream>
#include <vector>
#include <algorithm>

namespace compnal {
namespace solver {

using CRS = sparse_matrix::CRS;

using BraketVector = sparse_matrix::BraketVector;


int64_t CalculateLocalBasis(int64_t global_basis, int64_t site, int64_t dim_onsite) {
   for (int64_t i = 0; i < site; ++i) {
      global_basis = global_basis/dim_onsite;
   }
   return global_basis%dim_onsite;
}

template<typename RealType>
BraketVector<RealType> CalculateMatrixVectorProduct(const CRS<RealType> &M, int site,
                                                    const BraketVector<RealType> &ket_in,
                                                    const std::vector<int64_t> &bases_in,
                                                    const std::vector<int64_t> &bases_out
                                                    ) {
   
   if (bases_in.size() != ket_in.size()) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      ss << "The sizes of bases and ket_in does not match each other" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   int64_t dim_onsite    = M.GetRowDim();
   int64_t dim_in        = bases_in.size();
   int64_t dim_out       = bases_out.size();
   int64_t site_constant = std::pow(dim_onsite, site);
   BraketVector<RealType> ket_out(dim_out);

#pragma omp parallel for
   for (int64_t i = 0; i < dim_out; ++i) {
      int64_t global_basis_out = bases_out[i];
      int64_t local_basis_out  = CalculateLocalBasis(global_basis_out, site, dim_onsite);
      int64_t begin            = M.Row(local_basis_out);
      int64_t end              = M.Row(local_basis_out + 1);
      RealType temp_val = 0.0;
      for (int64_t j = begin; j < end; ++j) {
         std::int64_t a_basis_out = global_basis_out - (local_basis_out - M.Col(j))*site_constant;
         auto iter_find = std::lower_bound(bases_in.begin(), bases_in.end(), a_basis_out);
         if (iter_find != bases_in.end() && *iter_find == a_basis_out) {
            auto inv_in = std::distance(bases_in.begin(), iter_find);
            temp_val += ket_in.Val(inv_in)*M.Val(j);
         }
      }
      ket_out.Val(i) = temp_val;
   }
   return ket_out;
}



} // namespace solver
} // namespace compnal


#endif /* model_utility_hpp */
