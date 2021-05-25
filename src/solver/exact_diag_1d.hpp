//
//  exact_diag_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#ifndef exact_diag_1d_hpp
#define exact_diag_1d_hpp

#include "sparse_matrix.hpp"
#include "model.hpp"
#include <unordered_map>

namespace compnal {
namespace solver {

template<typename ModelClass1D>
class ExactDiag1D {
   
   using RealType = typename ModelClass1D::ValueType;
   
   using CRS = sparse_matrix::CRS<RealType>;
   
   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   
   const ModelClass1D model;
   
   explicit ExactDiag1D(const ModelClass1D &model_input): model(model_input) {}
   
   const BraketVector &GetGSVector() const {
      return gs_vector_;
   }
   
   const BraketVector &GetBasis() const {
      return basis_;
   }
   
   RealType CalculateExpectationValue(const CRS &M, int site) {
      return gs_vector_*CalculateMatrixVectorProduct(M, site, gs_vector_, basis_, basis_);
   }
   
   std::vector<RealType> CalculateAllExpectationValues(const CRS &M) {
      int system_size = model.GetSystemSize();
      std::vector<RealType> out(system_size);
      for (int site = 0; site < system_size; ++site) {
         out[site] = gs_vector_*CalculateMatrixVectorProduct(M, site, gs_vector_, basis_, basis_);
      }
      return out;
   }
   
   
private:
   std::vector<int64_t> basis_;
   BraketVector gs_vector_;
   
   int64_t CalculateLocalBasis(int64_t global_basis, int64_t site, int64_t dim_onsite) {
      for (int64_t i = 0; i < site; ++i) {
         global_basis = global_basis/dim_onsite;
      }
      return global_basis%dim_onsite;
   }
   
   
   BraketVector CalculateMatrixVectorProduct(const CRS &M, int site,
                                             const BraketVector &ket_in,
                                             const std::vector<int64_t> &bases_in,
                                             const std::vector<int64_t> &bases_out
                                             ) {
      
      if (bases_in.size() != ket_in.GetDim()) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The sizes of bases and ket_in does not match each other" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      int64_t dim_onsite    = M.GetRowDim();
      int64_t dim_in        = bases_in.size();
      int64_t dim_out       = bases_out.size();
      int64_t site_constant = std::pow(dim_onsite, site);
      BraketVector ket_out(dim_out);
      
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
   
};

} // namespace solver
} // namespace compnal



#endif /* exact_diag_1d_hpp */
