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
#include "exact_diag_utility.hpp"
#include <unordered_map>

namespace compnal {
namespace solver {

template<class ModelClass1D>
class ExactDiag1D {
   
   using RealType = typename ModelClass1D::ValueType;
   
   using CRS = sparse_matrix::CRS<RealType>;
   
   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   
   ModelClass1D model;
   
   explicit ExactDiag1D(const ModelClass1D &model_input): model(model_input) {}
   
   void SetBasis() {
      if (model.GetFlagRecalcBasis()) {
         model.GenerateBasis(&basis_, &basis_inv);
      }
   }
   
   
   
   void Diagonalize() {
      SetBasis();
      
   }
   
   const BraketVector &GetGSVector() const {
      return gs_vector_;
   }
   
   const BraketVector &GetBasis() const {
      return basis_;
   }
   
   RealType CalculateExpectationValue(const CRS &M, int site) {
      return gs_vector_*CalculateMatrixVectorProduct(M, site, gs_vector_, basis_, basis_);
   }
   
   
   
   
   
private:
   std::vector<int64_t> basis_;
   std::unordered_map<int64_t, int64_t> basis_inv;
   BraketVector gs_vector_;
   bool flag_recalc_gs = true;
   
   void GenerateMatrixElementsOnsite(ExactDiagMatrixElements<RealType> *edme,
                                     const int64_t basis,
                                     const int site,
                                     const CRS &matrix_onsite,
                                     const RealType coeef) const {
      
      if (std::abs(coeef) == 0.0) {
         return;
      }
      
      const int     basis_onsite  = edme->basis_onsite[site];
      const int64_t site_constant = edme->site_constant[site];
      
      for (int64_t i = matrix_onsite.Row(basis_onsite); i < matrix_onsite.Row(basis_onsite + 1); ++i) {
         const int64_t a_basis = basis + (matrix_onsite.Col(i) - basis_onsite)*site_constant;
         if (edme->inv_basis_affected.count(a_basis) == 0) {
            edme->inv_basis_affected[a_basis] = edme->basis_affected.size();
            edme->val.push_back(coeef*matrix_onsite.Val(i));
            edme->basis_affected.push_back(a_basis);
         }
         else {
            const RealType val = edme->val[edme->inv_basis_affected.at(a_basis)] + coeef*matrix_onsite.Val(i);
            if (std::abs(val) > edme->zero_precision) {
               edme->val[edme->inv_basis_affected.at(a_basis)] = val;
            }
         }
      }
   }
   
   void GenerateMatrixElementsIntersite(ExactDiagMatrixElements<RealType> *edme,
                                        const int64_t basis,
                                        const int site_1,
                                        const CRS &matrix_onsite_1,
                                        const int site_2,
                                        const CRS &matrix_onsite_2,
                                        const RealType coeef,
                                        const int fermion_sign) const {
      
      if (std::abs(coeef) == 0.0) {
         return;
      }
      
      const int     basis_onsite_1  = edme->basis_onsite[site_1];
      const int     basis_onsite_2  = edme->basis_onsite[site_2];
      const int64_t site_constant_1 = edme->site_constant[site_1];
      const int64_t site_constant_2 = edme->site_constant[site_1];

      for (int64_t i1 = matrix_onsite_1.Row(basis_onsite_1); i1 < matrix_onsite_1.Row(basis_onsite_1 + 1); ++i1) {
         const RealType val_1 = matrix_onsite_1.Val(i1);
         const int64_t  col_1 = matrix_onsite_1.Col(i1);
         for (int64_t i2 = matrix_onsite_2.Row(basis_onsite_2); i2 < matrix_onsite_2.Row(basis_onsite_2 + 1); ++i2) {
            const int64_t a_basis = (col_1 - basis_onsite_1)*site_constant_1 + (matrix_onsite_2.Col(i2) - basis_onsite_2)*site_constant_2;
            if (edme->inv_basis_affected.count(a_basis) == 0) {
               edme->inv_basis_affected[a_basis] = edme->basis_affected.size();
               edme->val.push_back(fermion_sign*coeef*val_1*matrix_onsite_2.Val(i2));
               edme->basis_affected.push_back(a_basis);
            }
            else {
               const RealType val = edme->val[edme->inv_basis_affected.at(a_basis)] + fermion_sign*coeef*val_1*matrix_onsite_2.Val(i2);
               if (std::abs(val) > edme->zero_precision) {
                  edme->val[edme->inv_basis_affected.at(a_basis)] = val;
               }
            }
         }
      }
   }
   
   int64_t CalculateLocalBasis(int64_t global_basis, int site, int dim_onsite) {
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
