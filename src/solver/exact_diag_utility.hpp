//
//  exact_diag_utility.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef exact_diag_utility_hpp
#define exact_diag_utility_hpp

#include "model.hpp"
#include "sparse_matrix.hpp"

#include <vector>
#include <unordered_map>

namespace compnal {
namespace solver {

template<typename RealType>
struct ExactDiagMatrixElements {
   std::vector<RealType> val;
   std::vector<int64_t>  basis_affected;
   std::vector<int>      basis_onsite;
   std::vector<int64_t>  site_constant;
   std::unordered_map<int64_t, int64_t> inv_basis_affected;
   RealType zero_precision = 0.0;
};

int CalculateLocalBasis(int64_t global_basis, const int site, const int dim_onsite) {
   for (int64_t i = 0; i < site; ++i) {
      global_basis = global_basis/dim_onsite;
   }
   return static_cast<int>(global_basis%dim_onsite);
}

template<typename RealType>
void GenerateMatrixElementsOnsite(ExactDiagMatrixElements<RealType> *edme,
                                  const int64_t basis,
                                  const int site,
                                  const sparse_matrix::CRS<RealType> &matrix_onsite,
                                  const RealType coeef) {
   
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

template<typename RealType>
void GenerateMatrixElementsIntersite(ExactDiagMatrixElements<RealType> *edme,
                                     const int64_t basis,
                                     const int site_1,
                                     const sparse_matrix::CRS<RealType> &matrix_onsite_1,
                                     const int site_2,
                                     const sparse_matrix::CRS<RealType> &matrix_onsite_2,
                                     const RealType coeef,
                                     const int fermion_sign) {
   
   if (std::abs(coeef) == 0.0) {
      return;
   }
   
   const int     basis_onsite_1  = edme->basis_onsite[site_1];
   const int     basis_onsite_2  = edme->basis_onsite[site_2];
   const int64_t site_constant_1 = edme->site_constant[site_1];
   const int64_t site_constant_2 = edme->site_constant[site_2];

   for (int64_t i1 = matrix_onsite_1.Row(basis_onsite_1); i1 < matrix_onsite_1.Row(basis_onsite_1 + 1); ++i1) {
      const RealType val_1 = matrix_onsite_1.Val(i1);
      const int64_t  col_1 = matrix_onsite_1.Col(i1);
      for (int64_t i2 = matrix_onsite_2.Row(basis_onsite_2); i2 < matrix_onsite_2.Row(basis_onsite_2 + 1); ++i2) {
         const int64_t a_basis = basis + (col_1 - basis_onsite_1)*site_constant_1 + (matrix_onsite_2.Col(i2) - basis_onsite_2)*site_constant_2;
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

template<typename RealType>
void GenerateMatrixElements(ExactDiagMatrixElements<RealType> *edme, const int64_t basis, const model::Heisenberg1D<RealType> &model) {
   
   for (int site = 0; site < model.GetSystemSize(); ++site) {
      edme->basis_onsite[site] = CalculateLocalBasis(basis, site, model.GetDimOnsite());
   }
   
   //Onsite elements
   for (int site = 0; site < model.GetSystemSize(); ++site) {
      GenerateMatrixElementsOnsite(edme, basis, site, model.GetOperatorHam(), 1.0);
   }
   
   //Intersite elements SzSz
   for (int distance = 1; distance <= model.GetJz().size(); ++distance) {
      for (int site = 0; site < model.GetSystemSize() - distance; ++site) {
         GenerateMatrixElementsIntersite(edme, basis, site, model.GetOperatorSz(), site + distance, model.GetOperatorSz(), model.GetJz(distance - 1), 1.0);
      }
   }
   
   //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
   for (int distance = 1; distance <= model.GetJxy().size(); ++distance) {
      for (int site = 0; site < model.GetSystemSize() - distance; ++site) {
         GenerateMatrixElementsIntersite(edme, basis, site, model.GetOperatorSp(), site + distance, model.GetOperatorSm(), 0.5*model.GetJxy(distance - 1), 1.0);
         GenerateMatrixElementsIntersite(edme, basis, site, model.GetOperatorSm(), site + distance, model.GetOperatorSp(), 0.5*model.GetJxy(distance - 1), 1.0);
      }
   }
   
   if (model.GetBoundaryCondition() == model::BoundaryCondition::PBC) {
      //Intersite elements SzSz
      for (int distance = 1; distance <= model.GetJz().size(); ++distance) {
         for (int i = 0; i < distance; ++i) {
            const auto d1 = model.GetSystemSize() - distance + i;
            const auto d2 = i;
            GenerateMatrixElementsIntersite(edme, basis, d1, model.GetOperatorSz(), d2, model.GetOperatorSz(), model.GetJz(distance - 1), 1.0);
         }
      }
      
      //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
      for (int distance = 1; distance <= model.GetJxy().size(); ++distance) {
         for (int i = 0; i < distance; ++i) {
            const auto d1 = model.GetSystemSize() - distance + i;
            const auto d2 = i;
            GenerateMatrixElementsIntersite(edme, basis, d1, model.GetOperatorSp(), d2, model.GetOperatorSm(), 0.5*model.GetJxy(distance - 1), 1.0);
            GenerateMatrixElementsIntersite(edme, basis, d1, model.GetOperatorSm(), d2, model.GetOperatorSp(), 0.5*model.GetJxy(distance - 1), 1.0);
         }
      }
   }
}


}
}


#endif /* exact_diag_utility_hpp */
