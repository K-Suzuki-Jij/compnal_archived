//
//  exact_diag.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/07.
//

#ifndef COMPNAL_SOLVER_EXACT_DIAG_HPP_
#define COMPNAL_SOLVER_EXACT_DIAG_HPP_

#include "../model/all.hpp"
#include "../sparse_matrix/all.hpp"

namespace compnal {
namespace solver {

template<typename RealType>
struct ExactDiagMatrixComponents {
   std::vector<RealType> val;
   std::vector<std::size_t>  basis_affected;
   std::vector<int>          basis_onsite;
   std::vector<std::size_t>  site_constant;
   std::unordered_map<std::size_t, std::size_t> inv_basis_affected;
   double zero_precision = std::pow(10, -15);
};

template<class ModelClass>
class ExactDiag {
   
   using RealType = typename ModelClass::ValueType;
   
   using CRS = sparse_matrix::CRS<RealType>;
   
   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   
   ModelClass model;
   sparse_matrix::ParametersAll params;
   
   explicit ExactDiag(const ModelClass &model_input): model(model_input) {}
   ExactDiag(const ModelClass &model_input, const sparse_matrix::ParametersAll &params_input): model(model_input), params(params_input) {}
   
   void GenerateBasis() {
      if (model.GetFlagRecalcBasis()) {
         model.GenerateBasis(&basis_, &basis_inv_);
      }
   }
   
   
   
private:
   std::vector<std::size_t> basis_;
   std::unordered_map<std::size_t, std::size_t> basis_inv_;
   BraketVector gs_vector_;
   RealType gs_value_ = 0.0;
   
   int CalculateLocalBasis(std::size_t global_basis, const int site, const int dim_onsite) {
      for (int i = 0; i < site; ++i) {
         global_basis = global_basis/dim_onsite;
      }
      return static_cast<int>(global_basis%dim_onsite);
   }
   
   void GenerateMatrixComponentsOnsite(ExactDiagMatrixComponents<RealType> *edmc,
                                       const std::size_t basis,
                                       const int site,
                                       const CRS &matrix_onsite,
                                       const RealType coeef) {
      
      if (std::abs(coeef) <= edmc->zero_precision) {
         return;
      }
      
      const int         basis_onsite  = edmc->basis_onsite[site];
      const std::size_t site_constant = edmc->site_constant[site];
      
      for (std::size_t i = matrix_onsite.row[basis_onsite]; i < matrix_onsite.row[basis_onsite + 1]; ++i) {
         const std::size_t a_basis = basis + (matrix_onsite.col[i] - basis_onsite)*site_constant;
         if (edmc->inv_basis_affected.count(a_basis) == 0) {
            edmc->inv_basis_affected[a_basis] = edmc->basis_affected.size();
            edmc->val.push_back(coeef*matrix_onsite.val[i]);
            edmc->basis_affected.push_back(a_basis);
         }
         else {
            const RealType val = edmc->val[edmc->inv_basis_affected.at(a_basis)] + coeef*matrix_onsite.val[i];
            if (std::abs(val) > edmc->zero_precision) {
               edmc->val[edmc->inv_basis_affected.at(a_basis)] = val;
            }
         }
      }
   }
   
   void GenerateMatrixComponentsIntersite(ExactDiagMatrixComponents<RealType> *edmc,
                                          const std::size_t basis,
                                          const int site_1,
                                          const CRS &matrix_onsite_1,
                                          const int site_2,
                                          const CRS &matrix_onsite_2,
                                          const RealType coeef,
                                          const int fermion_sign = 1.0) {
      
      if (std::abs(coeef) <= edmc->zero_precision) {
         return;
      }
      
      const int         basis_onsite_1  = edmc->basis_onsite[site_1];
      const int         basis_onsite_2  = edmc->basis_onsite[site_2];
      const std::size_t site_constant_1 = edmc->site_constant[site_1];
      const std::size_t site_constant_2 = edmc->site_constant[site_2];
      
      for (std::size_t i1 = matrix_onsite_1.row[basis_onsite_1]; i1 < matrix_onsite_1.row[basis_onsite_1 + 1]; ++i1) {
         const RealType    val_1 = matrix_onsite_1.val[i1];
         const std::size_t col_1 = matrix_onsite_1.col[i1];
         for (std::size_t i2 = matrix_onsite_2.row[basis_onsite_2]; i2 < matrix_onsite_2.row[basis_onsite_2 + 1]; ++i2) {
            const std::size_t a_basis = basis + (col_1 - basis_onsite_1)*site_constant_1 + (matrix_onsite_2.col[i2] - basis_onsite_2)*site_constant_2;
            if (edmc->inv_basis_affected.count(a_basis) == 0) {
               edmc->inv_basis_affected[a_basis] = edmc->basis_affected.size();
               edmc->val.push_back(fermion_sign*coeef*val_1*matrix_onsite_2.val[i2]);
               edmc->basis_affected.push_back(a_basis);
            }
            else {
               const RealType val = edmc->val[edmc->inv_basis_affected.at(a_basis)] + fermion_sign*coeef*val_1*matrix_onsite_2.val[i2];
               if (std::abs(val) > edmc->zero_precision) {
                  edmc->val[edmc->inv_basis_affected.at(a_basis)] = val;
               }
            }
         }
      }
   }
   
   void GenerateHamiltonian(CRS *ham) const {
      const std::size_t dim_target = basis_.size();
      std::size_t num_total_elements = 0;
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<ExactDiagMatrixComponents<RealType>> components(num_threads);
      
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         components[thread_num].site_constant.resize(model.GetSystemSize());
         for (int site = 0; site < model.GetSystemSize(); ++site) {
            components[thread_num].site_constant[site] = std::pow(model.GetDimOnsite(), site);
         }
         components[thread_num].basis_onsite.resize(model.GetSystemSize());
      }
      
      std::vector<std::size_t> num_row_element(dim_target + 1);
      
#pragma omp parallel for
      for (std::size_t row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixElements(&components[thread_num], basis_[row], model);
         for (const auto &a_basis: components[thread_num].basis_affected) {
            if (basis_inv_.count(a_basis) > 0) {
               const std::size_t inv = basis_inv_.at(a_basis);
               if (inv <= row) {
                  num_row_element[row + 1]++;
               }
            }
         }
         components[thread_num].val.clear();
         components[thread_num].basis_affected.clear();
         components[thread_num].inv_basis_affected.clear();
      }
      
#pragma omp parallel for reduction(+:num_total_elements) num_threads (num_threads)
      for (std::size_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      //Do not use openmp here
      for (std::size_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->ResizeRow(dim_target);
      ham->ResizeColVal(num_total_elements);
      
#pragma omp parallel for num_threads (num_threads)
      for (int row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixElements(&components[thread_num], basis_[row], model);
         for (std::size_t i = 0; i < components[thread_num].basis_affected.size(); ++i) {
            const std::size_t  a_basis = components[thread_num].basis_affected[i];
            const RealType val     = components[thread_num].val[i];
            if (basis_inv_.count(a_basis) > 0) {
               const std::size_t inv = basis_inv_.at(a_basis);
               if (inv <= row) {
                  ham->Col(num_row_element[row]) = inv;
                  ham->Val(num_row_element[row]) = val;
                  num_row_element[row]++;
               }
            }
         }
         ham->Row(row + 1) = num_row_element[row];
         components[thread_num].val.clear();
         components[thread_num].basis_affected.clear();
         components[thread_num].inv_basis_affected.clear();
      }
#else
      ExactDiagMatrixElements<RealType> components;
      components.site_constant.resize(model.GetSystemSize());
      for (int site = 0; site < model.GetSystemSize(); ++site) {
         components.site_constant[site] = model::CalculatePower(model.GetDimOnsite(), site);
      }
      components.basis_onsite.resize(model.GetSystemSize());
      
      std::vector<std::size_t> num_row_element(dim_target + 1);
      
      for (std::size_t row = 0; row < dim_target; ++row) {
         GenerateMatrixElements(&components, basis_[row], model);
         for (const auto &a_basis: components.basis_affected) {
            if (basis_inv_.count(a_basis) > 0) {
               const std::size_t inv = basis_inv_.at(a_basis);
               if (inv <= row) {
                  num_row_element[row + 1]++;
               }
            }
         }
         components.val.clear();
         components.basis_affected.clear();
         components.inv_basis_affected.clear();
      }
      
      for (std::size_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      for (std::size_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->row.resize(dim_target + 1);
      ham->col.resize(num_total_elements);
      ham->val.resize(num_total_elements);
      
      for (int row = 0; row < dim_target; ++row) {
         GenerateMatrixElements(&components, basis_[row], model);
         for (std::size_t i = 0; i < components.basis_affected.size(); ++i) {
            const std::size_t  a_basis = components.basis_affected[i];
            const RealType val     = components.val[i];
            if (basis_inv_.count(a_basis) > 0) {
               const std::size_t inv = basis_inv_.at(a_basis);
               if (inv <= row) {
                  ham->col[num_row_element[row]] = inv;
                  ham->val[num_row_element[row]] = val;
                  num_row_element[row]++;
               }
            }
         }
         ham->row[row + 1] = num_row_element[row];
         components.val.clear();
         components.basis_affected.clear();
         components.inv_basis_affected.clear();
      }
#endif
      
      ham->row_dim = dim_target;
      ham->col_dim = dim_target;
      
      const bool flag_check_1 = (ham->row[dim_target] != num_total_elements);
      const bool flag_check_2 = (ham->col.size() != num_total_elements);
      const bool flag_check_3 = (ham->val.size() != num_total_elements);
      
      if (flag_check_1 || flag_check_2 || flag_check_3) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__  << std::endl;
         throw std::runtime_error(ss.str());
      }
      ham->SortCol();
   }
   
   void GenerateMatrixElements(ExactDiagMatrixComponents<RealType> *edmc, const std::size_t basis, const model::XXZ_1D<RealType> &model) {
      
      for (int site = 0; site < model.GetSystemSize(); ++site) {
         edmc->basis_onsite[site] = CalculateLocalBasis(basis, site, model.GetDimOnsite());
      }
      
      //Onsite elements
      for (int site = 0; site < model.GetSystemSize(); ++site) {
         GenerateMatrixElementsOnsite(edmc, basis, site, model.GetOperatorHam(), 1.0);
      }
      
      //Intersite elements SzSz
      for (int distance = 1; distance <= model.GetJz().size(); ++distance) {
         for (int site = 0; site < model.GetSystemSize() - distance; ++site) {
            GenerateMatrixElementsIntersite(edmc, basis, site, model.GetOperatorSz(), site + distance, model.GetOperatorSz(), model.GetJz(distance - 1), 1.0);
         }
      }
      
      //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
      for (int distance = 1; distance <= model.GetJxy().size(); ++distance) {
         for (int site = 0; site < model.GetSystemSize() - distance; ++site) {
            GenerateMatrixElementsIntersite(edmc, basis, site, model.GetOperatorSp(), site + distance, model.GetOperatorSm(), 0.5*model.GetJxy(distance - 1), 1.0);
            GenerateMatrixElementsIntersite(edmc, basis, site, model.GetOperatorSm(), site + distance, model.GetOperatorSp(), 0.5*model.GetJxy(distance - 1), 1.0);
         }
      }
      
      if (model.GetBoundaryCondition() == utility::BoundaryCondition::PBC) {
         //Intersite elements SzSz
         for (int distance = 1; distance <= model.GetJz().size(); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model.GetSystemSize() - distance + i;
               const auto d2 = i;
               GenerateMatrixElementsIntersite(edmc, basis, d1, model.GetOperatorSz(), d2, model.GetOperatorSz(), model.GetJz(distance - 1), 1.0);
            }
         }
         
         //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
         for (int distance = 1; distance <= model.GetJxy().size(); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model.GetSystemSize() - distance + i;
               const auto d2 = i;
               GenerateMatrixElementsIntersite(edmc, basis, d1, model.GetOperatorSp(), d2, model.GetOperatorSm(), 0.5*model.GetJxy(distance - 1), 1.0);
               GenerateMatrixElementsIntersite(edmc, basis, d1, model.GetOperatorSm(), d2, model.GetOperatorSp(), 0.5*model.GetJxy(distance - 1), 1.0);
            }
         }
      }
      
      //Fill zero in the diagonal elements for symmetric matrix vector product calculation.
      if (edmc->inv_basis_affected.count(basis) == 0) {
         edmc->inv_basis_affected[basis] = edmc->basis_affected.size();
         edmc->val.push_back(0.0);
         edmc->basis_affected.push_back(basis);
      }
      
   }
   
   
   
};

}
}


#endif /* COMPNAL_SOLVER_EXACT_DIAG_HPP_ */
