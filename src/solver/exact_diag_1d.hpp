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
#include <sstream>
#include <ios>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace solver {

template<class ModelClass1D>
class ExactDiag1D {
   
   using RealType = typename ModelClass1D::ValueType;
   
   using CRS = sparse_matrix::CRS<RealType>;
   
   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   
   ModelClass1D model;
   
   sparse_matrix::DiagonalizationParameters diag_params;
   
   sparse_matrix::DiagonalizationInfomation diag_info;
   
   explicit ExactDiag1D(const ModelClass1D &model_input): model(model_input) {
      diag_params.flag_symmetric_matrix = true;
      diag_params.flag_output_converge_step = true;
      dim_ = model.GetDim();
   }
   
   void SetBasis() {
      if (model.GetFlagRecalcBasis()) {
         model.GenerateBasis(&basis_, &basis_inv);
      }
   }
   
   void Diagonalize() {
      SetBasis();
      CRS ham;
      GenerateHamiltonian(&ham);
      sparse_matrix::EigenvalueDecompositionLanczos(&gs_value_, &gs_vector_, ham, diag_params, &diag_info);
   }
   
   CRS GenerateHamiltonian() {
      SetBasis();
      CRS ham;
      GenerateHamiltonian(&ham);
      return ham;
   }
   
   const BraketVector &GetGSVector() const {
      return gs_vector_;
   }
   
   RealType GetGSValue() const {
      return gs_value_;
   }
   
   const std::vector<int64_t> &GetBasis() const {
      return basis_;
   }
   
   int64_t GetDim() const {
      return dim_;
   }
   
   RealType CalculateExpectationValue(const CRS &M, int site) {
      return gs_vector_*CalculateMatrixVectorProduct(M, site, gs_vector_, basis_, basis_);
   }
   
   void PrintBasis() const {
      if (model.GetFlagRecalcBasis()) {
         throw std::runtime_error("basis has not been calculated.");
      }
      const int dim_onsite = model.GetDimOnsite();
      for (std::size_t i = 0; i < basis_.size(); ++i) {
         const int64_t basis = basis_[i];
         std::cout << "basis[" << std::setw(3) << i << "]=" << std::setw(3) << basis << ":";
         for (int site = 0; site < model.GetSystemSize(); ++site) {
            int basis_onsite = CalculateLocalBasis(basis, site, dim_onsite);
            model.PrintBasisOnsite(basis_onsite, false);
         }
         std::cout << std::endl;
      }
   }

   
private:
   int64_t dim_;
   std::vector<int64_t> basis_;
   std::unordered_map<int64_t, int64_t> basis_inv;
   BraketVector gs_vector_;
   RealType gs_value_;
   bool flag_recalc_gs = true;
   
   void GenerateHamiltonian(CRS *ham) const {

      const int64_t dim_target = model.GetDim();
      int64_t num_total_elements = 0;
      
#ifdef _OPENMP
      const int num_threads = omp_get_num_threads();
      std::vector<ExactDiagMatrixElements<RealType>> components(num_threads);
      
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         components[thread_num].site_constant.resize(model.GetSystemSize());
         for (int site = 0; site < model.GetSystemSize(); ++site) {
            components[thread_num].site_constant[site] = CalculatePower(model.GetDimOnsite(), site);
         }
         components[thread_num].basis_onsite.resize(model.GetSystemSize());
      }
      
      std::vector<int64_t> num_row_element(dim_target + 1);
      
#pragma omp parallel for num_threads (num_threads)
      for (int row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixElements(&components[thread_num], basis_[row], model);
         for (const auto &a_basis: components[thread_num].basis_affected) {
            if (basis_inv.count(a_basis) > 0) {
               const int64_t inv = basis_inv.at(a_basis);
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
      for (int64_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      //Do not use openmp here
      for (int64_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->ResizeRow(dim_target);
      ham->ResizeColVal(num_total_elements);
      
#pragma omp parallel for num_threads (num_threads)
      for (int row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixElements(&components[thread_num], basis_[row], model);
         for (int64_t i = 0; i < components[thread_num].basis_affected.size(); ++i) {
            const int64_t  a_basis = components[thread_num].basis_affected[i];
            const RealType val     = components[thread_num].val[i];
            if (basis_inv.count(a_basis) > 0) {
               const int64_t inv = basis_inv.at(a_basis);
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
      
      std::vector<int64_t> num_row_element(dim_target + 1);
      
      for (int64_t row = 0; row < dim_target; ++row) {
         GenerateMatrixElements(&components, basis_[row], model);
         for (const auto &a_basis: components.basis_affected) {
            if (basis_inv.count(a_basis) > 0) {
               const int64_t inv = basis_inv.at(a_basis);
               if (inv <= row) {
                  num_row_element[row + 1]++;
               }
            }
         }
         components.val.clear();
         components.basis_affected.clear();
         components.inv_basis_affected.clear();
      }
      
      for (int64_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      for (int64_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->ResizeRow(dim_target);
      ham->ResizeColVal(num_total_elements);
      
      for (int row = 0; row < dim_target; ++row) {
         GenerateMatrixElements(&components, basis_[row], model);
         for (int64_t i = 0; i < components.basis_affected.size(); ++i) {
            const int64_t  a_basis = components.basis_affected[i];
            const RealType val     = components.val[i];
            if (basis_inv.count(a_basis) > 0) {
               const int64_t inv = basis_inv.at(a_basis);
               if (inv <= row) {
                  ham->Col(num_row_element[row]) = inv;
                  ham->Val(num_row_element[row]) = val;
                  num_row_element[row]++;
               }
            }
         }
         ham->Row(row + 1) = num_row_element[row];
         components.val.clear();
         components.basis_affected.clear();
         components.inv_basis_affected.clear();
      }
#endif
      
      ham->SetRowDim(dim_target);
      ham->SetColDim(dim_target);
      
      const bool flag_check_1 = (ham->Row(dim_target) != num_total_elements);
      const bool flag_check_2 = (ham->GetSizeCol() != num_total_elements);
      const bool flag_check_3 = (ham->GetSizeVal() != num_total_elements);

      if (flag_check_1 || flag_check_2 || flag_check_3) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__  << std::endl;
         throw std::runtime_error(ss.str());
      }
      ham->SortCol();
   }
   
   BraketVector CalculateMatrixVectorProduct(const CRS &M, int site,
                                             const BraketVector &ket_in,
                                             const std::vector<int64_t> &bases_in,
                                             const std::vector<int64_t> &bases_out) const {
      
      if (bases_in.size() != ket_in.GetDim()) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The sizes of bases and ket_in does not match each other" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      int64_t dim_onsite    = M.GetRowDim();
      //int64_t dim_in        = bases_in.size();
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
