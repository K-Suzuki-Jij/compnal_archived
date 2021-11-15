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
   
   explicit ExactDiag(const ModelClass &model_input): model(model_input) {
      params.lanczos.flag_symmetric_crs = true;
      params.ii.cg.flag_symmetric_crs   = true;
   }
   
   ExactDiag(const ModelClass &model_input, const sparse_matrix::ParametersAll &params_input): model(model_input), params(params_input) {
      params.lanczos.flag_symmetric_crs = true;
      params.ii.cg.flag_symmetric_crs   = true;
   }
      
   void CalculateGroundState(const std::string &diag_method = "Lanczos") {
      model.GenerateBasis();
      CRS ham;
      GenerateHamiltonian(&ham);
      if (eigenvalues_.size() == 0) {
         eigenvalues_.emplace_back();
      }
      if (eigenvectors_.size() == 0) {
         eigenvectors_.emplace_back();
      }
      if (diag_method == "Lanczos") {
         sparse_matrix::EigenvalueDecompositionLanczos(&eigenvalues_[0], &eigenvectors_[0], ham, params.lanczos);
      }
      else if (diag_method == "LOBPCG") {
         sparse_matrix::EigenvalueDecompositionLOBPCG(&eigenvalues_[0], &eigenvectors_[0], ham, params.lanczos);
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid diag_method: " << diag_method << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      sparse_matrix::InverseIteration(&ham, &eigenvectors_[0], eigenvalues_[0], params.ii);
      
      model.SetCalculatedEigenvectorSet(0);
   }
   
   RealType CalculateExpectationValue(const CRS &m, const std::size_t site, const std::size_t level = 0) const {
      if (model.GetCalculatedEigenvectorSet().count(level) == 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "An eigenvector of the energy level: " << level << " has not been calculated" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const auto &basis     = model.GetTargetBasis();
      const auto &basis_inv = model.GetTargetBasisInv();
      
      const std::size_t dim = eigenvectors_.at(level).val.size();
      if (basis.size() != dim || basis_inv.size() != dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Unknown error detected" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const int dim_onsite = static_cast<int>(m.row_dim);
      const std::size_t site_constant = static_cast<std::size_t>(std::pow(dim_onsite, site));
      const BraketVector &eigenvector = eigenvectors_.at(level);
      RealType val = 0.0;
      
#pragma omp parallel for reduction (+: val)
      for (std::size_t i = 0; i < dim; ++i) {
         const std::size_t global_basis = basis[i];
         const int         local_basis = CalculateLocalBasis(global_basis, site, dim_onsite);
         RealType temp_val = 0.0;
         for (std::size_t j = m.row[local_basis]; j < m.row[local_basis + 1]; ++j){
            const std::size_t a_basis = global_basis - (local_basis - m.col[j])*site_constant;
            if (basis_inv.count(a_basis) != 0) {
               temp_val += eigenvector.val[basis_inv.at(a_basis)]*m.val[j];
            }
         }
         val += eigenvector.val[i]*temp_val;
      }
      return val;
   }
   
   RealType CalculateCorrelationFunction(const CRS &m_1, const std::size_t site_1, const CRS &m_2, const std::size_t site_2, const std::size_t target_level = 0) {
      if (model.GetCalculatedEigenvectorSet().count(target_level) == 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "An eigenvector of the energy level: " << target_level << " has not been calculated" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (m_1.row_dim != m_2.row_dim || m_1.row_dim != m_1.col_dim || m_2.row_dim != m_2.col_dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid input of the local operators" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const CRS m1_dagger = sparse_matrix::CalculateTransposedMatrix(m_1);
      
      const auto level_set = model.GenerateTargetSector(m1_dagger, m_2);
      const auto &basis_inv = model.GetTargetBasisInv();
      const int dim_onsite = static_cast<int>(m1_dagger.row_dim);
      const std::size_t site_constant_m1 = static_cast<std::size_t>(std::pow(dim_onsite, site_1));
      const std::size_t site_constant_m2 = static_cast<std::size_t>(std::pow(dim_onsite, site_2));
      const BraketVector &eigenvector = eigenvectors_.at(target_level);
      BraketVector vector_work_m1;
      BraketVector vector_work_m2;
      RealType val = 0.0;
      
      for (const auto &level: level_set) {
         if (model.isValidQNumber(level)) {
            model.GenerateBasis(level);
            const auto &basis = model.GetBasis(level);
            const std::size_t dim_target = basis.size();
            vector_work_m1.val.resize(dim_target);
            vector_work_m2.val.resize(dim_target);
   #pragma omp parallel for
            for (std::size_t i = 0; i < dim_target; ++i) {
               const std::size_t global_basis   = basis[i];
               const int         local_basis_m1 = CalculateLocalBasis(global_basis, static_cast<int>(site_1), dim_onsite);
               const int         local_basis_m2 = CalculateLocalBasis(global_basis, static_cast<int>(site_2), dim_onsite);
               RealType temp_val_m1 = 0.0;
               RealType temp_val_m2 = 0.0;
               for (std::size_t j = m1_dagger.row[local_basis_m1]; j < m1_dagger.row[local_basis_m1 + 1]; ++j){
                  const std::size_t a_basis = global_basis - (local_basis_m1 - m1_dagger.col[j])*site_constant_m1;
                  if (basis_inv.count(a_basis) != 0) {
                     temp_val_m1 += eigenvector.val[basis_inv.at(a_basis)]*m1_dagger.val[j];
                  }
               }
               vector_work_m1.val[i] = temp_val_m1;
               
               for (std::size_t j = m_2.row[local_basis_m2]; j < m_2.row[local_basis_m2 + 1]; ++j){
                  const std::size_t a_basis = global_basis - (local_basis_m2 - m_2.col[j])*site_constant_m2;
                  if (basis_inv.count(a_basis) != 0) {
                     temp_val_m2 += eigenvector.val[basis_inv.at(a_basis)]*m_2.val[j];
                  }
               }
               vector_work_m2.val[i] = temp_val_m2;
            }
            val += sparse_matrix::CalculateInnerProduct(vector_work_m1, vector_work_m2);
         }
      }
      return val;
   }
   
   inline const std::vector<BraketVector> &GetEigenvectors() const { return eigenvectors_; }
   inline const std::vector<RealType>     &GetEigenvalues()  const { return eigenvalues_; }
   
private:
   std::vector<BraketVector> eigenvectors_;
   std::vector<RealType>     eigenvalues_;
      
   int CalculateLocalBasis(std::size_t global_basis, const int site, const int dim_onsite) const {
      for (int i = 0; i < site; ++i) {
         global_basis = global_basis/dim_onsite;
      }
      return static_cast<int>(global_basis%dim_onsite);
   }
   
   void GenerateMatrixComponentsOnsite(ExactDiagMatrixComponents<RealType> *edmc,
                                       const std::size_t basis,
                                       const int site,
                                       const CRS &matrix_onsite,
                                       const RealType coeef) const {
      
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
                                          const int fermion_sign = 1.0) const {
      
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
               edmc->val[edmc->inv_basis_affected.at(a_basis)] += fermion_sign*coeef*val_1*matrix_onsite_2.val[i2];
            }
         }
      }
   }
   
   void GenerateHamiltonian(CRS *ham) const {
      const auto start = std::chrono::system_clock::now();
      std::cout << "Generating Hamiltonian..." << std::flush;
      
      const auto &basis     = model.GetTargetBasis();
      const auto &basis_inv = model.GetTargetBasisInv();
      
      const std::size_t dim_target = basis.size();
      std::size_t num_total_elements = 0;
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<ExactDiagMatrixComponents<RealType>> components(num_threads);
      
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         components[thread_num].site_constant.resize(model.GetSystemSize());
         for (int site = 0; site < model.GetSystemSize(); ++site) {
            components[thread_num].site_constant[site] = static_cast<std::size_t>(std::pow(model.GetDimOnsite(), site));
         }
         components[thread_num].basis_onsite.resize(model.GetSystemSize());
      }
      
      std::vector<std::size_t> num_row_element(dim_target + 1);
      
#pragma omp parallel for
      for (std::size_t row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixComponents(&components[thread_num], basis[row], model);
         for (const auto &a_basis: components[thread_num].basis_affected) {
            if (basis_inv.count(a_basis) > 0) {
               const std::size_t inv = basis_inv.at(a_basis);
               if (inv <= row) {
                  num_row_element[row + 1]++;
               }
            }
         }
         components[thread_num].val.clear();
         components[thread_num].basis_affected.clear();
         components[thread_num].inv_basis_affected.clear();
      }
      
#pragma omp parallel for reduction(+:num_total_elements)
      for (std::size_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      //Do not use openmp here
      for (std::size_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->row.resize(dim_target + 1);
      ham->col.resize(num_total_elements);
      ham->val.resize(num_total_elements);

#pragma omp parallel for
      for (std::size_t row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixComponents(&components[thread_num], basis[row], model);
         for (std::size_t i = 0; i < components[thread_num].basis_affected.size(); ++i) {
            const std::size_t  a_basis = components[thread_num].basis_affected[i];
            const RealType     val     = components[thread_num].val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::size_t inv = basis_inv.at(a_basis);
               if (inv <= row) {
                  ham->col[num_row_element[row]] = inv;
                  ham->val[num_row_element[row]] = val;
                  num_row_element[row]++;
               }
            }
         }
         ham->row[row + 1] = num_row_element[row];
         components[thread_num].val.clear();
         components[thread_num].basis_affected.clear();
         components[thread_num].inv_basis_affected.clear();
      }
#else
      ExactDiagMatrixComponents<RealType> components;
      components.site_constant.resize(model.GetSystemSize());
      for (int site = 0; site < model.GetSystemSize(); ++site) {
         components.site_constant[site] = static_cast<std::size_t>(std::pow(model.GetDimOnsite(), site));
      }
      components.basis_onsite.resize(model.GetSystemSize());
      
      std::vector<std::size_t> num_row_element(dim_target + 1);
      
      for (std::size_t row = 0; row < dim_target; ++row) {
         GenerateMatrixComponents(&components, basis[row], model);
         for (const auto &a_basis: components.basis_affected) {
            if (basis_inv.count(a_basis) > 0) {
               const std::size_t inv = basis_inv.at(a_basis);
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
         GenerateMatrixComponents(&components, basis[row], model);
         for (std::size_t i = 0; i < components.basis_affected.size(); ++i) {
            const std::size_t  a_basis = components.basis_affected[i];
            const RealType val     = components.val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::size_t inv = basis_inv.at(a_basis);
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
      
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating Hamiltonian:" << time_sec << "[sec]" << std::endl;
      
   }
   
   void GenerateMatrixComponents(ExactDiagMatrixComponents<RealType> *edmc, const std::size_t basis, const model::XXZ_1D<RealType> &model_input) const {
      
      for (int site = 0; site < model_input.GetSystemSize(); ++site) {
         edmc->basis_onsite[site] = CalculateLocalBasis(basis, site, model_input.GetDimOnsite());
      }
      
      //Onsite elements
      for (int site = 0; site < model_input.GetSystemSize(); ++site) {
         GenerateMatrixComponentsOnsite(edmc, basis, site, model_input.GetOnsiteOperatorHam(), 1.0);
      }
      
      //Intersite elements SzSz
      for (int distance = 1; distance <= static_cast<int>(model_input.GetJz().size()); ++distance) {
         for (int site = 0; site < model_input.GetSystemSize() - distance; ++site) {
            GenerateMatrixComponentsIntersite(edmc, basis, site, model_input.GetOnsiteOperatorSz(), site + distance, model_input.GetOnsiteOperatorSz(), model_input.GetJz(distance - 1));
         }
      }
      
      //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
      for (int distance = 1; distance <= static_cast<int>(model_input.GetJxy().size()); ++distance) {
         for (int site = 0; site < model_input.GetSystemSize() - distance; ++site) {
            GenerateMatrixComponentsIntersite(edmc, basis, site, model_input.GetOnsiteOperatorSp(), site + distance, model_input.GetOnsiteOperatorSm(), 0.5*model_input.GetJxy(distance - 1));
            GenerateMatrixComponentsIntersite(edmc, basis, site, model_input.GetOnsiteOperatorSm(), site + distance, model_input.GetOnsiteOperatorSp(), 0.5*model_input.GetJxy(distance - 1));
         }
      }
      
      if (model_input.GetBoundaryCondition() == utility::BoundaryCondition::PBC) {
         //Intersite elements SzSz
         for (int distance = 1; distance <= static_cast<int>(model_input.GetJz().size()); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model_input.GetSystemSize() - distance + i;
               const auto d2 = i;
               GenerateMatrixComponentsIntersite(edmc, basis, d1, model_input.GetOnsiteOperatorSz(), d2, model_input.GetOnsiteOperatorSz(), model_input.GetJz(distance - 1));
            }
         }
         
         //Intersite elements 0.5*(SpSm + SmSp) = SxSx + SySy
         for (int distance = 1; distance <= static_cast<int>(model_input.GetJxy().size()); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model_input.GetSystemSize() - distance + i;
               const auto d2 = i;
               GenerateMatrixComponentsIntersite(edmc, basis, d1, model_input.GetOnsiteOperatorSp(), d2, model_input.GetOnsiteOperatorSm(), 0.5*model_input.GetJxy(distance - 1));
               GenerateMatrixComponentsIntersite(edmc, basis, d1, model_input.GetOnsiteOperatorSm(), d2, model_input.GetOnsiteOperatorSp(), 0.5*model.GetJxy(distance - 1));
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
