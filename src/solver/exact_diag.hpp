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
   std::vector<std::int64_t>  basis_affected;
   std::vector<int>           basis_onsite;
   std::vector<std::int64_t>  site_constant;
   std::unordered_map<std::int64_t, std::int64_t> inv_basis_affected;
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
      if (model.GetCalculatedEigenvectorSet().count(0) != 0) {
         return;
      }
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
   
   void CalculateTargetState(const int target_sector, const std::string &diag_method = "Lanczos") {
      if (target_sector < 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid target_sector: " << target_sector << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (target_sector == 0) {
         CalculateGroundState(diag_method);
         return;
      }
      if (model.GetCalculatedEigenvectorSet().count(target_sector) != 0) {
         return;
      }
      
      model.GenerateBasis();
      CRS ham;
      GenerateHamiltonian(&ham);
      if (diag_method == "Lanczos") {
         for (int sector = 1; sector <= target_sector; ++sector) {
            if (model.GetCalculatedEigenvectorSet().count(sector) == 0) {
               if (eigenvectors_.size() != sector) {
                  std::stringstream ss;
                  ss << "Unknown Error in " << __func__ << std::endl;
                  throw std::runtime_error(ss.str());
               }
               sparse_matrix::BraketVector<RealType> temp_vector(ham.row_dim);
               RealType temp_value = 0.0;
               sparse_matrix::EigenvalueDecompositionLanczos(&temp_value, &temp_vector, ham, sector, eigenvectors_);
               eigenvalues_.push_back(temp_value);
               eigenvectors_.push_back(temp_vector);
               model.SetCalculatedEigenvectorSet(sector);
            }
         }
      }
      else if (diag_method == "LOBPCG") {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "LOBPCG is under construction: " << target_sector << std::endl;
         throw std::runtime_error(ss.str());
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid diag_method: " << target_sector << std::endl;
         throw std::runtime_error(ss.str());
      }
      
   }
   
   RealType CalculateExpectationValue(const CRS &m, const int site, const int level = 0) const {
      if (model.GetCalculatedEigenvectorSet().count(level) == 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "An eigenvector of the energy level: " << level << " has not been calculated" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const auto &basis     = model.GetTargetBasis();
      const auto &basis_inv = model.GetTargetBasisInv();
      
      const std::int64_t dim = eigenvectors_.at(level).val.size();
      if (static_cast<std::int64_t>(basis.size()) != dim || static_cast<std::int64_t>(basis_inv.size()) != dim) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Unknown error detected" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const int dim_onsite = static_cast<int>(m.row_dim);
      const std::int64_t site_constant = static_cast<std::int64_t>(std::pow(dim_onsite, site));
      const BraketVector &eigenvector = eigenvectors_.at(level);
      RealType val = 0.0;
      
#pragma omp parallel for reduction (+: val)
      for (std::int64_t i = 0; i < dim; ++i) {
         const std::int64_t global_basis = basis[i];
         const int          local_basis = CalculateLocalBasis(global_basis, site, dim_onsite);
         RealType temp_val = 0.0;
         for (std::int64_t j = m.row[local_basis]; j < m.row[local_basis + 1]; ++j){
            const std::int64_t a_basis = global_basis - (local_basis - m.col[j])*site_constant;
            if (basis_inv.count(a_basis) != 0) {
               temp_val += eigenvector.val[basis_inv.at(a_basis)]*m.val[j];
            }
         }
         val += eigenvector.val[i]*temp_val;
      }
      return val;
   }
   
   RealType CalculateCorrelationFunction(const CRS &m_1, const int site_1, const CRS &m_2, const int site_2, const int target_level = 0) {
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
      
      const auto target_sectors  = model.GenerateTargetSector(m1_dagger, m_2);
      const auto &basis_inv = model.GetTargetBasisInv();
      const int dim_onsite = static_cast<int>(m1_dagger.row_dim);
      const std::int64_t site_constant_m1 = static_cast<std::int64_t>(std::pow(dim_onsite, site_1));
      const std::int64_t site_constant_m2 = static_cast<std::int64_t>(std::pow(dim_onsite, site_2));
      const BraketVector &eigenvector = eigenvectors_.at(target_level);
      BraketVector vector_work_m1;
      BraketVector vector_work_m2;
      RealType val = 0.0;
      
      for (const auto &sector: target_sectors) {
         if (model.isValidQNumber(sector)) {
            model.GenerateBasis(sector);
            const auto &basis = model.GetBasis(sector);
            const std::int64_t dim_target = static_cast<std::int64_t>(basis.size());
            vector_work_m1.val.resize(dim_target);
            vector_work_m2.val.resize(dim_target);
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target; ++i) {
               const std::int64_t global_basis = basis[i];
               const int local_basis_m1 = CalculateLocalBasis(global_basis, site_1, dim_onsite);
               const int local_basis_m2 = CalculateLocalBasis(global_basis, site_2, dim_onsite);
               RealType temp_val_m1 = 0.0;
               RealType temp_val_m2 = 0.0;
               for (std::int64_t j = m1_dagger.row[local_basis_m1]; j < m1_dagger.row[local_basis_m1 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m1 - m1_dagger.col[j])*site_constant_m1;
                  if (basis_inv.count(a_basis) != 0) {
                     temp_val_m1 += eigenvector.val[basis_inv.at(a_basis)]*m1_dagger.val[j];
                  }
               }
               vector_work_m1.val[i] = temp_val_m1;
               
               for (std::int64_t j = m_2.row[local_basis_m2]; j < m_2.row[local_basis_m2 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m2 - m_2.col[j])*site_constant_m2;
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
   
   RealType CalculateCorrelationFunction(const CRS &m_1, const int site_1,
                                         const CRS &m_2, const int site_2,
                                         const CRS &m_3, const int site_3,
                                         const int target_level = 0) {
      
      if (model.GetCalculatedEigenvectorSet().count(target_level) == 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "An eigenvector of the energy level: " << target_level << " has not been calculated" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const bool c1 = (m_1.row_dim == m_1.col_dim);
      const bool c2 = (m_2.row_dim == m_2.col_dim);
      const bool c3 = (m_3.row_dim == m_3.col_dim);
      const bool c4 = (m_1.row_dim == m_2.row_dim && m_2.row_dim == m_3.row_dim);
      
      if (!(c1 && c2 && c3 && c4)) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid input of the local operators" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (site_1 == site_2) {
         return CalculateCorrelationFunction(sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_1, 1.0, m_2), site_1, m_3, site_3, target_level);
      }
      if (site_2 == site_3) {
         return CalculateCorrelationFunction(m_1, site_1, sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_2, 1.0, m_3), site_3, target_level);
      }
      
      const CRS  m1_dagger      = sparse_matrix::CalculateTransposedMatrix(m_1);
      const auto target_sectors = model.GenerateTargetSector(m1_dagger, m_2, m_3);
      
      const auto &basis_gs_sector_inv = model.GetTargetBasisInv();
      const int dim_onsite = static_cast<int>(m1_dagger.row_dim);
      const std::int64_t site_constant_m1 = static_cast<std::int64_t>(std::pow(dim_onsite, site_1));
      const std::int64_t site_constant_m2 = static_cast<std::int64_t>(std::pow(dim_onsite, site_2));
      const std::int64_t site_constant_m3 = static_cast<std::int64_t>(std::pow(dim_onsite, site_3));
      const BraketVector &eigenvector = eigenvectors_.at(target_level);
      BraketVector vector_work_a;
      BraketVector vector_work_b;
      BraketVector vector_work_c;
      RealType val = 0.0;
      
      for (const auto &it: target_sectors) {
         const auto sector_m1 = std::get<0>(it);
         const auto sector_m3 = std::get<1>(it);
         if (model.isValidQNumber(sector_m1) && model.isValidQNumber(sector_m3)) {
            model.GenerateBasis(sector_m1);
            model.GenerateBasis(sector_m3);
            const auto &basis_m1 = model.GetBasis(sector_m1);
            const auto &basis_m3 = model.GetBasis(sector_m3);
            const std::int64_t dim_target_m1 = static_cast<std::int64_t>(basis_m1.size());
            const std::int64_t dim_target_m3 = static_cast<std::int64_t>(basis_m3.size());
            
            //m3|gs>
            vector_work_a.val.resize(dim_target_m3);
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target_m3; ++i) {
               const std::int64_t global_basis = basis_m3[i];
               const int local_basis = CalculateLocalBasis(global_basis, site_3, dim_onsite);
               RealType temp_val = 0.0;
               for (std::int64_t j = m_3.row[local_basis]; j < m_3.row[local_basis + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis - m_3.col[j])*site_constant_m3;
                  if (basis_gs_sector_inv.count(a_basis) != 0) {
                     temp_val += eigenvector.val[basis_gs_sector_inv.at(a_basis)]*m_3.val[j];
                  }
               }
               vector_work_a.val[i] = temp_val;
            }
            
            //m2 * m3|gs> and m1_dag|gs>
            vector_work_b.val.resize(dim_target_m1);
            vector_work_c.val.resize(dim_target_m1);
            const auto &basis_m3_sector_inv = model.GetBasisInv(sector_m3);
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target_m1; ++i) {
               const std::int64_t global_basis = basis_m1[i];
               const int local_basis_m1 = CalculateLocalBasis(global_basis, site_1, dim_onsite);
               const int local_basis_m2 = CalculateLocalBasis(global_basis, site_2, dim_onsite);
               RealType temp_val_m1 = 0.0;
               RealType temp_val_m2 = 0.0;
               
               for (std::int64_t j = m_2.row[local_basis_m2]; j < m_2.row[local_basis_m2 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m2 - m_2.col[j])*site_constant_m2;
                  if (basis_m3_sector_inv.count(a_basis) != 0) {
                     temp_val_m2 += vector_work_a.val[basis_m3_sector_inv.at(a_basis)]*m_2.val[j];
                  }
               }
               vector_work_b.val[i] = temp_val_m2;
               
               for (std::int64_t j = m_1.row[local_basis_m1]; j < m_1.row[local_basis_m1 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m1 - m_1.col[j])*site_constant_m1;
                  if (basis_gs_sector_inv.count(a_basis) != 0) {
                     temp_val_m1 += eigenvector.val[basis_gs_sector_inv.at(a_basis)]*m_1.val[j];
                  }
               }
               vector_work_c.val[i] = temp_val_m1;
            }
            val += sparse_matrix::CalculateInnerProduct(vector_work_b, vector_work_c);
         }
      }
      return val;
   }
   
   RealType CalculateCorrelationFunction(const CRS &m_1, const int site_1,
                                         const CRS &m_2, const int site_2,
                                         const CRS &m_3, const int site_3,
                                         const CRS &m_4, const int site_4,
                                         const int target_level = 0) {
      
      if (model.GetCalculatedEigenvectorSet().count(target_level) == 0) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "An eigenvector of the energy level: " << target_level << " has not been calculated" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const bool c1 = (m_1.row_dim == m_1.col_dim);
      const bool c2 = (m_2.row_dim == m_2.col_dim);
      const bool c3 = (m_3.row_dim == m_3.col_dim);
      const bool c4 = (m_4.row_dim == m_4.col_dim);
      const bool c5 = (m_1.row_dim == m_2.row_dim && m_2.row_dim == m_3.row_dim && m_3.row_dim == m_4.row_dim);
      
      if (!(c1 && c2 && c3 && c4 && c5)) {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "Invalid input of the local operators" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (site_1 == site_2) {
         return CalculateCorrelationFunction(sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_1, 1.0, m_2), site_1, m_3, site_3, m_4, site_4);
      }
      if (site_2 == site_3) {
         return CalculateCorrelationFunction(m_1, site_1, sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_2, 1.0, m_3), site_2, m_4, site_4);
      }
      if (site_3 == site_4) {
         return CalculateCorrelationFunction(m_1, site_1, m_2, site_2, sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_3, 1.0, m_4), site_3);
      }
      
      const CRS  m1_dagger = sparse_matrix::CalculateTransposedMatrix(m_1);
      const CRS  m2_dagger = sparse_matrix::CalculateTransposedMatrix(m_2);
      const auto target_sectors = model.GenerateTargetSector(m1_dagger, m2_dagger, m_3, m_4);
      
      const auto &basis_gs_sector_inv = model.GetTargetBasisInv();
      const int dim_onsite = static_cast<int>(m1_dagger.row_dim);
      const std::int64_t site_constant_m1 = static_cast<std::int64_t>(std::pow(dim_onsite, site_1));
      const std::int64_t site_constant_m2 = static_cast<std::int64_t>(std::pow(dim_onsite, site_2));
      const std::int64_t site_constant_m3 = static_cast<std::int64_t>(std::pow(dim_onsite, site_3));
      const std::int64_t site_constant_m4 = static_cast<std::int64_t>(std::pow(dim_onsite, site_4));
      const BraketVector &eigenvector = eigenvectors_.at(target_level);
      BraketVector vector_work_m1;
      BraketVector vector_work_m2m1;
      BraketVector vector_work_m4;
      BraketVector vector_work_m3m4;
      RealType val = 0.0;
      
      for (const auto &it: target_sectors) {
         const auto sector_bra_change_1 = std::get<0>(it);
         const auto sector_bra_change_2 = std::get<1>(it);
         const auto sector_ket_change_1 = std::get<2>(it);
         if (model.isValidQNumber(sector_bra_change_1) &&
             model.isValidQNumber(sector_bra_change_2) &&
             model.isValidQNumber(sector_ket_change_1)) {
            model.GenerateBasis(sector_bra_change_1);
            model.GenerateBasis(sector_bra_change_2);
            model.GenerateBasis(sector_ket_change_1);
            const auto &basis_bra_1 = model.GetBasis(sector_bra_change_1);
            const auto &basis_bra_2 = model.GetBasis(sector_bra_change_2);
            const auto &basis_ket_1 = model.GetBasis(sector_ket_change_1);
            const auto &basis_ket_1_inv = model.GetBasisInv(sector_ket_change_1);
            const auto &basis_bra_1_inv = model.GetBasisInv(sector_bra_change_1);
            const std::int64_t dim_target_bra_1 = static_cast<std::int64_t>(basis_bra_1.size());
            const std::int64_t dim_target_bra_2 = static_cast<std::int64_t>(basis_bra_2.size());
            const std::int64_t dim_target_ket_1 = static_cast<std::int64_t>(basis_ket_1.size());
            
            //m4|gs>
            vector_work_m4.val.resize(dim_target_ket_1);
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target_ket_1; ++i) {
               const std::int64_t global_basis = basis_ket_1[i];
               const int local_basis = CalculateLocalBasis(global_basis, site_4, dim_onsite);
               RealType temp_val = 0.0;
               for (std::int64_t j = m_4.row[local_basis]; j < m_4.row[local_basis + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis - m_4.col[j])*site_constant_m4;
                  if (basis_gs_sector_inv.count(a_basis) != 0) {
                     temp_val += eigenvector.val[basis_gs_sector_inv.at(a_basis)]*m_4.val[j];
                  }
               }
               vector_work_m4.val[i] = temp_val;
            }
            
            //m1_dag|gs>
            vector_work_m1.val.resize(dim_target_bra_1);
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target_bra_1; ++i) {
               const std::int64_t global_basis = basis_bra_1[i];
               const int local_basis = CalculateLocalBasis(global_basis, site_1, dim_onsite);
               RealType temp_val = 0.0;
               for (std::int64_t j = m1_dagger.row[local_basis]; j < m1_dagger.row[local_basis + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis - m1_dagger.col[j])*site_constant_m1;
                  if (basis_gs_sector_inv.count(a_basis) != 0) {
                     temp_val += eigenvector.val[basis_gs_sector_inv.at(a_basis)]*m1_dagger.val[j];
                  }
               }
               vector_work_m1.val[i] = temp_val;
            }
            
            //m3 * m4|gs> and m2_dag * m1_dag|gs>
            vector_work_m3m4.val.resize(dim_target_bra_2);
            vector_work_m2m1.val.resize(dim_target_bra_2);
            
#pragma omp parallel for
            for (std::int64_t i = 0; i < dim_target_bra_2; ++i) {
               const std::int64_t global_basis = basis_bra_2[i];
               const int local_basis_m2 = CalculateLocalBasis(global_basis, site_2, dim_onsite);
               const int local_basis_m3 = CalculateLocalBasis(global_basis, site_3, dim_onsite);
               RealType temp_val = 0.0;
               
               for (std::int64_t j = m_3.row[local_basis_m3]; j < m_3.row[local_basis_m3 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m3 - m_3.col[j])*site_constant_m3;
                  if (basis_ket_1_inv.count(a_basis) != 0) {
                     temp_val += vector_work_m4.val[basis_ket_1_inv.at(a_basis)]*m_3.val[j];
                  }
               }
               vector_work_m3m4.val[i] = temp_val;
               
               temp_val = 0.0;
               for (std::int64_t j = m2_dagger.row[local_basis_m2]; j < m2_dagger.row[local_basis_m2 + 1]; ++j){
                  const std::int64_t a_basis = global_basis - (local_basis_m2 - m2_dagger.col[j])*site_constant_m2;
                  if (basis_bra_1_inv.count(a_basis) != 0) {
                     temp_val += vector_work_m1.val[basis_bra_1_inv.at(a_basis)]*m2_dagger.val[j];
                  }
               }
               vector_work_m2m1.val[i] = temp_val;
            }
            val += sparse_matrix::CalculateInnerProduct(vector_work_m2m1, vector_work_m3m4);
         }
      }
      return val;
   }
   
   inline const std::vector<BraketVector> &GetEigenvectors() const { return eigenvectors_; }
   inline const std::vector<RealType>     &GetEigenvalues()  const { return eigenvalues_; }
   
private:
   std::vector<BraketVector> eigenvectors_;
   std::vector<RealType>     eigenvalues_;
   
   int CalculateLocalBasis(std::int64_t global_basis, const int site, const int dim_onsite) const {
      for (int i = 0; i < site; ++i) {
         global_basis = global_basis/dim_onsite;
      }
      return static_cast<int>(global_basis%dim_onsite);
   }
   
   void GenerateMatrixComponentsOnsite(ExactDiagMatrixComponents<RealType> *edmc,
                                       const std::int64_t basis,
                                       const int site,
                                       const CRS &matrix_onsite,
                                       const RealType coeef) const {
      
      if (std::abs(coeef) <= edmc->zero_precision) {
         return;
      }
      
      const int          basis_onsite  = edmc->basis_onsite[site];
      const std::int64_t site_constant = edmc->site_constant[site];
      
      for (std::int64_t i = matrix_onsite.row[basis_onsite]; i < matrix_onsite.row[basis_onsite + 1]; ++i) {
         const std::int64_t a_basis = basis + (matrix_onsite.col[i] - basis_onsite)*site_constant;
         if (edmc->inv_basis_affected.count(a_basis) == 0) {
            edmc->inv_basis_affected[a_basis] = edmc->basis_affected.size();
            edmc->val.push_back(coeef*matrix_onsite.val[i]);
            edmc->basis_affected.push_back(a_basis);
         }
         else {
            edmc->val[edmc->inv_basis_affected.at(a_basis)] += coeef*matrix_onsite.val[i];
         }
      }
   }
   
   void GenerateMatrixComponentsIntersite(ExactDiagMatrixComponents<RealType> *edmc,
                                          const std::int64_t basis,
                                          const int site_1,
                                          const CRS &matrix_onsite_1,
                                          const int site_2,
                                          const CRS &matrix_onsite_2,
                                          const RealType coeef,
                                          const int fermion_sign = 1.0) const {
      
      if (std::abs(coeef) <= edmc->zero_precision) {
         return;
      }
      
      const int          basis_onsite_1  = edmc->basis_onsite[site_1];
      const int          basis_onsite_2  = edmc->basis_onsite[site_2];
      const std::int64_t site_constant_1 = edmc->site_constant[site_1];
      const std::int64_t site_constant_2 = edmc->site_constant[site_2];
      
      for (std::int64_t i1 = matrix_onsite_1.row[basis_onsite_1]; i1 < matrix_onsite_1.row[basis_onsite_1 + 1]; ++i1) {
         const RealType     val_1 = matrix_onsite_1.val[i1];
         const std::int64_t col_1 = matrix_onsite_1.col[i1];
         for (std::int64_t i2 = matrix_onsite_2.row[basis_onsite_2]; i2 < matrix_onsite_2.row[basis_onsite_2 + 1]; ++i2) {
            const std::int64_t a_basis = basis + (col_1 - basis_onsite_1)*site_constant_1 + (matrix_onsite_2.col[i2] - basis_onsite_2)*site_constant_2;
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
      
      const std::int64_t dim_target = basis.size();
      std::int64_t num_total_elements = 0;
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<ExactDiagMatrixComponents<RealType>> components(num_threads);
      
      for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
         components[thread_num].site_constant.resize(model.GetSystemSize());
         for (int site = 0; site < model.GetSystemSize(); ++site) {
            components[thread_num].site_constant[site] = static_cast<std::int64_t>(std::pow(model.GetDimOnsite(), site));
         }
         components[thread_num].basis_onsite.resize(model.GetSystemSize());
      }
      
      std::vector<std::int64_t> num_row_element(dim_target + 1);
      
#pragma omp parallel for
      for (std::int64_t row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixComponents(&components[thread_num], basis[row], model);
         const std::size_t size = components[thread_num].basis_affected.size();
         for (std::size_t i = 0; i < size; ++i) {
            const std::int64_t a_basis = components[thread_num].basis_affected[i];
            const RealType     val     = components[thread_num].val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::int64_t inv = basis_inv.at(a_basis);
               if ((inv < row && std::abs(val) > components[thread_num].zero_precision) || inv == row) {
                  num_row_element[row + 1]++;
               }
            }
            else if (basis_inv.count(a_basis) == 0 && std::abs(val) > components[thread_num].zero_precision) {
               throw std::runtime_error("Matrix elements are not in the target space");
            }
         }
         components[thread_num].val.clear();
         components[thread_num].basis_affected.clear();
         components[thread_num].inv_basis_affected.clear();
      }
      
#pragma omp parallel for reduction(+:num_total_elements)
      for (std::int64_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      //Do not use openmp here
      for (std::int64_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->row.resize(dim_target + 1);
      ham->col.resize(num_total_elements);
      ham->val.resize(num_total_elements);
      
#pragma omp parallel for
      for (std::int64_t row = 0; row < dim_target; ++row) {
         const int thread_num = omp_get_thread_num();
         GenerateMatrixComponents(&components[thread_num], basis[row], model);
         const std::size_t size = components[thread_num].basis_affected.size();
         for (std::size_t i = 0; i < size; ++i) {
            const std::int64_t a_basis = components[thread_num].basis_affected[i];
            const RealType     val     = components[thread_num].val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::int64_t inv = basis_inv.at(a_basis);
               if ((inv < row && std::abs(val) > components[thread_num].zero_precision) || inv == row) {
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
         components.site_constant[site] = static_cast<std::int64_t>(std::pow(model.GetDimOnsite(), site));
      }
      components.basis_onsite.resize(model.GetSystemSize());
      
      std::vector<std::int64_t> num_row_element(dim_target + 1);
      
      for (std::int64_t row = 0; row < dim_target; ++row) {
         GenerateMatrixComponents(&components, basis[row], model);
         const std::size_t size = components.basis_affected.size();
         for (std::size_t i = 0; i < size; ++i) {
            const std::int64_t a_basis = components.basis_affected[i];
            const RealType     val     = components.val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::int64_t inv = basis_inv.at(a_basis);
               if ((inv <= row && std::abs(val) > components.zero_precision) || inv == row) {
                  num_row_element[row + 1]++;
               }
            }
            else if (basis_inv.count(a_basis) == 0 && std::abs(val) > components.zero_precision) {
               throw std::runtime_error("Matrix elements are not in the target space");
            }
         }
         components.val.clear();
         components.basis_affected.clear();
         components.inv_basis_affected.clear();
      }
      
      for (std::int64_t row = 0; row <= dim_target; ++row) {
         num_total_elements += num_row_element[row];
      }
      
      for (std::int64_t row = 0; row < dim_target; ++row) {
         num_row_element[row + 1] += num_row_element[row];
      }
      
      ham->row.resize(dim_target + 1);
      ham->col.resize(num_total_elements);
      ham->val.resize(num_total_elements);
      
      for (int row = 0; row < dim_target; ++row) {
         GenerateMatrixComponents(&components, basis[row], model);
         const std::size_t size = components.basis_affected.size();
         for (std::size_t i = 0; i < size; ++i) {
            const std::int64_t a_basis = components.basis_affected[i];
            const RealType     val     = components.val[i];
            if (basis_inv.count(a_basis) > 0) {
               const std::int64_t inv = basis_inv.at(a_basis);
               if ((inv <= row && std::abs(val) > components.zero_precision) || inv == row) {
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
      const bool flag_check_2 = (static_cast<std::int64_t>(ham->col.size()) != num_total_elements);
      const bool flag_check_3 = (static_cast<std::int64_t>(ham->val.size()) != num_total_elements);
      
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
   
   void GenerateMatrixComponents(ExactDiagMatrixComponents<RealType> *edmc, const std::int64_t basis, const model::XXZ_1D<RealType> &model_input) const {
      
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
   
   void GenerateMatrixComponents(ExactDiagMatrixComponents<RealType> *edmc, const std::int64_t basis, const model::Hubbard_1D<RealType> &model_input) const {
      
      const auto &nc            = model_input.GetOnsiteOperatorNC();
      const auto &c_up          = model_input.GetOnsiteOperatorCUp();
      const auto &c_up_dagger   = model_input.GetOnsiteOperatorCUpDagger();
      const auto &c_down        = model_input.GetOnsiteOperatorCDown();
      const auto &c_down_dagger = model_input.GetOnsiteOperatorCDownDagger();
      
      for (int site = 0; site < model_input.GetSystemSize(); ++site) {
         edmc->basis_onsite[site] = CalculateLocalBasis(basis, site, model_input.GetDimOnsite());
      }
      
      //Onsite elements
      for (int site = 0; site < model_input.GetSystemSize(); ++site) {
         GenerateMatrixComponentsOnsite(edmc, basis, site, model_input.GetOnsiteOperatorHam(), 1.0);
      }
      
      //Intersite Coulomb
      for (int distance = 1; distance <= static_cast<int>(model_input.GetIntersiteCoulomb().size()); ++distance) {
         for (int site = 0; site < model_input.GetSystemSize() - distance; ++site) {
            GenerateMatrixComponentsIntersite(edmc, basis, site, nc, site + distance, nc, model_input.GetIntersiteCoulomb(distance - 1));
         }
      }
      
      //Intersite Hopping
      for (int distance = 1; distance <= static_cast<int>(model_input.GetHopping().size()); ++distance) {
         for (int site = 0; site < model_input.GetSystemSize() - distance; ++site) {
            const auto t = model_input.GetHopping(distance - 1);
            int n_electron = 0;
            for (int s = site; s < site + distance; ++s) {
               n_electron += model_input.GetNumElectrons(edmc->basis_onsite[s]);
            }
            GenerateMatrixComponentsIntersite(edmc, basis, site, c_up_dagger  , site + distance, c_up         , +1.0*t, 2*(n_electron%2) - 1);
            GenerateMatrixComponentsIntersite(edmc, basis, site, c_up         , site + distance, c_up_dagger  , -1.0*t, 2*(n_electron%2) - 1);
            GenerateMatrixComponentsIntersite(edmc, basis, site, c_down_dagger, site + distance, c_down       , +1.0*t, 2*(n_electron%2) - 1);
            GenerateMatrixComponentsIntersite(edmc, basis, site, c_down       , site + distance, c_down_dagger, -1.0*t, 2*(n_electron%2) - 1);
         }
      }
      
      if (model_input.GetBoundaryCondition() == utility::BoundaryCondition::PBC) {
         //Intersite Coulomb
         for (int distance = 1; distance <= static_cast<int>(model_input.GetIntersiteCoulomb().size()); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model_input.GetSystemSize() - distance + i;
               const auto d2 = i;
               GenerateMatrixComponentsIntersite(edmc, basis, d1, nc, d2, nc, model_input.GetIntersiteCoulomb(distance - 1));
            }
         }
         
         //Intersite Hopping
         for (int distance = 1; distance <= static_cast<int>(model_input.GetHopping().size()); ++distance) {
            for (int i = 0; i < distance; ++i) {
               const auto d1 = model_input.GetSystemSize() - distance + i;
               const auto d2 = i;
               const auto t = model_input.GetHopping(distance - 1);
               int n_electron = 0;
               for (int s = d2; s < d2 + d1; ++s) {
                  n_electron += model_input.GetNumElectrons(edmc->basis_onsite[s]);
               }
               GenerateMatrixComponentsIntersite(edmc, basis, d1, c_up_dagger  , d2, c_up         , +1.0*t, 2*(n_electron%2) - 1);
               GenerateMatrixComponentsIntersite(edmc, basis, d1, c_up         , d2, c_up_dagger  , -1.0*t, 2*(n_electron%2) - 1);
               GenerateMatrixComponentsIntersite(edmc, basis, d1, c_down_dagger, d2, c_down       , +1.0*t, 2*(n_electron%2) - 1);
               GenerateMatrixComponentsIntersite(edmc, basis, d1, c_down       , d2, c_down_dagger, -1.0*t, 2*(n_electron%2) - 1);
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
   
   template<class BaseClass>
   void GenerateMatrixComponents(ExactDiagMatrixComponents<RealType> *edmc, const std::int64_t basis, const model::GeneralModel_1D<BaseClass> &model_input) const {
      
      for (int site = 0; site < model_input.GetSystemSize(); ++site) {
         edmc->basis_onsite[site] = CalculateLocalBasis(basis, site, model_input.GetDimOnsite());
      }
      
      //Onsite elements
      for (const auto &it: model_input.GetOnsiteOperatorList()) {
         GenerateMatrixComponentsOnsite(edmc, basis, it.first, it.second, 1.0);
      }
      
      //Intersite elements
      for (const auto &it1: model.GetIntersiteOperatorList()) {
         const int site_1 = it1.first;
         for (const auto &it2: it1.second) {
            const int site_2 = it2.first;
            for (const auto &it3: it2.second) {
               const auto &m_1 = it3.first;
               const auto &m_2 = it3.second;
               GenerateMatrixComponentsIntersite(edmc, basis, site_1, m_1, site_2, m_2, 1.0);
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
