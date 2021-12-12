//
//  base_u1_spin_electron_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/11.
//

#ifndef COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_1D_HPP_
#define COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"
#include "base_u1_spin_1d.hpp"
#include "base_u1_electron_1d.hpp"

#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

template<typename RealType>
class BaseU1SpinElectron_1D {
   
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   using ValueType = RealType;
   
   BaseU1SpinElectron_1D() {
      SetOnsiteOperator();
   }
   
   explicit BaseU1SpinElectron_1D(const int system_size): BaseU1SpinElectron_1D() {
      SetSystemSize(system_size);
   }
   
   BaseU1SpinElectron_1D(const int system_size, const double magnitude_lspin): BaseU1SpinElectron_1D(system_size) {
      SetMagnitudeLSpin(magnitude_lspin);
   }
   
   BaseU1SpinElectron_1D(const int system_size, const int total_electron): BaseU1SpinElectron_1D(system_size) {
      SetTotalElectron(total_electron);
   }
   
   BaseU1SpinElectron_1D(const int system_size,
                         const double magnitude_lspin,
                         const int total_electron): BaseU1SpinElectron_1D(system_size, magnitude_lspin) {
      SetTotalElectron(total_electron);
   }
   
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "system_size must be a non-negative integer" << std::endl;
         ss << "system_size=" << system_size << "is not allowed" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (system_size_ != system_size) {
         system_size_ = system_size;
         bases_.clear();
         bases_inv_.clear();
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetTotalSz(const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetTotalElectron(const int total_electron) {
      if (total_electron_ != total_electron) {
         total_electron_ = total_electron;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetMagnitudeLSpin(const double magnitude_lspin) {
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      if (magnitude_2lspin <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Please set magnitude_2lspin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_2lspin_ != magnitude_2lspin) {
         magnitude_2lspin_ = magnitude_2lspin;
         dim_onsite_lspin_ = magnitude_2lspin + 1;
         dim_onsite_       = dim_onsite_lspin_*dim_onsite_electron_;
         SetOnsiteOperator();
         bases_.clear();
         bases_inv_.clear();
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetCalculatedEigenvectorSet(const int level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   bool isValidQNumber() const {
      return isValidQNumber(system_size_, 0.5*magnitude_2lspin_, total_electron_, 0.5*total_2sz_);
   }
   
   bool isValidQNumber(const int total_electron, const double total_sz) const {
      return isValidQNumber(system_size_, 0.5*magnitude_2lspin_, total_electron, total_sz);
   }
   
   bool isValidQNumber(const std::pair<int, double> &quantum_number) const {
      return isValidQNumber(system_size_, 0.5*magnitude_2lspin_, quantum_number.first, quantum_number.second);
   }
   
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, total_electron_, 0.5*total_2sz_);
   }
   
   std::int64_t CalculateTargetDim(const int total_electron, const double total_sz) const {
      return CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, total_electron, total_sz);
   }
   
   std::int64_t CalculateTargetDim(const std::pair<int, double> &quantum_number) const {
      return CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, quantum_number.first, quantum_number.second);
   }
   
   void GenerateBasis() {
      GenerateBasis(total_electron_, 0.5*total_2sz_);
   }
   
   void GenerateBasis(const std::pair<int, double> &quantum_number) {
      GenerateBasis(quantum_number.first, quantum_number.second);
   }
   
   void GenerateBasis(const int total_electron, const double total_sz) {
      if (!isValidQNumber(total_electron, total_sz)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const auto start     = std::chrono::system_clock::now();
      const int  total_2sz = utility::DoubleTheNumber(total_sz);
      
      if (bases_.count({total_electron, total_2sz}) != 0) {
         return;
      }
      
      std::cout << "Generating Basis..." << std::flush;
      const int max_n_up_down = static_cast<int>(total_electron/2);
      const std::int64_t dim_target_global = CalculateTargetDim(total_electron, total_sz);
      
      std::vector<std::int64_t>().swap(bases_[{total_electron, total_2sz}]);
      auto &global_basis_ref = bases_.at({total_electron, total_2sz});
      
      std::vector<std::int64_t> site_constant_global(system_size_);
      std::vector<std::int64_t> site_constant_electron(system_size_);
      std::vector<std::int64_t> site_constant_lspin(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant_global[site]   = static_cast<std::int64_t>(std::pow(dim_onsite_         , site));
         site_constant_electron[site] = static_cast<std::int64_t>(std::pow(dim_onsite_electron_, site));
         site_constant_lspin[site]    = static_cast<std::int64_t>(std::pow(dim_onsite_lspin_   , site));
      }

      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size_);
      std::vector<int> q_number_spin_vec;
      std::vector<int> q_number_n_vac_vec;
      std::vector<int> q_number_n_up_vec;
      std::vector<int> q_number_n_down_vec;
      std::vector<int> q_number_n_up_down_vec;
      std::vector<std::int64_t> bias_basis;
      bias_basis.push_back(0);
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         for (int n_up = 0; n_up <= total_electron - 2*max_n_up_down; ++n_up) {
            const int n_down = total_electron - 2*max_n_up_down - n_up;
            const int n_vac  = system_size_ - n_up - n_down - n_up_down;
            if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
               const double total_sz_electron = 0.5*(n_up - n_down);
               if (BaseU1Spin_1D<RealType>::isValidQNumber(system_size_, 0.5*magnitude_2lspin_, total_sz - total_sz_electron)) {
                  q_number_spin_vec.push_back(total_2sz - (n_up - n_down));
                  q_number_n_vac_vec.push_back(n_vac);
                  q_number_n_up_vec.push_back(n_up );
                  q_number_n_down_vec.push_back(n_down);
                  q_number_n_up_down_vec.push_back(n_up_down);
                  const std::int64_t dim_electron = binom[system_size_][n_up]*binom[system_size_ - n_up][n_down]*binom[system_size_ - n_up - n_down][n_up_down];
                  const std::int64_t dim_lspin    = BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, total_sz - total_sz_electron);
                  bias_basis.push_back(dim_electron*dim_lspin);
               }
            }
         }
      }
      
      for (std::size_t i = 0; i < bias_basis.size(); ++i) {
         bias_basis[i+1] += bias_basis[i];
      }
      
      if (static_cast<std::int64_t>(bias_basis[bias_basis.size()]) != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }

      //Generate spin bases
      std::vector<int> temp_q_number_spin_vec = q_number_spin_vec;
      std::sort(temp_q_number_spin_vec.begin(), temp_q_number_spin_vec.end());
      temp_q_number_spin_vec.erase(std::unique(temp_q_number_spin_vec.begin(), temp_q_number_spin_vec.end()), temp_q_number_spin_vec.end());
      std::unordered_map<int, std::vector<std::int64_t>> spin_bases;
#pragma omp parallel for
      for (std::int64_t i = 0; i < static_cast<std::int64_t>(temp_q_number_spin_vec.size()); ++i) {
         const int total_2_sz_lspin = temp_q_number_spin_vec[i];
         const int shifted_2sz      = (system_size_*magnitude_2lspin_ - total_2_sz_lspin)/2;
         const std::int64_t dim_target_lspin = BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size_, magnitude_2lspin_/2, total_2_sz_lspin);
         std::vector<std::vector<int>> partition_integers;
         utility::GenerateIntegerPartition(&partition_integers, shifted_2sz, magnitude_2lspin_);
         auto &spin_basis = spin_bases[total_2_sz_lspin];
         spin_basis.reserve(dim_target_lspin);
         for (auto &&integer_list: partition_integers) {
            const bool condition1 = (0 < integer_list.size()) && (static_cast<int>(integer_list.size()) <= system_size_);
            const bool condition2 = (integer_list.size() == 0) && (shifted_2sz  == 0);
            if (condition1 || condition2) {
               for (int j = static_cast<int>(integer_list.size()); j < system_size_; ++j) {
                  integer_list.push_back(0);
               }
               std::sort(integer_list.begin(), integer_list.end());
               do {
                  std::int64_t basis_global = 0;
                  for (std::size_t j = 0; j < integer_list.size(); ++j) {
                     basis_global += integer_list[j]*site_constant_global[j];
                  }
                  spin_basis.push_back(basis_global);
               } while (std::next_permutation(integer_list.begin(), integer_list.end()));
            }
         }
         if (static_cast<std::int64_t>(spin_basis.size()) != dim_target_lspin) {
            std::stringstream ss;
            ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
            throw std::runtime_error(ss.str());
         }
         std::sort(spin_basis.begin(), spin_basis.end());
      }
      
      //Generate global bases
      const std::int64_t loop_size = static_cast<std::int64_t>(q_number_spin_vec.size());
      global_basis_ref.resize(dim_target_global);
#pragma omp parallel for
      for (std::int64_t i = 0; i < loop_size; ++i) {
         const int n_vac     = q_number_n_vac_vec[i];
         const int n_up      = q_number_n_up_vec[i];
         const int n_down    = q_number_n_down_vec[i];
         const int n_up_down = q_number_n_up_down_vec[i];
         const int total_2sz_lspin = q_number_spin_vec[i];
         std::int64_t count = bias_basis[i];
         std::vector<int> basis_list_electron(system_size_);
         for (int s = 0; s < n_vac; ++s) {
            basis_list_electron[s] = 0;
         }
         for (int s = 0; s < n_up; ++s) {
            basis_list_electron[s + n_vac] = 1;
         }
         for (int s = 0; s < n_down; ++s) {
            basis_list_electron[s + n_vac + n_up] = 2;
         }
         for (int s = 0; s < n_up_down; ++s) {
            basis_list_electron[s + n_vac + n_up + n_down] = 3;
         }
         
         do {
            std::int64_t basis_global_electron = 0;
            for (std::size_t j = 0; j < basis_list_electron.size(); ++j) {
               basis_global_electron += basis_list_electron[j]*site_constant_global[j]*dim_onsite_lspin_;
            }
            for (const auto &basis_global_lspin: spin_bases.at(total_2sz_lspin)) {
               global_basis_ref[count++] = basis_global_electron + basis_global_lspin;
               
            }
         } while (std::next_permutation(basis_list_electron.begin(), basis_list_electron.end()));
      }
      
      if (static_cast<std::int64_t>(bases_.at({total_electron, total_2sz}).size()) != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(global_basis_ref.begin(), global_basis_ref.end());
      
      bases_inv_[{total_electron, total_2sz}].clear();
      
      auto &basis_inv_ref = bases_inv_.at({total_electron, total_2sz});
      for (std::int64_t i = 0; i < dim_target_global; ++i) {
         basis_inv_ref[global_basis_ref[i]] = i;
      }
      
      if (basis_inv_ref.size() != global_basis_ref.size()) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         ss << "The same basis has been detected" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;

   }
   
   static bool isValidQNumber(const int system_size, const double magnitude_lspin, const int total_electron, const double total_sz) {
      const int total_2sz        = utility::DoubleTheNumber(total_sz);
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const bool c1 = (0 <= total_electron && total_electron <= 2*system_size);
      const bool c2 = ((total_electron + system_size*magnitude_2lspin - total_2sz)%2 == 0);
      const bool c3 = (-(total_electron + system_size*magnitude_2lspin) <= total_2sz);
      const bool c4 = (total_2sz <= total_electron + system_size*magnitude_2lspin);
      if (c1 && c2 && c3 && c4) {
         return true;
      }
      else {
         return false;
      }
   }
   
   static std::int64_t CalculateTargetDim(const int system_size, const double magnitude_lspin, const int total_electron, const double total_sz) {
      if (!isValidQNumber(system_size, magnitude_lspin, total_electron, total_sz)) {
         return 0;
      }
      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size);
      const int max_n_up_down = static_cast<int>(total_electron/2);
      std::int64_t dim = 0;
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         for (int n_up = 0; n_up <= total_electron - 2*max_n_up_down; ++n_up) {
            const int n_down = total_electron - 2*max_n_up_down - n_up;
            const int n_vac  = system_size - n_up - n_down - n_up_down;
            if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
               const double total_sz_electron  = 0.5*(n_up - n_down);
               const std::int64_t dim_electron = binom[system_size][n_up]*binom[system_size - n_up][n_down]*binom[system_size - n_up - n_down][n_up_down];
               const std::int64_t dim_lspin    = BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size, magnitude_lspin, total_sz - total_sz_electron);
               dim += dim_electron*dim_lspin;
            }
         }
      }
      return dim;
   }
   
   std::vector<std::pair<int, double>> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      // TODO: Check input matrics
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m2;
      for (std::int64_t i = 0; i < m_1.row_dim; ++i) {
         for (std::int64_t j = m_1.row[i]; j < m_1.row[i + 1]; ++j) {
            if (m_1.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(i, m_1.col[j]));
            }
         }
      }
      for (std::int64_t i = 0; i < m_2.row_dim; ++i) {
         for (std::int64_t j = m_2.row[i]; j < m_2.row[i + 1]; ++j) {
            if (m_2.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(i, m_2.col[j]));
            }
         }
      }
      std::vector<std::pair<int, double>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            if (del_sec_m1 == del_sec_m2) {
               target_sector_set.push_back({del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_});
            }
         }
      }
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   std::vector<std::pair<std::pair<int, double>, std::pair<int, double>>> GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_ket, const CRS &m_3_ket) const {
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m2;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m3;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(i, m_1_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_ket.row_dim; ++i) {
         for (std::int64_t j = m_2_ket.row[i]; j < m_2_ket.row[i + 1]; ++j) {
            if (m_2_ket.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(i, m_2_ket.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateQuntumNumberDifference(i, m_3_ket.col[j]));
            }
         }
      }
      
      std::vector<std::pair<std::pair<int, double>, std::pair<int, double>>> target_sector_set;
      
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               const std::pair<int, double> del_sec_m2_m3 = {del_sec_m2.first + del_sec_m3.first, del_sec_m2.second + del_sec_m3.second};
               if (del_sec_m1 == del_sec_m2_m3) {
                  target_sector_set.push_back({
                     {del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_},
                     {del_sec_m3.first + total_electron_, del_sec_m3.second + 0.5*total_2sz_}
                  });
               }
            }
         }
      }
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   std::vector<std::tuple<std::pair<int, double>, std::pair<int, double>, std::pair<int, double>>>
   GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_bra, const CRS &m_3_ket, const CRS &m_4_ket) const {
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m2;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m3;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m4;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(i, m_1_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_bra.row_dim; ++i) {
         for (std::int64_t j = m_2_bra.row[i]; j < m_2_bra.row[i + 1]; ++j) {
            if (m_2_bra.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(i, m_2_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateQuntumNumberDifference(i, m_3_ket.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_4_ket.row_dim; ++i) {
         for (std::int64_t j = m_4_ket.row[i]; j < m_4_ket.row[i + 1]; ++j) {
            if (m_4_ket.val[j] != 0.0) {
               delta_sector_set_m4.emplace(CalculateQuntumNumberDifference(i, m_4_ket.col[j]));
            }
         }
      }
      
      std::vector<std::tuple<std::pair<int, double>, std::pair<int, double>, std::pair<int, double>>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               for (const auto &del_sec_m4: delta_sector_set_m4) {
                  const std::pair<int, double> del_sec_m1_m2 = {del_sec_m1.first + del_sec_m2.first, del_sec_m1.second + del_sec_m2.second};
                  const std::pair<int, double> del_sec_m3_m4 = {del_sec_m3.first + del_sec_m4.first, del_sec_m3.second + del_sec_m4.second};
                  if (del_sec_m1_m2 == del_sec_m3_m4) {
                     target_sector_set.push_back({
                        {del_sec_m1.first    + total_electron_, del_sec_m1.second    + 0.5*total_2sz_},
                        {del_sec_m1_m2.first + total_electron_, del_sec_m1_m2.second + 0.5*total_2sz_},
                        {del_sec_m4.first    + total_electron_, del_sec_m4.second    + 0.5*total_2sz_}
                     });
                  }
               }
            }
         }
      }
      
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   static CRS CreateOnsiteOperatorCUp(const double magnitude_lspin) {
      
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*4;
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(1.0);
         matrix.col.push_back(dim_onsite_lspin + element);
         matrix.row[element + 1] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(1.0);
         matrix.col.push_back(3*dim_onsite_lspin + element);
         matrix.row[element + 1 + 2*dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 3*dim_onsite_lspin] = matrix.col.size();
      }
      matrix.tag = sparse_matrix::CRSTag::FERMION;
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorCDown(const double magnitude_lspin) {
      
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*4;
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(1.0);
         matrix.col.push_back(element + 2*dim_onsite_lspin);
         matrix.row[element + 1] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(-1.0);
         matrix.col.push_back(element + 3*dim_onsite_lspin);
         matrix.row[element + 1 + dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 2*dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 2*dim_onsite_lspin] = matrix.col.size();
      }
      matrix.tag = sparse_matrix::CRSTag::FERMION;
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorSzL(const double magnitude_lspin) {
      const int magnitude_2lspin    = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         for (int row_l = 0; row_l < dim_onsite_lspin; ++row_l) {
            const RealType val = magnitude_lspin - row_l;
            if (val != 0.0) {
               matrix.val.push_back(val);
               matrix.col.push_back(row_l + row_c*dim_onsite_lspin);
            }
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
      }
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorSpL(const double magnitude_lspin) {
      const int magnitude_2lspin    = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         for (int row_l = 1; row_l < dim_onsite_lspin; ++row_l) {
            matrix.val.push_back(std::sqrt((magnitude_lspin + 1)*2*row_l - row_l*(row_l + 1)));
            matrix.col.push_back(row_l + row_c*dim_onsite_lspin);
            matrix.row[row_l + row_c*dim_onsite_lspin] = matrix.col.size();
         }
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorSmL(const double magnitude_lspin) {
      const int magnitude_2lspin    = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         for (int row_l = 1; row_l < dim_onsite_lspin; ++row_l) {
            matrix.val.push_back(std::sqrt((magnitude_lspin + 1)*2*row_l - row_l*(row_l + 1)));
            matrix.col.push_back(row_l - 1 + row_c*dim_onsite_lspin);
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorSxL(const double magnitude_lspin) {
      const int magnitude_2lspin    = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         int a = 0;
         int b = 1;
         
         matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[1 + row_c*dim_onsite_lspin] = matrix.col.size();
         
         for (int row_l = 1; row_l < dim_onsite_lspin - 1; ++row_l) {
            a = row_l;
            b = row_l - 1;
            matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            a = row_l;
            b = row_l + 1;
            matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
         
         a = dim_onsite_lspin - 1;
         b = dim_onsite_lspin - 2;
         
         matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[dim_onsite_lspin + row_c*dim_onsite_lspin] = matrix.col.size();
      }
      return matrix;
   }
   
   static CRS CreateOnsiteOperatoriSyL(const double magnitude_lspin) {
      const int magnitude_2lspin    = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         int a = 0;
         int b = 1;
         
         matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[1 + row_c*dim_onsite_lspin] = matrix.col.size();
         
         for (int row_l = 1; row_l < dim_onsite_lspin - 1; ++row_l) {
            a = row_l;
            b = row_l - 1;
            
            matrix.val.push_back(-0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            a = row_l;
            b = row_l + 1;
            matrix.val.push_back(0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
         
         a = dim_onsite - 1;
         b = dim_onsite - 2;
         
         matrix.val.push_back(-0.5*std::sqrt((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1)));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[dim_onsite_lspin + row_c*dim_onsite_lspin] = matrix.col.size();
      }
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorCUpDagger(const double magnitude_lspin) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCUp(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorCDownDagger(const double magnitude_lspin) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCDown(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorNCUp(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(magnitude_lspin),
                                                         1.0, CreateOnsiteOperatorCUp(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorNCDown(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(magnitude_lspin),
                                                         1.0, CreateOnsiteOperatorCDown(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorNC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixSum(1.0, CreateOnsiteOperatorNCUp(magnitude_lspin),
                                                     1.0, CreateOnsiteOperatorNCDown(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorSxC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorSpC(magnitude_lspin),
                                                     0.5, CreateOnsiteOperatorSmC(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatoriSyC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixSum(+0.5, CreateOnsiteOperatorSpC(magnitude_lspin),
                                                     -0.5, CreateOnsiteOperatorSmC(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorSzC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixSum(+0.5, CreateOnsiteOperatorNCUp(magnitude_lspin),
                                                     -0.5,CreateOnsiteOperatorNCDown(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorSpC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(magnitude_lspin),
                                                         1.0, CreateOnsiteOperatorCDown(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorSmC(const double magnitude_lspin) {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(magnitude_lspin),
                                                         1.0, CreateOnsiteOperatorCUp(magnitude_lspin));
   }
   
   static CRS CreateOnsiteOperatorSCSL(const double magnitude_lspin) {
      const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin);
      const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin);
      const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin);
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin);
      const CRS spcsml = sparse_matrix::CalculateMatrixMatrixProduct(0.5, spc, 1.0, sml);
      const CRS smcspl = sparse_matrix::CalculateMatrixMatrixProduct(0.5, smc, 1.0, spl);
      const CRS szcszl = sparse_matrix::CalculateMatrixMatrixProduct(1.0, szc, 1.0, szl);
      const CRS xy = sparse_matrix::CalculateMatrixMatrixSum(1.0, spcsml, 1.0, smcspl);
      return sparse_matrix::CalculateMatrixMatrixSum(1.0, szcszl, 1.0, xy);
   }
   
   inline std::pair<int, double> CalculateQuntumNumberDifference(const int row, const int col) const {
      const int row_electron = CalculateBasisOnsiteElectron(row);
      const int col_electron = CalculateBasisOnsiteElectron(col);
      const int row_lspin    = CalculateBasisOnsiteLSpin(row);
      const int col_lspin    = CalculateBasisOnsiteLSpin(col);
      const auto diff_electron = BaseU1Electron_1D<RealType>::CalculateQuntumNumberDifference(row_electron, col_electron);
      const auto diff_lspin    = BaseU1Spin_1D<RealType>::CalculateQuntumNumberDifference(row_lspin, col_lspin);
      return {diff_electron.first, diff_electron.second + diff_lspin};
   }
      
   inline int    GetSystemSize()           const { return system_size_;            }
   inline int    GetDimOnsite()            const { return dim_onsite_;             }
   inline int    GetTotal2Sz()             const { return total_2sz_;              }
   inline double GetTotalSz()              const { return 0.5*total_2sz_;          }
   inline int    GetNumConservedQuantity() const { return num_conserved_quantity_; }
   inline double GetMagnitudeSpin()        const { return 0.5*magnitude_2lspin_;    }
   inline int    GetTotalElectron()        const { return total_electron_; }
   
   inline const CRS &GetOnsiteOperatorCUp()         const { return onsite_operator_c_up_;   }
   inline const CRS &GetOnsiteOperatorCDown()       const { return onsite_operator_c_down_; }
   inline const CRS &GetOnsiteOperatorCUpDagger()   const { return onsite_operator_c_up_dagger_; }
   inline const CRS &GetOnsiteOperatorCDownDagger() const { return onsite_operator_c_down_dagger_; }
   inline const CRS &GetOnsiteOperatorNCUp()        const { return onsite_operator_nc_up_; }
   inline const CRS &GetOnsiteOperatorNCDown()      const { return onsite_operator_nc_down_; }
   inline const CRS &GetOnsiteOperatorNC()          const { return onsite_operator_nc_; }
   inline const CRS &GetOnsiteOperatorSxC ()        const { return onsite_operator_sxc_ ; }
   inline const CRS &GetOnsiteOperatoriSyC()        const { return onsite_operator_isyc_; }
   inline const CRS &GetOnsiteOperatorSzC ()        const { return onsite_operator_szc_ ; }
   inline const CRS &GetOnsiteOperatorSpC ()        const { return onsite_operator_spc_ ; }
   inline const CRS &GetOnsiteOperatorSmC ()        const { return onsite_operator_smc_ ; }
   
   inline const CRS &GetOnsiteOperatorSxL()  const { return onsite_operator_sxl_; }
   inline const CRS &GetOnsiteOperatoriSyL() const { return onsite_operator_isyl_; }
   inline const CRS &GetOnsiteOperatorSzL()  const { return onsite_operator_szl_; }
   inline const CRS &GetOnsiteOperatorSpL()  const { return onsite_operator_spl_; }
   inline const CRS &GetOnsiteOperatorSmL()  const { return onsite_operator_sml_; }
   inline const CRS &GetOnsiteOperatorSCSL() const { return onsite_operator_scsl_; }
   
   inline const std::unordered_set<int> &GetCalculatedEigenvectorSet() const {
      return calculated_eigenvector_set_;
   }
   
   inline const std::vector<std::int64_t> &GetBasis(const std::pair<int, double> &quantum_number) const {
      return bases_.at({quantum_number.first, utility::DoubleTheNumber(quantum_number.second)});
   }
   
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetBasisInv(const std::pair<int, double> &quantum_number) const {
      return bases_inv_.at({quantum_number.first, utility::DoubleTheNumber(quantum_number.second)});
   }
   
   inline const std::vector<std::int64_t> &GetTargetBasis() const {
      return bases_.at({total_electron_, total_2sz_});
   }
   
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetTargetBasisInv() const {
      return bases_inv_.at({total_electron_, total_2sz_});
   }
   
protected:
   CRS onsite_operator_c_up_;
   CRS onsite_operator_c_down_;
   CRS onsite_operator_c_up_dagger_;
   CRS onsite_operator_c_down_dagger_;
   CRS onsite_operator_nc_up_;
   CRS onsite_operator_nc_down_;
   CRS onsite_operator_nc_;
   CRS onsite_operator_sxc_;
   CRS onsite_operator_isyc_;
   CRS onsite_operator_szc_;
   CRS onsite_operator_spc_;
   CRS onsite_operator_smc_;
   CRS onsite_operator_sxl_;
   CRS onsite_operator_isyl_;
   CRS onsite_operator_szl_;
   CRS onsite_operator_spl_;
   CRS onsite_operator_sml_;
   CRS onsite_operator_scsl_;
   
   const int dim_onsite_electron_    = 4;
   const int magnitude_2cspin_       = 1;
   const int num_conserved_quantity_ = 2;
   int system_size_      = 0;
   int total_2sz_        = 0;
   int total_electron_   = 0;
   int magnitude_2lspin_ = 1;
   int dim_onsite_lspin_ = magnitude_2lspin_ + 1;
   int dim_onsite_       = dim_onsite_electron_*dim_onsite_lspin_;
   
   std::unordered_set<int> calculated_eigenvector_set_;
   
   std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, utility::pair_hash> bases_;
   std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, utility::pair_hash> bases_inv_;
   
   void SetOnsiteOperator() {
      onsite_operator_c_up_   = CreateOnsiteOperatorCUp(0.5*magnitude_2lspin_);
      onsite_operator_c_down_ = CreateOnsiteOperatorCDown(0.5*magnitude_2lspin_);
      onsite_operator_c_up_dagger_   = CreateOnsiteOperatorCUpDagger(0.5*magnitude_2lspin_);
      onsite_operator_c_down_dagger_ = CreateOnsiteOperatorCDownDagger(0.5*magnitude_2lspin_);
      onsite_operator_nc_up_   = CreateOnsiteOperatorNCUp(0.5*magnitude_2lspin_);
      onsite_operator_nc_down_ = CreateOnsiteOperatorNCDown(0.5*magnitude_2lspin_);
      onsite_operator_nc_   = CreateOnsiteOperatorNC(0.5*magnitude_2lspin_);
      onsite_operator_sxc_  = CreateOnsiteOperatorSxC(0.5*magnitude_2lspin_);
      onsite_operator_isyc_ = CreateOnsiteOperatoriSyC(0.5*magnitude_2lspin_);
      onsite_operator_szc_  = CreateOnsiteOperatorSzC(0.5*magnitude_2lspin_);
      onsite_operator_spc_  = CreateOnsiteOperatorSpC(0.5*magnitude_2lspin_);
      onsite_operator_smc_  = CreateOnsiteOperatorSmC(0.5*magnitude_2lspin_);
      onsite_operator_sxl_  = CreateOnsiteOperatorSxL(0.5*magnitude_2lspin_);
      onsite_operator_isyl_ = CreateOnsiteOperatoriSyL(0.5*magnitude_2lspin_);
      onsite_operator_szl_  = CreateOnsiteOperatorSzL(0.5*magnitude_2lspin_);
      onsite_operator_spl_  = CreateOnsiteOperatorSpL(0.5*magnitude_2lspin_);
      onsite_operator_sml_  = CreateOnsiteOperatorSmL(0.5*magnitude_2lspin_);
      onsite_operator_scsl_ = CreateOnsiteOperatorSCSL(0.5*magnitude_2lspin_);
   }
   
   inline int CalculateBasisOnsiteElectron(const int basis_onsite) const {
      return basis_onsite/dim_onsite_lspin_;
   }
   
   inline int CalculateBasisOnsiteLSpin(const int basis_onsite) const {
      return basis_onsite%dim_onsite_lspin_;
   }
   
};




}
}



#endif /* COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_1D_HPP_ */
