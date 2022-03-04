//  Copyright 2022 Kohei Suzuki
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Created by Kohei Suzuki on 2021/12/11.
//

#ifndef COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_HPP_
#define COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_HPP_

#include "../blas/all.hpp"
#include "../utility/all.hpp"
#include "../type/all.hpp"
#include "base_u1_spin.hpp"
#include "base_u1_electron.hpp"

#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

//! @brief The base class for spin-electron systems with the U(1) symmetry.
//! @tparam RealType The type of real values.
template<typename RealType>
class BaseU1SpinElectron {
   
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
   //------------------------------------------------------------------
   //------------------------Private Type Alias------------------------
   //------------------------------------------------------------------
   //! @brief Alias of HalfInt type.
   using HalfInt = type::HalfInt;
      
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = type::CRS<RealType>;
      
public:
   //------------------------------------------------------------------
   //------------------------Public Type Alias-------------------------
   //------------------------------------------------------------------
   //! @brief Alias of quantum number (total electron, total sz) pair.
   using QType = std::pair<int, HalfInt>;
   
   //! @brief Alias of quantum number hash.
   using QHash = utility::IntHalfIntHash;
   
   //! @brief Alias of RealType.
   using ValueType = RealType;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of BaseU1SpinElectron class.
   BaseU1SpinElectron() {
      SetOnsiteOperator();
   }
      
   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   BaseU1SpinElectron(const HalfInt magnitude_lspin, const int total_electron) {
      SetTotalElectron(total_electron);
      SetMagnitudeLSpin(magnitude_lspin);
   }
   

   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   //! @param total_sz The total sz
   //! \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{s}^{z}_{i}+\hat{S}^{z}_{i}\rangle \f$
   BaseU1SpinElectron(const HalfInt magnitude_lspin,
                      const int total_electron,
                      const HalfInt total_sz) {
      SetTotalElectron(total_electron);
      SetTotalSz(total_sz);
      SetMagnitudeLSpin(magnitude_lspin);
   }
   
   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------   
   //! @brief Set target Hilbert space specified by the total sz to be diagonalized.
   //! @param total_sz The total sz
   //! \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{s}^{z}_{i}+\hat{S}^{z}_{i}\rangle \f$
   void SetTotalSz(const HalfInt total_sz) {
      total_sz_ = total_sz;
   }
   
   //! @brief Set the number of total electrons.
   //! @param total_electron The number of total electrons
   //! \f$ \hat{N}_{\rm e}=\sum^{N}_{i=1}\hat{n}_{i} \f$
   void SetTotalElectron(const int total_electron) {
      if (total_electron < 0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "total_electron must be a non-negative integer" << std::endl;
         throw std::runtime_error(ss.str());
      }
      total_electron_ = total_electron;
   }
   
   //! @brief Set the magnitude of the spin \f$ S \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   void SetMagnitudeLSpin(const HalfInt magnitude_lspin) {
      if (magnitude_lspin <= 0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Please set magnitude_spin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_lspin_ != magnitude_lspin) {
         magnitude_lspin_  = magnitude_lspin;
         dim_onsite_lspin_ = 2*magnitude_lspin + 1;
         dim_onsite_       = dim_onsite_lspin_*dim_onsite_electron_;
         SetOnsiteOperator();
      }
   }
   
   //! @brief Calculate the number of electrons from the input onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The number of electrons.
   int CalculateNumElectron(const int basis_onsite) const {
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      const int basis_onsite_electron = basis_onsite/dim_onsite_lspin_;
      
      if (basis_onsite_electron == 0) {
         return 0;
      }
      else if (basis_onsite_electron == 1 || basis_onsite_electron == 2) {
         return 1;
      }
      else if (basis_onsite_electron == 3) {
         return 2;
      }
      else {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Invalid onsite basis" << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Print the onsite bases.
   void PrintBasisOnsite() const {
      for (int row = 0; row < dim_onsite_; ++row) {
         std::string b_ele = "None";
         if (CalculateBasisOnsiteElectron(row) == 0) {
            b_ele = "|vac>";
         }
         else if (CalculateBasisOnsiteElectron(row) == 1) {
            b_ele = "|↑>";
         }
         else if (CalculateBasisOnsiteElectron(row) == 2) {
            b_ele = "|↓>";
         }
         else if (CalculateBasisOnsiteElectron(row) == 3) {
            b_ele = "|↑↓>";
         }
         std::cout << "row " << row << ": " << b_ele << "|Sz=";
         std::cout << magnitude_lspin_ - CalculateBasisOnsiteLSpin(row) << ">" << std::endl;
      }
   }
   
   //! @brief Calculate sectors generated by an onsite operator.
   //! @tparam IntegerType Value type int row and col.
   //! @param row The row in the matrix representation of an onsite operator.
   //! @param col The column in the matrix representation of an onsite operator.
   //! @return Sectors generated by an onsite operator.
   template<typename IntegerType>
   QType CalculateQNumber(const IntegerType row, const IntegerType col) const {
      static_assert(std::is_integral<IntegerType>::value, "Template parameter IntegerType must be integer type");
      if (row < 0 || col < 0 || dim_onsite_ <= row || dim_onsite_ <= col) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Invalid parameters" << std::endl;
         throw std::runtime_error(ss.str());
      }
      // Calculate total sz from local spin.
      const int row_lspin = row%dim_onsite_lspin_;
      const int col_lspin = col%dim_onsite_lspin_;
      const int total_sz_from_lspin = col_lspin - row_lspin;
      
      // Calculate total electron and total sz from electrons.
      const int row_electron = row/dim_onsite_lspin_;
      const int col_electron = col/dim_onsite_lspin_;
      if (row_electron < 0 || col_electron < 0 || dim_onsite_electron_ <= row_electron || dim_onsite_electron_ <= col_electron) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Invalid parameters" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (row_electron == col_electron && 0 <= row_electron && row_electron < 4 && 0 <= col_electron && col_electron < 4) {
         return {+0 + total_electron_, +0.0 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 0 && col_electron == 1) {
         return {-1 + total_electron_, -0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 0 && col_electron == 2) {
         return {-1 + total_electron_, +0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 0 && col_electron == 3) {
         return {-2 + total_electron_, +0.0 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 1 && col_electron == 0) {
         return {+1 + total_electron_, +0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 1 && col_electron == 2) {
         return {+0 + total_electron_, +1.0 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 1 && col_electron == 3) {
         return {-1 + total_electron_, +0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 2 && col_electron == 0) {
         return {+1 + total_electron_, -0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 2 && col_electron == 1) {
         return {+0 + total_electron_, -1.0 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 2 && col_electron == 3) {
         return {-1 + total_electron_, -0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 3 && col_electron == 0) {
         return {+2 + total_electron_, +0.0 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 3 && col_electron == 1) {
         return {+1 + total_electron_, -0.5 + total_sz_from_lspin + total_sz_};
      }
      else if (row_electron == 3 && col_electron == 2) {
         return {+1 + total_electron_, +0.5 + total_sz_from_lspin + total_sz_};
      }
      else {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "The dimenstion of the matrix must be 4";
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$
   //! and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @param flag_display_info If true, display the progress status. Set ture by default.
   std::vector<std::int64_t> GenerateBasis(const int system_size,
                                           const QType &quantum_number,
                                           const bool flag_display_info = true) const {
      const auto start = std::chrono::system_clock::now();
      
      if (!ValidateQNumber(quantum_number)) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "Invalid parameters (system_size or total_electron or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (flag_display_info) {
         std::cout << "Generating Basis..." << std::flush;
      }
      
      const int total_electron = quantum_number.first;
      const HalfInt total_sz   = quantum_number.second;

      std::vector<std::int64_t> site_constant_global(system_size);
      std::vector<std::int64_t> site_constant_electron(system_size);
      std::vector<std::int64_t> site_constant_lspin(system_size);
      for (int site = 0; site < system_size; ++site) {
         site_constant_global[site]   = static_cast<std::int64_t>(std::pow(dim_onsite_         , site));
         site_constant_electron[site] = static_cast<std::int64_t>(std::pow(dim_onsite_electron_, site));
         site_constant_lspin[site]    = static_cast<std::int64_t>(std::pow(dim_onsite_lspin_   , site));
      }
      
      const int max_n_up_down = static_cast<int>(total_electron/2);
      const std::vector<std::vector<std::int64_t>> binom = utility::GenerateBinomialTable(system_size);
      std::vector<int> q_number_spin_vec;
      std::vector<int> q_number_n_vac_vec;
      std::vector<int> q_number_n_up_vec;
      std::vector<int> q_number_n_down_vec;
      std::vector<int> q_number_n_up_down_vec;
      std::vector<std::int64_t> bias_basis;
      bias_basis.push_back(0);

      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         for (int n_up = 0; n_up <= total_electron - 2*n_up_down; ++n_up) {
            const int n_down = total_electron - 2*n_up_down - n_up;
            const int n_vac  = system_size - n_up - n_down - n_up_down;
            if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
               const HalfInt total_sz_electron = 0.5*(n_up - n_down);
               if (BaseU1Spin<RealType>::ValidateQNumber(system_size, magnitude_lspin_, total_sz - total_sz_electron)) {
                  q_number_spin_vec.push_back(2*total_sz - (n_up - n_down));
                  q_number_n_vac_vec.push_back(n_vac);
                  q_number_n_up_vec.push_back(n_up );
                  q_number_n_down_vec.push_back(n_down);
                  q_number_n_up_down_vec.push_back(n_up_down);
                  const std::int64_t dim_electron = binom[system_size][n_up]*binom[system_size - n_up][n_down]*binom[system_size - n_up - n_down][n_up_down];
                  const std::int64_t dim_lspin    = BaseU1Spin<RealType>::CalculateTargetDim(system_size, magnitude_lspin_, total_sz - total_sz_electron);
                  bias_basis.push_back(dim_electron*dim_lspin);
               }
            }
         }
      }
      
      for (std::size_t i = 1; i < bias_basis.size(); ++i) {
         bias_basis[i] += bias_basis[i - 1];
      }
      
      const std::int64_t dim_target_global = CalculateTargetDim(total_electron, total_sz);

      if (bias_basis.back() != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         throw std::runtime_error(ss.str());
      }

      //Generate spin bases
      std::vector<int> temp_q_number_spin_vec = q_number_spin_vec;
      std::sort(temp_q_number_spin_vec.begin(), temp_q_number_spin_vec.end());
      temp_q_number_spin_vec.erase(std::unique(temp_q_number_spin_vec.begin(), temp_q_number_spin_vec.end()), temp_q_number_spin_vec.end());
      std::unordered_map<int, std::vector<std::int64_t>> spin_bases;
#pragma omp parallel for
      for (std::size_t i = 0; i < temp_q_number_spin_vec.size(); ++i) {
         const int total_2sz_lspin    = temp_q_number_spin_vec[i];
         const int shifted_2sz        = system_size*magnitude_lspin_ - total_2sz_lspin/2;
         const std::int64_t dim_target_lspin = BaseU1Spin<RealType>::CalculateTargetDim(system_size, magnitude_lspin_, total_2sz_lspin/2);
         std::vector<std::vector<int>> partition_integers = utility::GenerateIntegerPartition(shifted_2sz, 2*magnitude_lspin_);
         auto &spin_basis = spin_bases[total_2sz_lspin];
         spin_basis.reserve(dim_target_lspin);
         for (auto &&integer_list: partition_integers) {
            const bool condition1 = (0 < integer_list.size()) && (static_cast<int>(integer_list.size()) <= system_size);
            const bool condition2 = (integer_list.size() == 0) && (shifted_2sz == 0);
            if (condition1 || condition2) {
               for (int j = static_cast<int>(integer_list.size()); j < system_size; ++j) {
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
            ss << "Unknown error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
            throw std::runtime_error(ss.str());
         }
         std::sort(spin_basis.begin(), spin_basis.end());
      }
      
      //Generate global bases
      const std::int64_t loop_size = static_cast<std::int64_t>(q_number_spin_vec.size());
      std::vector<std::int64_t> basis;
      basis.reserve(dim_target_global);

#pragma omp parallel for
      for (std::int64_t i = 0; i < loop_size; ++i) {
         const int n_vac     = q_number_n_vac_vec[i];
         const int n_up      = q_number_n_up_vec[i];
         const int n_down    = q_number_n_down_vec[i];
         const int n_up_down = q_number_n_up_down_vec[i];
         const int total_2sz_lspin = 2*q_number_spin_vec[i];
         std::int64_t count = bias_basis[i];
         std::vector<int> basis_list_electron(system_size);
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
               basis[count++] = basis_global_electron + basis_global_lspin;
            }
         } while (std::next_permutation(basis_list_electron.begin(), basis_list_electron.end()));
      }
      
      basis.shrink_to_fit();
      if (static_cast<std::int64_t>(basis.size()) != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(basis.begin(), basis.end());
      
      if (flag_display_info) {
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/blas::TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
      }
      return basis;
   }

   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool ValidateQNumber(const int system_size, const QType &quantum_number) const {
      return ValidateQNumber(system_size, magnitude_lspin_, quantum_number.first, quantum_number.second);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const int system_size, const QType &quantum_number) const {
      return CalculateTargetDim(system_size, magnitude_lspin_, quantum_number.first, quantum_number.second);
   }
   

   //------------------------------------------------------------------
   //----------------------Static Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool ValidateQNumber(const int system_size, const HalfInt magnitude_lspin, const int total_electron, const HalfInt total_sz) {
      if (system_size <= 0 || magnitude_lspin <= 0 || total_electron < 0) {
         return false;
      }
      const int total_2sz        = 2*total_sz;
      const int magnitude_2lspin = 2*magnitude_lspin;
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
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The dimension of the target Hilbert space.
   static std::int64_t CalculateTargetDim(const int system_size, const HalfInt magnitude_lspin, const int total_electron, const HalfInt total_sz) {
      if (!ValidateQNumber(system_size, magnitude_lspin, total_electron, total_sz)) {
         return 0;
      }
      const std::vector<std::vector<std::int64_t>> binom = utility::GenerateBinomialTable(system_size);
      const int max_n_up_down = static_cast<int>(total_electron/2);
      std::int64_t dim = 0;
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         for (int n_up = 0; n_up <= total_electron - 2*n_up_down; ++n_up) {
            const int n_down = total_electron - 2*n_up_down - n_up;
            const int n_vac  = system_size - n_up - n_down - n_up_down;
            if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
               const double total_sz_electron  = 0.5*(n_up - n_down);
               const std::int64_t dim_electron = binom[system_size][n_up]*binom[system_size - n_up][n_down]*binom[system_size - n_up - n_down][n_up_down];
               const std::int64_t dim_lspin    = BaseU1Spin<RealType>::CalculateTargetDim(system_size, magnitude_lspin, total_sz - total_sz_electron);
               dim += dim_electron*dim_lspin;
            }
         }
      }
      return dim;
   }

   
   //! @brief Generate the annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorCUp(const HalfInt magnitude_lspin) {
      
      const int magnitude_2lspin = 2*magnitude_lspin;
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
         matrix.val.push_back(RealType{1.0});
         matrix.col.push_back(dim_onsite_lspin + element);
         matrix.row[element + 1] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(RealType{1.0});
         matrix.col.push_back(3*dim_onsite_lspin + element);
         matrix.row[element + 1 + 2*dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 3*dim_onsite_lspin] = matrix.col.size();
      }
      matrix.tag = type::CRSTag::FERMION;
      return matrix;
   }
   
   //! @brief Generate the annihilation operator for the electrons with the down spin \f$ \hat{c}_{\downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorCDown(const HalfInt magnitude_lspin) {
      
      const int magnitude_2lspin = 2*magnitude_lspin;
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
         matrix.val.push_back(RealType{1.0});
         matrix.col.push_back(element + 2*dim_onsite_lspin);
         matrix.row[element + 1] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.val.push_back(-RealType{1.0});
         matrix.col.push_back(element + 3*dim_onsite_lspin);
         matrix.row[element + 1 + dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 2*dim_onsite_lspin] = matrix.col.size();
      }
      for (int element = 0; element < dim_onsite_lspin; ++element) {
         matrix.row[element + 1 + 3*dim_onsite_lspin] = matrix.col.size();
      }
      matrix.tag = type::CRSTag::FERMION;
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{z}\f$.
   static CRS CreateOnsiteOperatorSzL(const HalfInt magnitude_lspin) {
      const int magnitude_2lspin    = 2*magnitude_lspin;
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
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{+}\f$.
   static CRS CreateOnsiteOperatorSpL(const HalfInt magnitude_lspin) {
      const int magnitude_2lspin    = 2*magnitude_lspin;
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         for (int row_l = 1; row_l < dim_onsite_lspin; ++row_l) {
            matrix.val.push_back(std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*2*row_l - row_l*(row_l + 1))));
            matrix.col.push_back(row_l + row_c*dim_onsite_lspin);
            matrix.row[row_l + row_c*dim_onsite_lspin] = matrix.col.size();
         }
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{-}\f$.
   static CRS CreateOnsiteOperatorSmL(const HalfInt magnitude_lspin) {
      const int magnitude_2lspin    = 2*magnitude_lspin;
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         for (int row_l = 1; row_l < dim_onsite_lspin; ++row_l) {
            matrix.val.push_back(std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*2*row_l - row_l*(row_l + 1))));
            matrix.col.push_back(row_l - 1 + row_c*dim_onsite_lspin);
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{x}\f$.
   static CRS CreateOnsiteOperatorSxL(const HalfInt magnitude_lspin) {
      const int magnitude_2lspin    = 2*magnitude_lspin;
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         int a = 0;
         int b = 1;
         
         matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[1 + row_c*dim_onsite_lspin] = matrix.col.size();
         
         for (int row_l = 1; row_l < dim_onsite_lspin - 1; ++row_l) {
            a = row_l;
            b = row_l - 1;
            matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            a = row_l;
            b = row_l + 1;
            matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
         
         a = dim_onsite_lspin - 1;
         b = dim_onsite_lspin - 2;
         
         matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[dim_onsite_lspin + row_c*dim_onsite_lspin] = matrix.col.size();
      }
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{S}^{y}\f$.
   static CRS CreateOnsiteOperatoriSyL(const HalfInt magnitude_lspin) {
      const int magnitude_2lspin    = 2*magnitude_lspin;
      const int dim_onsite_lspin    = magnitude_2lspin + 1;
      const int dim_onsite_electron = 4;
      const int dim_onsite          = dim_onsite_lspin*dim_onsite_electron;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row_c = 0; row_c < dim_onsite_electron; ++row_c) {
         int a = 0;
         int b = 1;
         
         matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[1 + row_c*dim_onsite_lspin] = matrix.col.size();
         
         for (int row_l = 1; row_l < dim_onsite_lspin - 1; ++row_l) {
            a = row_l;
            b = row_l - 1;
            
            matrix.val.push_back(-RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            a = row_l;
            b = row_l + 1;
            matrix.val.push_back(RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
            matrix.col.push_back(b + row_c*dim_onsite_lspin);
            
            matrix.row[row_l + 1 + row_c*dim_onsite_lspin] = matrix.col.size();
         }
         
         a = dim_onsite - 1;
         b = dim_onsite - 2;
         
         matrix.val.push_back(-RealType{0.5}*std::sqrt(static_cast<RealType>((magnitude_lspin + 1)*(a + b + 1) - (a + 1)*(b + 1))));
         matrix.col.push_back(b + row_c*dim_onsite_lspin);
         matrix.row[dim_onsite_lspin + row_c*dim_onsite_lspin] = matrix.col.size();
      }
      return matrix;
   }
   
   //! @brief Generate the creation operator for the electrons with the up spin
   //! \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorCUpDagger(const HalfInt magnitude_lspin) {
      return type::CalculateTransposedMatrix(CreateOnsiteOperatorCUp(magnitude_lspin));
   }
   
   //! @brief Generate the creation operator for the electrons with the down spin
   //! \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorCDownDagger(const HalfInt magnitude_lspin) {
      return type::CalculateTransposedMatrix(CreateOnsiteOperatorCDown(magnitude_lspin));
   }
   
   //! @brief Generate the number operator for the electrons with the up spin
   //! \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorNCUp(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin)*CreateOnsiteOperatorCUp(magnitude_lspin);
   }
   
   //! @brief Generate the number operator for the electrons with the down spin
   //! \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorNCDown(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin)*CreateOnsiteOperatorCDown(magnitude_lspin);
   }
   
   //! @brief Generate the number operator for the electrons
   //! \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}\f$.
   static CRS CreateOnsiteOperatorNC(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorNCUp(magnitude_lspin) + CreateOnsiteOperatorNCDown(magnitude_lspin);
   }
   
   //! @brief Generate the spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   static CRS CreateOnsiteOperatorSxC(const HalfInt magnitude_lspin) {
      return RealType{0.5}*(CreateOnsiteOperatorSpC(magnitude_lspin) + CreateOnsiteOperatorSmC(magnitude_lspin));
   }
   
   //! @brief Generate the spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   static CRS CreateOnsiteOperatoriSyC(const HalfInt magnitude_lspin) {
      return RealType{0.5}*(CreateOnsiteOperatorSpC(magnitude_lspin) - CreateOnsiteOperatorSmC(magnitude_lspin));
   }
   
   //! @brief Generate the spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   static CRS CreateOnsiteOperatorSzC(const HalfInt magnitude_lspin) {
      return RealType{0.5}*(CreateOnsiteOperatorNCUp(magnitude_lspin) - CreateOnsiteOperatorNCDown(magnitude_lspin));
   }
   
   //! @brief Generate the raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   static CRS CreateOnsiteOperatorSpC(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin)*CreateOnsiteOperatorCDown(magnitude_lspin);
   }
   
   //! @brief Generate the lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   static CRS CreateOnsiteOperatorSmC(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin)*CreateOnsiteOperatorCUp(magnitude_lspin);
   }
   
   //! @brief Generate \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}=\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}+\hat{s}^{z}\hat{S}^{z}\f$
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}\f$.
   static CRS CreateOnsiteOperatorSCSL(const HalfInt magnitude_lspin) {
      const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin);
      const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin);
      const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin);
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin);
      return szc*szl + RealType{0.5}*(spc*sml + smc*spl);
   }
   
   static CRS CreateOnsiteOperatorSp(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorSpC(magnitude_lspin) + CreateOnsiteOperatorSpL(magnitude_lspin);
   }
   
   static CRS CreateOnsiteOperatorSm(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorSmC(magnitude_lspin) + CreateOnsiteOperatorSmL(magnitude_lspin);
   }
   
   static CRS CreateOnsiteOperatorSz(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorSzC(magnitude_lspin) + CreateOnsiteOperatorSzL(magnitude_lspin);
   }
   
   static CRS CreateOnsiteOperatorSx(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatorSxC(magnitude_lspin) + CreateOnsiteOperatorSxL(magnitude_lspin);
   }
   
   static CRS CreateOnsiteOperatoriSy(const HalfInt magnitude_lspin) {
      return CreateOnsiteOperatoriSyC(magnitude_lspin) + CreateOnsiteOperatoriSyL(magnitude_lspin);
   }

   //! @brief Get dimension of the local Hilbert space, \f$ 4*(2S+1)\f$.
   //! @return The dimension of the local Hilbert space, \f$ 4*(2S+1)\f$.
   inline int GetDimOnsite() const { return dim_onsite_; }
   
   inline int GetDimOnsiteELectron() const { return dim_onsite_electron_; }
   
   inline int GetDimOnsiteLSpin() const { return dim_onsite_lspin_; }
   
   
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline HalfInt GetTotalSz() const { return total_sz_; }

   //! @brief Get the magnitude of the local spin \f$ S\f$.
   //! @return The magnitude of the spin \f$ S\f$.
   inline HalfInt GetMagnitudeLSpin() const { return magnitude_lspin_; }
   
   //! @brief Get the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @return The total electrons.
   inline int GetTotalElectron() const { return total_electron_; }
   
   //! @brief Get the annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}_{\uparrow}\f$.
   inline const CRS &GetOnsiteOperatorCUp() const { return onsite_operator_c_up_; }
   
   //! @brief Get the annihilation operator for the electrons with the down spin \f$ \hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{c}_{\downarrow}\f$.
   inline const CRS &GetOnsiteOperatorCDown() const { return onsite_operator_c_down_; }
   
   //! @brief Get the creation operator for the electrons with the up spin
   //! \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   inline const CRS &GetOnsiteOperatorCUpDagger() const { return onsite_operator_c_up_dagger_; }
   
   //! @brief Get the creation operator for the electrons with the down spin
   //! \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   inline const CRS &GetOnsiteOperatorCDownDagger() const { return onsite_operator_c_down_dagger_; }
   
   //! @brief Get the number operator for the electrons with the up spin
   //! \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\uparrow}\f$.
   inline const CRS &GetOnsiteOperatorNCUp() const { return onsite_operator_nc_up_; }
   
   //! @brief Get the number operator for the electrons with the down spin
   //! \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\downarrow}\f$.
   inline const CRS &GetOnsiteOperatorNCDown() const { return onsite_operator_nc_down_; }
   
   //! @brief Get the number operator for the electrons
   //! \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}\f$.
   inline const CRS &GetOnsiteOperatorNC() const { return onsite_operator_nc_; }
   
   //! @brief Get the spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   inline const CRS &GetOnsiteOperatorSxC() const { return onsite_operator_sxc_ ; }
   
   //! @brief Get the spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   inline const CRS &GetOnsiteOperatoriSyC() const { return onsite_operator_isyc_; }
   
   //! @brief Get the spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   inline const CRS &GetOnsiteOperatorSzC() const { return onsite_operator_szc_ ; }
   
   //! @brief Get the raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   inline const CRS &GetOnsiteOperatorSpC() const { return onsite_operator_spc_ ; }
   
   //! @brief Get the lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   inline const CRS &GetOnsiteOperatorSmC() const { return onsite_operator_smc_ ; }
   
   //! @brief Get the spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   //! @return The matrix of \f$ \hat{S}^{x}\f$.
   inline const CRS &GetOnsiteOperatorSxL()  const { return onsite_operator_sxl_; }
   
   //! @brief Get the spin-\f$ S\f$ operator of the local spin for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @return The matrix of \f$ i\hat{S}^{y}\f$.
   inline const CRS &GetOnsiteOperatoriSyL() const { return onsite_operator_isyl_; }
   
   //! @brief Get the spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   //! @return The matrix of \f$ \hat{S}^{z}\f$.
   inline const CRS &GetOnsiteOperatorSzL()  const { return onsite_operator_szl_; }
   
   //! @brief Get the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   //! @return The matrix of \f$ \hat{S}^{+}\f$.
   inline const CRS &GetOnsiteOperatorSpL()  const { return onsite_operator_spl_; }
   
   //! @brief Get the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   //! @return The matrix of \f$ \hat{S}^{-}\f$.
   inline const CRS &GetOnsiteOperatorSmL()  const { return onsite_operator_sml_; }
   
   //! @brief Get \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}=\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}+\hat{s}^{z}\hat{S}^{z}\f$
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}\f$.
   inline const CRS &GetOnsiteOperatorSCSL() const { return onsite_operator_scsl_; }
   
   inline const CRS &GetOnsiteOperatorSp()  const { return onsite_operator_sp ; }
   inline const CRS &GetOnsiteOperatorSm()  const { return onsite_operator_sm ; }
   inline const CRS &GetOnsiteOperatorSz()  const { return onsite_operator_sz ; }
   inline const CRS &GetOnsiteOperatorSx()  const { return onsite_operator_sx ; }
   inline const CRS &GetOnsiteOperatoriSy() const { return onsite_operator_isy; }

   
protected:
   //------------------------------------------------------------------
   //---------------------Private Member Variables---------------------
   //------------------------------------------------------------------
   //! @brief The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   HalfInt total_sz_ = 0;
   
   //! @brief The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   int total_electron_ = 0;
   
   //! @brief The dimension of the local Hilbert space for the electrons, 4.
   const int dim_onsite_electron_ = 4;
   
   //! @brief The magnitude of the local spin \f$ S\f$.
   HalfInt magnitude_lspin_ = 0.5;
   
   //! @brief The dimension of the local Hilbert space, \f$ 2S + 1\f$.
   int dim_onsite_lspin_ = 2*magnitude_lspin_ + 1;
   
   //! @brief The dimension of the local Hilbert space, \f$ 4\times (2S + 1) \f$.
   int dim_onsite_ = dim_onsite_electron_*dim_onsite_lspin_;
   
   //! @brief The annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   CRS onsite_operator_c_up_;
   
   //! @brief The annihilation operator for the electrons with the down spin \f$ \hat{c}_{\downarrow}\f$.
   CRS onsite_operator_c_down_;
   
   //! @brief The creation operator for the electrons with the up spin \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   CRS onsite_operator_c_up_dagger_;
   
   //! @brief The the creation operator for the electrons with the down spin \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   CRS onsite_operator_c_down_dagger_;
   
   //! @brief The number operator for the electrons with the up spin \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   CRS onsite_operator_nc_up_;
   
   //! @brief The number operator for the electrons with the down spin \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   CRS onsite_operator_nc_down_;
   
   //! @brief The number operator for the electrons \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   CRS onsite_operator_nc_;
   
   //! @brief The spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   CRS onsite_operator_sxc_;
   
   //! @brief The spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   CRS onsite_operator_isyc_;
   
   //! @brief The spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   CRS onsite_operator_szc_;
   
   //! @brief The raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   CRS onsite_operator_spc_;
   
   //! @brief The lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   CRS onsite_operator_smc_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   CRS onsite_operator_sxl_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   CRS onsite_operator_isyl_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   CRS onsite_operator_szl_;
   
   //! @brief The spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   CRS onsite_operator_spl_;
   
   //! @brief The spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   CRS onsite_operator_sml_;
   
   //! @brief The correlation between the electron spin and local spin
   //! \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}=\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}+\hat{s}^{z}\hat{S}^{z}\f$
   CRS onsite_operator_scsl_;
   
   CRS onsite_operator_sp ;
   CRS onsite_operator_sm ;
   CRS onsite_operator_sz ;
   CRS onsite_operator_sx ;
   CRS onsite_operator_isy;

   
   //------------------------------------------------------------------
   //----------------------Private Member Functions--------------------
   //------------------------------------------------------------------
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_c_up_          = CreateOnsiteOperatorCUp(magnitude_lspin_);
      onsite_operator_c_down_        = CreateOnsiteOperatorCDown(magnitude_lspin_);
      onsite_operator_c_up_dagger_   = CreateOnsiteOperatorCUpDagger(magnitude_lspin_);
      onsite_operator_c_down_dagger_ = CreateOnsiteOperatorCDownDagger(magnitude_lspin_);
      onsite_operator_nc_up_         = CreateOnsiteOperatorNCUp(magnitude_lspin_);
      onsite_operator_nc_down_       = CreateOnsiteOperatorNCDown(magnitude_lspin_);
      onsite_operator_nc_            = CreateOnsiteOperatorNC(magnitude_lspin_);
      onsite_operator_sxc_           = CreateOnsiteOperatorSxC(magnitude_lspin_);
      onsite_operator_isyc_          = CreateOnsiteOperatoriSyC(magnitude_lspin_);
      onsite_operator_szc_           = CreateOnsiteOperatorSzC(magnitude_lspin_);
      onsite_operator_spc_           = CreateOnsiteOperatorSpC(magnitude_lspin_);
      onsite_operator_smc_           = CreateOnsiteOperatorSmC(magnitude_lspin_);
      onsite_operator_sxl_           = CreateOnsiteOperatorSxL(magnitude_lspin_);
      onsite_operator_isyl_          = CreateOnsiteOperatoriSyL(magnitude_lspin_);
      onsite_operator_szl_           = CreateOnsiteOperatorSzL(magnitude_lspin_);
      onsite_operator_spl_           = CreateOnsiteOperatorSpL(magnitude_lspin_);
      onsite_operator_sml_           = CreateOnsiteOperatorSmL(magnitude_lspin_);
      onsite_operator_scsl_          = CreateOnsiteOperatorSCSL(magnitude_lspin_);
      onsite_operator_sp             = CreateOnsiteOperatorSp(magnitude_lspin_);
      onsite_operator_sm             = CreateOnsiteOperatorSm(magnitude_lspin_);
      onsite_operator_sz             = CreateOnsiteOperatorSz(magnitude_lspin_);
      onsite_operator_sx             = CreateOnsiteOperatorSx(magnitude_lspin_);
      onsite_operator_isy            = CreateOnsiteOperatoriSy(magnitude_lspin_);
   }
   
   //! @brief Calculate onsite basis for the electrons from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The onsite basis for the electrons.
   inline int CalculateBasisOnsiteElectron(const int basis_onsite) const {
      return basis_onsite/dim_onsite_lspin_;
   }
   
   //! @brief Calculate onsite basis for the loca spins from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The onsite basis for the local spins.
   inline int CalculateBasisOnsiteLSpin(const int basis_onsite) const {
      return basis_onsite%dim_onsite_lspin_;
   }
   
};

} //namespace model
} //namespace compnal



#endif /* COMPNAL_MODEL_BASE_U1SPIN_ELECTRON_HPP_ */
