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
//  Created by Kohei Suzuki on 2021/12/20.
//

#ifndef COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_
#define COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"
#include "base_u1_spin_1d.hpp"
#include "base_u1_electron_1d.hpp"

#include <unordered_map>
#include <unordered_set>

namespace compnal {
namespace model {

//! @brief The base class for one-dimensional spin-multiband-electron systems with the U(1) symmetry.
//! @tparam RealType The type of real values.
template<typename RealType>
class BaseU1SpinMultiElectrons_1D {
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = type::CRS<RealType>;
   
   //! @brief Alias of quantum number (total electron list, total sz) pair.
   using QType = std::pair<std::vector<int>, double>;
   
public:
   
   //! @brief The type of real values.
   using ValueType = RealType;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of BaseU1SpinMultiElectrons_1D class.
   BaseU1SpinMultiElectrons_1D() {
      SetOnsiteOperator();
   }
   
   //! @brief Constructor of BaseU1SpinMultiElectrons_1D class.
   //! @param system_size The system size \f$ N \f$.
   explicit BaseU1SpinMultiElectrons_1D(const int system_size): BaseU1SpinMultiElectrons_1D() {
      SetSystemSize(system_size);
   }
   
   //! @brief Constructor of BaseU1SpinMultiElectrons_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   BaseU1SpinMultiElectrons_1D(const int system_size, const double magnitude_lspin): BaseU1SpinMultiElectrons_1D(system_size) {
      SetMagnitudeLSpin(magnitude_lspin);
   }
   
   //! @brief Constructor of BaseU1SpinMultiElectrons_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   BaseU1SpinMultiElectrons_1D(const int system_size,
                               const std::vector<int> &total_electron): BaseU1SpinMultiElectrons_1D(system_size) {
      SetTotalElectron(total_electron);
   }
   
   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   BaseU1SpinMultiElectrons_1D(const int system_size,
                               const double magnitude_lspin,
                               const std::vector<int> &total_electron): BaseU1SpinMultiElectrons_1D(system_size, magnitude_lspin) {
      SetTotalElectron(total_electron);
   }
   
   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz
   //! \f$ \langle\hat{S}^{z}_{\rm tot}\rangle =
   //! \sum^{N}_{i=1}\left(\hat{S}^{z}_{i} + \sum_{\alpha}\hat{s}^{z}_{i,\alpha}\right)\f$.
   BaseU1SpinMultiElectrons_1D(const int system_size,
                               const double magnitude_lspin,
                               const std::vector<int> &total_electron,
                               const double total_sz): BaseU1SpinMultiElectrons_1D(system_size, magnitude_lspin, total_electron) {
      SetTotalSz(total_sz);
   }
   
   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------
   //! @brief Set system size.
   //! @param system_size The system size \f$ N \f$.
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
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
   
   //! @brief Set target Hilbert space specified by the total sz to be diagonalized.
   //! @param total_sz The total sz
   //! \f$ \langle\hat{S}^{z}_{\rm tot}\rangle =
   //! \sum^{N}_{i=1}\left(\hat{S}^{z}_{i} + \sum_{\alpha}\hat{s}^{z}_{i,\alpha}\right)\f$.
   void SetTotalSz(const double total_sz) {
      const int total_2sz = utility::DoubleHalfInteger(total_sz);
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the number of total electrons.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   void SetTotalElectron(const std::vector<int> &total_electron) {
      if (num_electron_orbital_ != static_cast<int>(total_electron.size())) {
         num_electron_orbital_ = static_cast<int>(total_electron.size());
         dim_onsite_all_electrons_ = static_cast<int>(std::pow(dim_onsite_electron_, num_electron_orbital_));
         dim_onsite_ = dim_onsite_lspin_*dim_onsite_all_electrons_;
         SetOnsiteOperator();
      }
      
      if (total_electron_.size() < total_electron.size()) {
         for (std::size_t i = 0; i < total_electron_.size(); ++i) {
            if (total_electron_.at(i) != total_electron.at(i)) {
               total_electron_.at(i) = total_electron.at(i);
            }
         }
         for (std::size_t i = total_electron_.size(); i < total_electron.size(); ++i) {
            total_electron_.push_back(total_electron.at(i));
         }
         calculated_eigenvector_set_.clear();
      }
      else {
         for (std::size_t i = 0; i < total_electron.size(); ++i) {
            if (total_electron_.at(i) != total_electron.at(i)) {
               total_electron_.at(i) = total_electron.at(i);
               calculated_eigenvector_set_.clear();
            }
         }
      }
   }
   
   //! @brief Set the magnitude of the spin \f$ S \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   void SetMagnitudeLSpin(const double magnitude_lspin) {
      const int magnitude_2lspin = utility::DoubleHalfInteger(magnitude_lspin);
      if (magnitude_2lspin <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Please set magnitude_2lspin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_2lspin_ != magnitude_2lspin) {
         magnitude_2lspin_ = magnitude_2lspin;
         dim_onsite_lspin_ = magnitude_2lspin + 1;
         dim_onsite_       = dim_onsite_lspin_*dim_onsite_all_electrons_;
         SetOnsiteOperator();
         bases_.clear();
         bases_inv_.clear();
         calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @param level Energy level.
   void SetCalculatedEigenvectorSet(const int level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   //! @brief Calculate the number of electrons from the onsite electron basis.
   //! @param basis_onsite The onsite basis.
   //! @return The number of electrons
   int CalculateNumElectron(const int basis_onsite) const {
      int num_electron = 0;
      for (int o = 0; o < num_electron_orbital_; ++o) {
         const int basis_electron_onsite = CalculateBasisOnsiteElectron(basis_onsite, o);
         if (basis_electron_onsite == 0) {
            num_electron += 0;
         }
         else if (basis_electron_onsite == 1 || basis_electron_onsite == 2) {
            num_electron += 1;
         }
         else if (basis_electron_onsite == 3) {
            num_electron += 2;
         }
         else {
            std::stringstream ss;
            ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
            throw std::runtime_error(ss.str());
         }
      }
      return num_electron;
   }
   
   //! @brief Print the onsite bases.
   void PrintBasisOnsite() const {
      const double magnitude_lspin = magnitude_2lspin_/2.0;
      for (int row = 0; row < dim_onsite_; ++row) {
         std::vector<std::string> b_ele;
         for (int o = 0; o < num_electron_orbital_; ++o) {
            if (CalculateBasisOnsiteElectron(row, o) == 0) {
               b_ele.push_back("|vac>");
            }
            else if (CalculateBasisOnsiteElectron(row, o) == 1) {
               b_ele.push_back("| ↑ >");
            }
            else if (CalculateBasisOnsiteElectron(row, o) == 2) {
               b_ele.push_back("| ↓ >");
            }
            else if (CalculateBasisOnsiteElectron(row, o) == 3) {
               b_ele.push_back("|↑↓ >");
            }
            else {
               std::stringstream ss;
               ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
               throw std::runtime_error(ss.str());
            }
         }
         std::cout << "row " << std::setw(2) << row << ": ";
         for (const auto &it: b_ele) {
            std::cout << it;
         }
         std::cout << "|Sz=" << std::showpos << magnitude_lspin - CalculateBasisOnsiteLSpin(row) << ">" << std::noshowpos << std::endl;
      }
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const QType &quantum_number) const {
      return isValidQNumber(system_size_, 0.5*magnitude_2lspin_, quantum_number.first, quantum_number.second);
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const std::vector<int> &total_electron, const double total_sz) const {
      return isValidQNumber(system_size_, 0.5*magnitude_2lspin_, total_electron, total_sz);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim(0.5*total_2sz_);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const double total_sz) const {
      return CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, total_electron_, total_sz);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const std::vector<int> &total_electron, const double total_sz) const {
      return CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, total_electron, total_sz);
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1 The matrix of an onsite operator.
   //! @param m_2 The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<QType> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      // TODO: Check input matrics
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m1;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m2;
      for (std::int64_t i = 0; i < m_1.row_dim; ++i) {
         for (std::int64_t j = m_1.row[i]; j < m_1.row[i + 1]; ++j) {
            if (m_1.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateTargetQuantumNumber(i, m_1.col[j]));
            }
         }
      }
      for (std::int64_t i = 0; i < m_2.row_dim; ++i) {
         for (std::int64_t j = m_2.row[i]; j < m_2.row[i + 1]; ++j) {
            if (m_2.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateTargetQuantumNumber(i, m_2.col[j]));
            }
         }
      }
      std::vector<QType> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            const bool c1 = isValidQNumber(del_sec_m1.first, del_sec_m1.second);
            if (del_sec_m1 == del_sec_m2 && c1) {
               target_sector_set.push_back(del_sec_m1);
            }
         }
      }
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1_bra The matrix of an onsite operator.
   //! @param m_2_ket The matrix of an onsite operator.
   //! @param m_3_ket The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<std::pair<QType, QType>> GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_ket, const CRS &m_3_ket) const {
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m1;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m2;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m3;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateTargetQuantumNumber(i, m_1_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_ket.row_dim; ++i) {
         for (std::int64_t j = m_2_ket.row[i]; j < m_2_ket.row[i + 1]; ++j) {
            if (m_2_ket.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateTargetQuantumNumber(i, m_2_ket.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateTargetQuantumNumber(i, m_3_ket.col[j]));
            }
         }
      }
      
      std::vector<std::pair<QType, QType>> target_sector_set;
      
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               QType del_sec_m2_m3;
               for (std::size_t i = 0; i < del_sec_m2.first.size(); ++i) {
                  del_sec_m2_m3.first.push_back(del_sec_m2.first[i] + del_sec_m3.first[i]);
               }
               del_sec_m2_m3.second = del_sec_m2.second + del_sec_m3.second;
               const bool c1 = isValidQNumber(del_sec_m1.first, del_sec_m1.second);
               const bool c2 = isValidQNumber(del_sec_m3.first, del_sec_m3.second);
               if (del_sec_m1 == del_sec_m2_m3 && c1 && c2) {
                  target_sector_set.push_back({del_sec_m1, del_sec_m3});
               }
            }
         }
      }
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1_bra The matrix of an onsite operator.
   //! @param m_2_bra The matrix of an onsite operator.
   //! @param m_3_ket The matrix of an onsite operator.
   //! @param m_4_ket The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<std::tuple<QType, QType, QType>>
   GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_bra, const CRS &m_3_ket, const CRS &m_4_ket) const {
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m1;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m2;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m3;
      std::unordered_set<QType, utility::VecIntHash> delta_sector_set_m4;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateTargetQuantumNumber(i, m_1_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_bra.row_dim; ++i) {
         for (std::int64_t j = m_2_bra.row[i]; j < m_2_bra.row[i + 1]; ++j) {
            if (m_2_bra.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateTargetQuantumNumber(i, m_2_bra.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateTargetQuantumNumber(i, m_3_ket.col[j]));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_4_ket.row_dim; ++i) {
         for (std::int64_t j = m_4_ket.row[i]; j < m_4_ket.row[i + 1]; ++j) {
            if (m_4_ket.val[j] != 0.0) {
               delta_sector_set_m4.emplace(CalculateTargetQuantumNumber(i, m_4_ket.col[j]));
            }
         }
      }
      
      std::vector<std::tuple<QType, QType, QType>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               for (const auto &del_sec_m4: delta_sector_set_m4) {
                  QType del_sec_m1_m2;
                  QType del_sec_m3_m4;
                  for (std::size_t i = 0; i < del_sec_m1.first.size(); ++i) {
                     del_sec_m1_m2.first.push_back(del_sec_m1.first[i] + del_sec_m2.first[i]);
                     del_sec_m3_m4.first.push_back(del_sec_m3.first[i] + del_sec_m4.first[i]);
                  }
                  del_sec_m1_m2.second = del_sec_m1.second + del_sec_m2.second;
                  del_sec_m3_m4.second = del_sec_m3.second + del_sec_m4.second;
                  const bool c1 = isValidQNumber(del_sec_m1.first   , del_sec_m1.second   );
                  const bool c2 = isValidQNumber(del_sec_m1_m2.first, del_sec_m1_m2.second);
                  const bool c3 = isValidQNumber(del_sec_m4.first   , del_sec_m4.second   );
                  if (del_sec_m1_m2 == del_sec_m3_m4 && c1 && c2 && c3) {
                     target_sector_set.push_back({del_sec_m1, del_sec_m1_m2, del_sec_m4});
                  }
               }
            }
         }
      }
      
      std::sort(target_sector_set.begin(), target_sector_set.end());
      target_sector_set.erase(std::unique(target_sector_set.begin(), target_sector_set.end()), target_sector_set.end());
      return target_sector_set;
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   void GenerateBasis() {
      GenerateBasis(total_electron_, 0.5*total_2sz_);
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   void GenerateBasis(const QType &quantum_number) {
      GenerateBasis(quantum_number.first, quantum_number.second);
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   void GenerateBasis(const std::vector<int> &total_electron, const double total_sz) {
      const auto start     = std::chrono::system_clock::now();
      const int  total_2sz = utility::DoubleHalfInteger(total_sz);
      if (bases_.count({total_electron, total_2sz}) != 0) {
         return;
      }
      std::cout << "Generating Basis..." << std::flush;
            
      std::vector<std::vector<std::vector<int>>> electron_configuration_list;
      std::vector<int> length_list;
      std::int64_t length = 1;
      for (const auto num_electron: total_electron) {
         const auto electron_configuration = GenerateElectronConfigurations(system_size_, num_electron);
         electron_configuration_list.push_back(electron_configuration);
         length_list.push_back(static_cast<int>(electron_configuration[0].size()));
         length *= electron_configuration[0].size();
      }
      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size_);
      const std::int64_t dim_target_global = CalculateTargetDim(total_electron, total_sz);

      std::vector<std::vector<int>> q_number_n_up_down;
      std::vector<std::vector<int>> q_number_n_up;
      std::vector<std::vector<int>> q_number_n_down;
      std::vector<std::vector<int>> q_number_n_vac;
      std::vector<std::vector<int>> q_number_ele_list;
      std::vector<int> q_number_spin_vec;
      std::vector<std::int64_t> total_electron_dim_list;
      std::vector<std::vector<std::int64_t>> each_electron_basis_dim_prod;
      std::vector<std::vector<std::int64_t>> each_electron_basis_dim;
      std::vector<std::int64_t> bias_basis;
      bias_basis.push_back(0);
      
      for (std::int64_t i = 0; i < length; ++i) {
         int electron_2sz = 0;
         std::vector<int> temp_n_up_down;
         std::vector<int> temp_n_up;
         std::vector<int> temp_n_down;
         std::vector<int> temp_n_vac;
         std::int64_t electron_dim = 1;
         for (std::size_t j = 0; j < length_list.size(); ++j) {
            std::int64_t prod = 1;
            for (std::size_t k = j + 1; k < length_list.size(); ++k) {
               prod *= length_list[k];
            }
            const std::size_t index = (i/prod)%length_list[j];
            const int n_up_down = electron_configuration_list[j][0][index];
            const int n_up      = electron_configuration_list[j][1][index];
            const int n_down    = electron_configuration_list[j][2][index];
            const int n_vac     = electron_configuration_list[j][3][index];
            electron_2sz += n_up - n_down;
            temp_n_up_down.push_back(n_up_down);
            temp_n_up     .push_back(n_up     );
            temp_n_down   .push_back(n_down   );
            temp_n_vac    .push_back(n_vac    );
            electron_dim *= binom[system_size_][n_up]*binom[system_size_ - n_up][n_down]*binom[system_size_ - n_up - n_down][n_up_down];
         }
         const int spin_2sz = total_2sz - electron_2sz;
         if (BaseU1Spin_1D<RealType>::isValidQNumber(system_size_, 0.5*magnitude_2lspin_, 0.5*spin_2sz)) {
            q_number_n_up_down.push_back(temp_n_up_down);
            q_number_n_up     .push_back(temp_n_up     );
            q_number_n_down   .push_back(temp_n_down   );
            q_number_n_vac    .push_back(temp_n_vac    );
            q_number_spin_vec .push_back(spin_2sz);
            total_electron_dim_list.push_back(electron_dim);
            each_electron_basis_dim.emplace_back();
            bias_basis.push_back(electron_dim*BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, 0.5*spin_2sz));
            for (std::size_t j = 0; j < length_list.size(); ++j) {
               q_number_ele_list.push_back({temp_n_up_down[j], temp_n_up[j], temp_n_down[j], temp_n_vac[j]});
               each_electron_basis_dim.back().push_back(binom[system_size_][temp_n_up[j]]*
                                                        binom[system_size_ - temp_n_up[j]][temp_n_down[j]]*
                                                        binom[system_size_ - temp_n_up[j] - temp_n_down[j]][temp_n_up_down[j]]);
               
            }
         }
      }
      
      each_electron_basis_dim_prod.resize(each_electron_basis_dim.size());
      for (std::size_t i = 0; i < each_electron_basis_dim_prod.size(); ++i) {
         each_electron_basis_dim_prod[i].resize(total_electron.size());
         for (std::size_t j = 0; j < total_electron.size(); ++j) {
            std::int64_t prod = 1;
            for (std::size_t k = j + 1; k < total_electron.size(); ++k) {
               prod *= each_electron_basis_dim[i][k];
            }
            each_electron_basis_dim_prod[i][j] = prod;
         }
      }
      
      for (std::size_t i = 1; i < bias_basis.size(); ++i) {
         bias_basis[i] += bias_basis[i - 1];
      }
      
      if (bias_basis.back() != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::vector<std::int64_t> site_constant_global(system_size_);
      std::vector<std::int64_t> ele_constant_global(system_size_);

      for (int site = 0; site < system_size_; ++site) {
         site_constant_global[site] = static_cast<std::int64_t>(std::pow(dim_onsite_, site));
      }
      
      for (int o = 0; o < static_cast<int>(total_electron.size()); ++o) {
         ele_constant_global[o] = static_cast<std::int64_t>(std::pow(dim_onsite_electron_, o));
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
         const std::int64_t dim_target_lspin = BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size_, 0.5*magnitude_2lspin_, 0.5*total_2_sz_lspin);
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
            ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
            throw std::runtime_error(ss.str());
         }
         std::sort(spin_basis.begin(), spin_basis.end());
      }
      
      //Generate electron bases
      const std::int64_t loop_size = static_cast<std::int64_t>(q_number_spin_vec.size());
      std::vector<std::vector<int>> temp_q_number_ele_list = q_number_ele_list;
      std::unordered_map<std::vector<int>, std::vector<std::int64_t>, utility::VecHash> electron_bases;
      std::sort(temp_q_number_ele_list.begin(), temp_q_number_ele_list.end());
      temp_q_number_ele_list.erase(std::unique(temp_q_number_ele_list.begin(), temp_q_number_ele_list.end()), temp_q_number_ele_list.end());
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<std::unordered_map<std::vector<int>, std::vector<std::int64_t>, utility::VecHash>> temp_electron_bases(num_threads);
      
#pragma omp parallel for
      for (std::size_t i = 0; i < temp_q_number_ele_list.size(); ++i) {
         const int n_up_down = temp_q_number_ele_list[i][0];
         const int n_down    = temp_q_number_ele_list[i][1];
         const int n_up      = temp_q_number_ele_list[i][2];
         const int n_vac     = temp_q_number_ele_list[i][3];
         const int thread_num = omp_get_thread_num();
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
         auto &temp_basis = temp_electron_bases[thread_num][{n_up_down, n_up, n_down, n_vac}];
         do {
            std::int64_t basis_global_electron = 0;
            for (std::size_t j = 0; j < basis_list_electron.size(); ++j) {
               basis_global_electron += basis_list_electron[j]*site_constant_global[j]*dim_onsite_lspin_;
            }
            temp_basis.push_back(basis_global_electron);
         } while (std::next_permutation(basis_list_electron.begin(), basis_list_electron.end()));
      }
      
      for (auto &&basis: temp_electron_bases) {
         electron_bases.merge(basis);
      }
      
#else
      for (std::size_t i = 0; i < temp_q_number_ele_list.size(); ++i) {
         const int n_up_down = temp_q_number_ele_list[i][0];
         const int n_down    = temp_q_number_ele_list[i][1];
         const int n_up      = temp_q_number_ele_list[i][2];
         const int n_vac     = temp_q_number_ele_list[i][3];
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
            electron_bases[{n_up_down, n_up, n_down, n_vac}].push_back(basis_global_electron);
         } while (std::next_permutation(basis_list_electron.begin(), basis_list_electron.end()));
      }
#endif
      
      //Generate global bases
      std::vector<std::int64_t>().swap(bases_[{total_electron, total_2sz}]);
      auto &global_basis_ref = bases_.at({total_electron, total_2sz});
      global_basis_ref.resize(dim_target_global);
      
#pragma omp parallel for
      for (std::int64_t i = 0; i < loop_size; ++i) {
         std::int64_t count = bias_basis[i];
         const int total_2sz_lspin = q_number_spin_vec[i];
         const std::int64_t total_electron_dim = total_electron_dim_list[i];
         for (const auto &spin_basis: spin_bases.at(total_2sz_lspin)) {
            
            for (std::int64_t j = 0; j < total_electron_dim; ++j) {
               std::int64_t global_basis = spin_basis;
               for (std::size_t k = 0; k < total_electron.size(); ++k) {
                  const int n_up_down = q_number_n_up_down[i][k];
                  const int n_up      = q_number_n_up     [i][k];
                  const int n_down    = q_number_n_down   [i][k];
                  const int n_vac     = q_number_n_vac    [i][k];
                  std::size_t index = (j/each_electron_basis_dim_prod[i][k])%each_electron_basis_dim[i][k];
                  global_basis += electron_bases.at({n_up_down, n_up, n_down, n_vac})[index];
               }
               global_basis_ref[count++] = global_basis;
            }
         }
      }
      
      if (static_cast<std::int64_t>(bases_.at({total_electron, total_2sz}).size()) != dim_target_global) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
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
         ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "The same basis has been detected" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
      
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The dimension of the target Hilbert space.
   static std::int64_t CalculateTargetDim(const int system_size,
                                          const double magnitude_lspin,
                                          const std::vector<int> &total_electron,
                                          const double total_sz) {
      std::vector<std::vector<std::vector<int>>> electron_configuration_list;
      std::vector<int> length_list;
      std::int64_t length = 1;
      for (const auto num_electron: total_electron) {
         const auto electron_configuration = GenerateElectronConfigurations(system_size, num_electron);
         electron_configuration_list.push_back(electron_configuration);
         length_list.push_back(static_cast<int>(electron_configuration[0].size()));
         length *= electron_configuration[0].size();
      }
      const int total_2sz = utility::DoubleHalfInteger(total_sz);
      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size);
      std::int64_t dim = 0;
      for (std::int64_t i = 0; i < length; ++i) {
         int electron_2sz = 0;
         std::int64_t electron_dim = 1;
         for (std::size_t j = 0; j < length_list.size(); ++j) {
            std::int64_t prod = 1;
            for (std::size_t k = j + 1; k < length_list.size(); ++k) {
               prod *= length_list[k];
            }
            const std::size_t index = (i/prod)%length_list[j];
            const int n_up_down = electron_configuration_list[j][0][index];
            const int n_up      = electron_configuration_list[j][1][index];
            const int n_down    = electron_configuration_list[j][2][index];
            electron_2sz += n_up - n_down;
            electron_dim *= binom[system_size][n_up]*binom[system_size - n_up][n_down]*binom[system_size - n_up - n_down][n_up_down];
         }
         const int spin_2sz = total_2sz - electron_2sz;
         if (BaseU1Spin_1D<RealType>::isValidQNumber(system_size, magnitude_lspin, 0.5*spin_2sz)) {
            dim += electron_dim*BaseU1Spin_1D<RealType>::CalculateTargetDim(system_size, magnitude_lspin, 0.5*spin_2sz);
         }
      }
      return dim;
   }
      
   //! @brief Calculate difference of the number of total electrons and the total sz
   //! from the rows and columns in the matrix representation of an onsite operator.
   //! @param row The row in the matrix representation of an onsite operator.
   //! @param col The column in the matrix representation of an onsite operator.
   //! @return The differences of the total electron and the total sz.
   inline QType CalculateTargetQuantumNumber(const int row, const int col) const {
      std::vector<int> diff_electron_list;
      int target_total_sz = 0;
      for (int o = 0; o < num_electron_orbital_; ++o) {
         const int row_electron   = CalculateBasisOnsiteElectron(row, o);
         const int col_electron   = CalculateBasisOnsiteElectron(col, o);
         const auto diff_electron = BaseU1Electron_1D<RealType>::CalculateQuntumNumberDifference(row_electron, col_electron);
         diff_electron_list.push_back(diff_electron.first + total_electron_[o]);
         target_total_sz += diff_electron.second;
      }
      const int row_lspin = CalculateBasisOnsiteLSpin(row);
      const int col_lspin = CalculateBasisOnsiteLSpin(col);
      target_total_sz += BaseU1Spin_1D<RealType>::CalculateQuntumNumberDifference(row_lspin, col_lspin) + 0.5*total_2sz_;
      
      return {diff_electron_list, target_total_sz};
   }
   
   //! @brief Generate the annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}_{\alpha, \uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{c}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorCUp(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleHalfInteger(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         int num_electron = 0;
         for (int o = 0; o < orbital; ++o) {
            const int basis_electron_onsite = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
            num_electron += CalculateNumElectronFromElectronBasis(basis_electron_onsite);
         }
         int sign = 1;
         if (num_electron%2 == 1) {
            sign = -1;
         }
         
         for (int col = 0; col < dim_onsite; col++) {
            bool c_1 = true;
            for (int o = 0; o < num_orbital; ++o) {
               const int basis_1 = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
               const int basis_2 = CalculateBasisOnsiteElectron(col, magnitude_lspin, o, num_orbital);
               if (o != orbital && basis_1 != basis_2) {
                  c_1 = false;
                  break;
               }
            }
            
            int basis_row_electron = CalculateBasisOnsiteElectron(row, magnitude_lspin, orbital, num_orbital);
            int basis_col_electron = CalculateBasisOnsiteElectron(col, magnitude_lspin, orbital, num_orbital);
            
            const bool c_2 = (basis_row_electron == 0 && basis_col_electron == 1);
            const bool c_3 = (basis_row_electron == 2 && basis_col_electron == 3);
            
            if (c_1 && row%dim_onsite_lspin == col%dim_onsite_lspin && (c_2 || c_3)) {
               matrix.val.push_back(sign);
               matrix.col.push_back(col);
            }
      
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      matrix.tag = type::CRSTag::FERMION;
      return matrix;
   }
   
   //! @brief Generate the annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}_{\alpha, \downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{c}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorCDown(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleHalfInteger(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         int num_electron = 0;
         for (int o = 0; o < orbital; ++o) {
            const int basis_electron_onsite = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
            num_electron += CalculateNumElectronFromElectronBasis(basis_electron_onsite);
         }
         int sign_1 = 1;
         if (num_electron%2 == 1) {
            sign_1 = -1;
         }
         
         for (int col = 0; col < dim_onsite; col++) {
            bool c_1 = true;
            for (int o = 0; o < num_orbital; ++o) {
               const int basis_1 = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
               const int basis_2 = CalculateBasisOnsiteElectron(col, magnitude_lspin, o, num_orbital);
               if (o != orbital && basis_1 != basis_2) {
                  c_1 = false;
                  break;
               }
            }
            
            int basis_row_electron = CalculateBasisOnsiteElectron(row, magnitude_lspin, orbital, num_orbital);
            int basis_col_electron = CalculateBasisOnsiteElectron(col, magnitude_lspin, orbital, num_orbital);
            
            const bool c_2 = (basis_row_electron == 0 && basis_col_electron == 2);
            const bool c_3 = (basis_row_electron == 1 && basis_col_electron == 3);
            
            int sign_2 = 1;
            if (c_3) {
               sign_2 = -1;
            }
            
            if (c_1 && row%dim_onsite_lspin == col%dim_onsite_lspin && (c_2 || c_3)) {
               matrix.val.push_back(sign_1*sign_2);
               matrix.col.push_back(col);
            }
      
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      matrix.tag = type::CRSTag::FERMION;
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{z}\f$.
   static CRS CreateOnsiteOperatorSzL(const double magnitude_lspin, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleHalfInteger(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; ++row) {
         for (int col = 0; col < dim_onsite; ++col) {
            const int basis_row_lspin = row%dim_onsite_lspin;
            const int basis_col_lspin = col%dim_onsite_lspin;
            const int basis_row_total_electron = row/dim_onsite_lspin;
            const int basis_col_total_electron = col/dim_onsite_lspin;
            if (basis_row_total_electron == basis_col_total_electron && basis_row_lspin == basis_col_lspin) {
               matrix.val.push_back(magnitude_lspin - basis_row_lspin);
               matrix.col.push_back(col);
            }
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{+}\f$.
   static CRS CreateOnsiteOperatorSpL(const double magnitude_lspin, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleHalfInteger(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; ++row) {
         for (int col = 0; col < dim_onsite; ++col) {
            const int a = row%dim_onsite_lspin + 1;
            const int b = col%dim_onsite_lspin + 1;
            const int basis_row_total_electron = row/dim_onsite_lspin;
            const int basis_col_total_electron = col/dim_onsite_lspin;
            if (basis_row_total_electron == basis_col_total_electron && a + 1 == b) {
               matrix.val.push_back(std::sqrt( (magnitude_lspin + 1.0)*(a + b - 1) - a*b ) );
               matrix.col.push_back(col);
            }
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorCUpDagger(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorCDownDagger(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin
   //! \f$ \hat{n}_{\alpha, \uparrow}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorNCUp(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin
   //! \f$ \hat{n}_{\alpha, \downarrow}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorNCDown(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the number operator for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{n}_{\alpha}=\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorNC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorNCUp(magnitude_lspin, orbital, num_orbital) + CreateOnsiteOperatorNCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the number operator for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{n}=\sum_{\alpha}\left(\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\right)\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \sum_{\alpha}\hat{n}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorNCTot(const double magnitude_lspin, const int num_orbital) {
      const int dim = (utility::DoubleHalfInteger(magnitude_lspin) + 1)*static_cast<int>(std::pow(4, num_orbital));
      CRS out(dim, dim);
      for (int o = 0; o < num_orbital; ++o) {
         out = out + CreateOnsiteOperatorNC(magnitude_lspin, o, num_orbital);
      }
      return out;
   }
   
   //! @brief Generate the spin operator for the x-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{x}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //!  + \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{x}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSxC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital) + CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the spin operator for the y-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ i\hat{s}^{y}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //!  - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{s}^{y}_{\alpha}\f$.
   static CRS CreateOnsiteOperatoriSyC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital) - CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the spin operator for the z-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{z}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}
   //!  - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow})\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{z}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSzC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorNCUp(magnitude_lspin, orbital, num_orbital) - CreateOnsiteOperatorNCDown(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the raising operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{+}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{+}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSpC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the lowering operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{-}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{-}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSmC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{-}\f$.
   static CRS CreateOnsiteOperatorSmL(const double magnitude_lspin, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{x}\f$.
   static CRS CreateOnsiteOperatorSxL(const double magnitude_lspin, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital) + CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin
   //!  for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{S}^{y}\f$.
   static CRS CreateOnsiteOperatoriSyL(const double magnitude_lspin, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital) - CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}=
   //! \hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}+\hat{s}^{z}_{\alpha}\hat{S}^{z}\f$
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}\f$.
   static CRS CreateOnsiteOperatorSCSL(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital);
      const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital);
      const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin, orbital, num_orbital);
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin, num_orbital);
      return szc*szl + 0.5*(spc*sml + smc*spl);
   }
   
   //! @brief Generate \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}=
   //! \sum_{\alpha}\left(\hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}
   //!  +\hat{s}^{z}_{\alpha}\hat{S}^{z}\right)\f$
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}\f$.
   static CRS CreateOnsiteOperatorSCSLTot(const double magnitude_lspin, const int num_orbital) {
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin, num_orbital);
      CRS out(szl.row_dim, szl.col_dim);
      for (int o = 0; o < num_orbital; ++o) {
         const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin, o, num_orbital);
         const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin, o, num_orbital);
         const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin, o, num_orbital);
         out = out + szc*szl + 0.5*(spc*sml + smc*spl);
      }
      return out;
   }
   
   //---------------------------Access Private Member Functions---------------------------
   inline const std::vector<CRS> &GetOnsiteOperatorCUp() const { return onsite_operator_c_up_; }
   inline const std::vector<CRS> &GetOnsiteOperatorCDown() const { return onsite_operator_c_down_; }
   inline const std::vector<CRS> &GetOnsiteOperatorCUpDagger() const { return onsite_operator_c_up_dagger_; }
   inline const std::vector<CRS> &GetOnsiteOperatorCDownDagger() const { return onsite_operator_c_down_dagger_; }
   inline const std::vector<CRS> &GetOnsiteOperatorNCUp() const { return onsite_operator_nc_up_; }
   inline const std::vector<CRS> &GetOnsiteOperatorNCDown() const { return onsite_operator_nc_up_; }
   inline const std::vector<CRS> &GetOnsiteOperatorNC() const { return onsite_operator_nc_; }
   inline const std::vector<CRS> &GetOnsiteOperatorSxC() const { return onsite_operator_sxc_; }
   inline const std::vector<CRS> &GetOnsiteOperatoriSyC() const { return onsite_operator_isyc_; }
   inline const std::vector<CRS> &GetOnsiteOperatorSzC() const { return onsite_operator_szc_; }
   inline const std::vector<CRS> &GetOnsiteOperatorSpC() const { return onsite_operator_spc_; }
   inline const std::vector<CRS> &GetOnsiteOperatorSmC() const { return onsite_operator_smc_; }
   inline const std::vector<CRS> &GetOnsiteOperatorSCSL() const { return onsite_operator_scsl_;}

   inline const CRS &GetOnsiteOperatorNCTot() { return onsite_operator_nc_tot_; }
   inline const CRS &GetOnsiteOperatorCUp(const int orbital) const { return onsite_operator_c_up_.at(orbital); }
   inline const CRS &GetOnsiteOperatorCDown(const int orbital) const { return onsite_operator_c_down_.at(orbital); }
   inline const CRS &GetOnsiteOperatorCUpDagger(const int orbital) const { return onsite_operator_c_up_dagger_.at(orbital); }
   inline const CRS &GetOnsiteOperatorCDownDagger(const int orbital) const { return onsite_operator_c_down_dagger_.at(orbital); }
   inline const CRS &GetOnsiteOperatorNCUp(const int orbital) const { return onsite_operator_nc_up_.at(orbital); }
   inline const CRS &GetOnsiteOperatorNCDown(const int orbital) const { return onsite_operator_nc_up_.at(orbital); }
   inline const CRS &GetOnsiteOperatorNC(const int orbital) const { return onsite_operator_nc_.at(orbital); }
   inline const CRS &GetOnsiteOperatorSxC(const int orbital) const { return onsite_operator_sxc_.at(orbital); }
   inline const CRS &GetOnsiteOperatoriSyC(const int orbital) const { return onsite_operator_isyc_.at(orbital); }
   inline const CRS &GetOnsiteOperatorSzC(const int orbital) const { return onsite_operator_szc_.at(orbital); }
   inline const CRS &GetOnsiteOperatorSpC(const int orbital) const { return onsite_operator_spc_.at(orbital); }
   inline const CRS &GetOnsiteOperatorSmC(const int orbital) const { return onsite_operator_smc_.at(orbital); }
   inline const CRS &GetOnsiteOperatorSCSL(const int orbital) const { return onsite_operator_scsl_.at(orbital);}
   
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
   
   inline int GetDimOnsiteElectron() const { return dim_onsite_electron_; }
   
   inline const std::vector<int> &GetTotalElectron() const { return total_electron_; }
   
   inline int GetTotalElectron(const int orbital) const { return total_electron_.at(orbital); }

   inline int GetNumElectronOrbital() const { return num_electron_orbital_; }
   
   inline int GetDimOnsiteAllElectrons() const { return dim_onsite_all_electrons_; }
   
   //! @brief Get the system size \f$ N\f$.
   //! @return The system size \f$ N\f$.
   inline int GetSystemSize() const { return system_size_; }
   
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline double GetTotalSz() const { return 0.5*total_2sz_; }
   
   //! @brief Get the magnitude of the local spin \f$ S\f$.
   //! @return The magnitude of the spin \f$ S\f$.
   inline double GetMagnitudeLSpin() const { return 0.5*magnitude_2lspin_; }
   
   inline int GetDimOnsiteLSpin() const { return magnitude_2lspin_ + 1; }
   
   //! @brief Get dimension of the local Hilbert space, \f$ 4^{n_{\rm o}}*(2S+1)\f$.
   //! @return The dimension of the local Hilbert space, \f$ 4^{n_{\rm o}}*(2S+1)\f$.
   inline int GetDimOnsite() const { return dim_onsite_; }
   
   //! @brief Get calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @return calculated_eigenvector_set_.
   inline const std::unordered_set<int> &GetCalculatedEigenvectorSet() const {
      return calculated_eigenvector_set_;
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetBasis(const QType &quantum_number) const {
      return bases_.at({quantum_number.first, utility::DoubleHalfInteger(quantum_number.second)});
   }
   
   //! @brief Get inverse basis of the target Hilbert space space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return Inverse basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetBasisInv(const QType &quantum_number) const {
      return bases_inv_.at({quantum_number.first, utility::DoubleHalfInteger(quantum_number.second)});
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetTargetBasis() const {
      return bases_.at({total_electron_, total_2sz_});
   }
   
   //! @brief Get inverse basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Inverse basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetTargetBasisInv() const {
      return bases_inv_.at({total_electron_, total_2sz_});
   }
   
protected:
   //! @brief The annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin, \f$ \hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_c_up_;
   
   //! @brief The annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_c_down_;
   
   //! @brief The creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_c_up_dagger_;
   
   //! @brief The the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_c_down_dagger_;
   
   //! @brief The number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin,
   //! \f$ \hat{n}_{\alpha, \uparrow}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_nc_up_;
   
   //! @brief The number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin,
   //! \f$ \hat{n}_{\alpha, \downarrow}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_nc_down_;
   
   //! @brief The number operator for the electrons with the orbital
   //! \f$ \alpha \f$, \f$ \hat{n}_{\alpha}=\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_nc_;
   
   //! @brief The number operator for the electrons,
   //! \f$ \hat{n}=\sum_{\alpha}\left(\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\right)\f$.
   CRS onsite_operator_nc_tot_;
   
   //! @brief The spin operator for the x-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{x}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //! + \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   std::vector<CRS> onsite_operator_sxc_;
   
   //! @brief The spin operator for the y-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ i\hat{s}^{y}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //! - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   std::vector<CRS> onsite_operator_isyc_;
   
   //! @brief The spin operator for the z-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{z}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}
   //! - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow})\f$.
   std::vector<CRS> onsite_operator_szc_;
   
   //! @brief The raising operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{+}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_spc_;
   
   //! @brief The lowering operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{-}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_smc_;
   
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
   
   //! @brief The correlation between the electron with the orbital \f$ \alpha \f$ spin and local spin
   //! \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}=
   //! \hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}+\hat{s}^{z}_{\alpha}\hat{S}^{z}\f$
   std::vector<CRS> onsite_operator_scsl_;
   
   //! @brief The dimension of the local Hilbert space for the electrons, 4.
   const int dim_onsite_electron_ = 4;
   
   //! @brief The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   std::vector<int> total_electron_ = {0};
   
   //! @brief The number of electron orbitals \f$ n_{\rm o}\f$.
   int num_electron_orbital_ = static_cast<int>(total_electron_.size());
   
   //! @brief The dimension of the local Hilbert space for the total electrons, \f$ 4^{n_{\rm o}}\f$.
   int dim_onsite_all_electrons_ = static_cast<int>(std::pow(dim_onsite_electron_, num_electron_orbital_));
   
   //! @brief The system size.
   int system_size_ = 0;
   
   //! @brief Twice the number of the total sz
   //! \f$ 2\langle\hat{S}^{z}_{\rm tot}\rangle =
   //! 2\sum^{N}_{i=1}\left(\hat{S}^{z}_{i} + \sum_{\alpha}\hat{s}^{z}_{i,\alpha}\right)\f$.
   int total_2sz_ = 0;
   
   //! @brief The magnitude of the local spin \f$ S\f$.
   int magnitude_2lspin_ = 1;
   
   //! @brief The dimension of the local Hilbert space, \f$ 2S + 1\f$.
   int dim_onsite_lspin_ = magnitude_2lspin_ + 1;
   
   //! @brief The dimension of the local Hilbert space, \f$ 4^{n_{\rm o}}\times (2S + 1) \f$.
   int dim_onsite_ = dim_onsite_all_electrons_*dim_onsite_lspin_;
   
   //! @brief The calculated eigenvectors and eigenvalues.
   std::unordered_set<int> calculated_eigenvector_set_;
   
   //! @brief Bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons at each orbital \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! The first value of std::vector<int> stores twice the number of the total sz and remaining values correspond to the orbitals of the electrons.
   std::unordered_map<std::pair<std::vector<int>, int>, std::vector<std::int64_t>, utility::VecIntHash> bases_;
   
   //! @brief Inverse bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons at each orbital \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! The first value of std::vector<int> stores twice the number of the total sz and remaining values correspond to the orbitals of the electrons.
   std::unordered_map<std::pair<std::vector<int>, int>, std::unordered_map<std::int64_t, std::int64_t>, utility::VecIntHash> bases_inv_;
   
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_c_up_         .resize(num_electron_orbital_);
      onsite_operator_c_down_       .resize(num_electron_orbital_);
      onsite_operator_c_up_dagger_  .resize(num_electron_orbital_);
      onsite_operator_c_down_dagger_.resize(num_electron_orbital_);
      onsite_operator_nc_up_        .resize(num_electron_orbital_);
      onsite_operator_nc_down_      .resize(num_electron_orbital_);
      onsite_operator_nc_           .resize(num_electron_orbital_);
      onsite_operator_sxc_          .resize(num_electron_orbital_);
      onsite_operator_isyc_         .resize(num_electron_orbital_);
      onsite_operator_szc_          .resize(num_electron_orbital_);
      onsite_operator_spc_          .resize(num_electron_orbital_);
      onsite_operator_smc_          .resize(num_electron_orbital_);
      onsite_operator_scsl_         .resize(num_electron_orbital_);
      for (int o = 0; o < num_electron_orbital_; ++o) {
         onsite_operator_c_up_[o]          = CreateOnsiteOperatorCUp        (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_c_down_[o]        = CreateOnsiteOperatorCDown      (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_c_up_dagger_[o]   = CreateOnsiteOperatorCUpDagger  (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_c_down_dagger_[o] = CreateOnsiteOperatorCDownDagger(0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_nc_up_[o]         = CreateOnsiteOperatorNCUp       (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_nc_down_[o]       = CreateOnsiteOperatorNCDown     (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_nc_[o]            = CreateOnsiteOperatorNC         (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_sxc_[o]           = CreateOnsiteOperatorSxC        (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_isyc_[o]          = CreateOnsiteOperatoriSyC       (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_szc_[o]           = CreateOnsiteOperatorSzC        (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_spc_[o]           = CreateOnsiteOperatorSpC        (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_smc_[o]           = CreateOnsiteOperatorSmC        (0.5*magnitude_2lspin_, o, num_electron_orbital_);
         onsite_operator_scsl_[o]          = CreateOnsiteOperatorSCSL       (0.5*magnitude_2lspin_, o, num_electron_orbital_);
      }
      onsite_operator_nc_tot_ = CreateOnsiteOperatorNCTot(0.5*magnitude_2lspin_, num_electron_orbital_);
      onsite_operator_sxl_    = CreateOnsiteOperatorSxL  (0.5*magnitude_2lspin_, num_electron_orbital_);
      onsite_operator_isyl_   = CreateOnsiteOperatoriSyL (0.5*magnitude_2lspin_, num_electron_orbital_);
      onsite_operator_szl_    = CreateOnsiteOperatorSzL  (0.5*magnitude_2lspin_, num_electron_orbital_);
      onsite_operator_spl_    = CreateOnsiteOperatorSpL  (0.5*magnitude_2lspin_, num_electron_orbital_);
      onsite_operator_sml_    = CreateOnsiteOperatorSmL  (0.5*magnitude_2lspin_, num_electron_orbital_);
   }
      
   //! @brief Calculate onsite basis for the electrons from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @return The onsite basis for the electrons.
   inline int CalculateBasisOnsiteElectron(const int basis_onsite, const int orbital) const {
      const int num_inner_electron = num_electron_orbital_ - orbital - 1;
      int basis_electron_onsite = basis_onsite/dim_onsite_lspin_;
      for (int i = 0; i < num_inner_electron; ++i) {
         basis_electron_onsite /= dim_onsite_electron_;
      }
      return basis_electron_onsite%dim_onsite_electron_;
   }
   
   //! @brief Calculate onsite basis for the electrons from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The onsite basis for the electrons.
   inline static int CalculateBasisOnsiteElectron(const int basis_onsite,
                                                  const double magnitude_lspin,
                                                  const int orbital,
                                                  const int num_orbital) {
      const int dim_onsite_lspin = utility::DoubleHalfInteger(magnitude_lspin) + 1;
      const int num_inner_electron = num_orbital - orbital - 1;
      const int dim_onsite_electron = 4;
      int basis_electron_onsite = basis_onsite/dim_onsite_lspin;
      for (int i = 0; i < num_inner_electron; ++i) {
         basis_electron_onsite /= dim_onsite_electron;
      }
      return basis_electron_onsite%dim_onsite_electron;
   }
   
   //! @brief Calculate onsite basis for the loca spins from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The onsite basis for the local spins.
   inline int CalculateBasisOnsiteLSpin(const int basis_onsite) const {
      return basis_onsite%dim_onsite_lspin_;
   }
   
   //! @brief Calculate the number of electrons from the onsite electron basis.
   //! @param basis_electron_onsite The onsite electron basis.
   //! @return The number of electrons
   inline static int CalculateNumElectronFromElectronBasis(const int basis_electron_onsite) {
      if (basis_electron_onsite == 0) {
         return 0;
      }
      else if (basis_electron_onsite == 1 || basis_electron_onsite == 2) {
         return 1;
      }
      else if (basis_electron_onsite == 3) {
         return 2;
      }
      else {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
   static std::vector<std::vector<int>> GenerateElectronConfigurations(const int system_size, const int total_electron) {
      const int max_n_up_down = static_cast<int>(total_electron/2);
      std::vector<std::vector<int>> ele_configurations(4);
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         for (int n_up = 0; n_up <= total_electron - 2*n_up_down; ++n_up) {
            const int n_down = total_electron - 2*n_up_down - n_up;
            const int n_vac  = system_size - n_up - n_down - n_up_down;
            if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
               ele_configurations[0].push_back(n_up_down);
               ele_configurations[1].push_back(n_up     );
               ele_configurations[2].push_back(n_down   );
               ele_configurations[3].push_back(n_vac    );
            }
         }
      }
      return ele_configurations;
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool isValidQNumber(const int system_size, const double magnitude_lspin, const std::vector<int> &total_electron, const double total_sz) {
      const int total_2sz            = utility::DoubleHalfInteger(total_sz);
      const int magnitude_2lspin     = utility::DoubleHalfInteger(magnitude_lspin);
      const int total_total_electron = std::accumulate(total_electron.begin(), total_electron.end(), 0);
      
      for (const auto &electron: total_electron) {
         if (electron < 0 || 2*system_size < electron) {
            return false;
         }
      }
      
      const bool c1 = ((total_total_electron + system_size*magnitude_2lspin - total_2sz)%2 == 0);
      const bool c2 = (-(total_total_electron + system_size*magnitude_2lspin) <= total_2sz);
      const bool c3 = (total_2sz <= total_total_electron + system_size*magnitude_2lspin);
      
      if (!c1 && !c2 && !c3) {
         return false;
      }
      
      return true;
      
   }
   
};

}
}


#endif /* COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_ */
