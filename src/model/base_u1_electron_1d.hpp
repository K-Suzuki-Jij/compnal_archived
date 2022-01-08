//
//  base_u1_electron_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/27.
//

#ifndef COMPNAL_MODEL_BASE_U1ELECTRON_1D_HPP_
#define COMPNAL_MODEL_BASE_U1ELECTRON_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"

#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

//! @brief The base class for one-dimensional electron systems with the U(1) symmetry.
//! @tparam RealType The type of real values.
template<typename RealType>
class BaseU1Electron_1D {
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   //! @brief The type of real values.
   using ValueType = RealType;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of BaseU1Electron_1D class.
   BaseU1Electron_1D() {
      SetOnsiteOperator();
   }
   
   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param system_size The system size \f$ N \f$.
   explicit BaseU1Electron_1D(const int system_size): BaseU1Electron_1D() {
      SetSystemSize(system_size);
   }
   
   //! @brief Constructor of BaseU1Electron_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   BaseU1Electron_1D(const int system_size, const int total_electron): BaseU1Electron_1D(system_size) {
      SetTotalElectron(total_electron);
   }
   
   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------
   //! @brief Set system size.
   //! @param system_size The system size \f$ N \f$.
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
   
   //! @brief Set target Hilbert space specified by the total sz to be diagonalized.
   //! @param total_sz The total sz is the expectation value of the following operator:
   //! \f[ \hat{S}^{z}_{\rm tot}=\sum^{N}_{i=1}\hat{S}^{z}_{i} \f]
   void SetTotalSz(const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the number of total electrons.
   //! @param total_electron The number of total electrons is represented by the expectation value of the following operator:
   //! \f[ \hat{N}_{\rm e}=\sum^{N}_{i=1}\hat{n}_{i} \f]
   void SetTotalElectron(const int total_electron) {
      if (total_electron_ != total_electron) {
         total_electron_ = total_electron;
         calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @param level Energy level.
   void SetCalculatedEigenvectorSet(const std::int64_t level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const std::pair<int, double> &quantum_number) const {
      return isValidQNumber(system_size_, quantum_number.first, quantum_number.second);
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const int total_electron, const double total_sz) const {
      return isValidQNumber(system_size_, total_electron, total_sz);
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
      
      if (basis_onsite == 0) {
         return 0;
      }
      else if (basis_onsite == 1 || basis_onsite == 2) {
         return 1;
      }
      else if (basis_onsite == 3) {
         return 2;
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "Invalid onsite basis" << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Print the onsite bases.
   void PrintBasisOnsite() const {
      std::cout << "row " << 0 << ": |vac>" << std::endl;
      std::cout << "row " << 1 << ": |↑>"   << std::endl;
      std::cout << "row " << 2 << ": |↓>"   << std::endl;
      std::cout << "row " << 3 << ": |↑↓>"  << std::endl;
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim(system_size_, total_electron_, 0.5*total_2sz_);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const int total_electron, const double total_sz) const {
      return CalculateTargetDim(system_size_, total_electron, total_sz);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const std::pair<int, double> &quantum_number) const {
      return CalculateTargetDim(system_size_, quantum_number.first, quantum_number.second);
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   void GenerateBasis() {
      GenerateBasis({total_electron_, 0.5*total_2sz_});
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   void GenerateBasis(const std::pair<int, double> &quantum_number) {
      if (!isValidQNumber(quantum_number)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const auto start         = std::chrono::system_clock::now();
      const int total_electron = quantum_number.first;
      const double total_sz    = quantum_number.second;
      const int total_2sz      = utility::DoubleTheNumber(total_sz);
      
      if (bases_.count({total_electron, total_2sz}) != 0) {
         return;
      }
      
      std::cout << "Generating Basis..." << std::flush;
      const int max_n_up_down = static_cast<int>(total_electron/2);
      const std::int64_t dim_target = CalculateTargetDim({total_electron, total_sz});
      
      std::vector<std::int64_t>().swap(bases_[{total_electron, total_2sz}]);
      auto &basis_ref = bases_.at({total_electron, total_2sz});
      
      std::vector<std::int64_t> site_constant(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant[site] = static_cast<std::int64_t>(std::pow(dim_onsite_, site));
      }
      
      std::vector<int> basis_list(system_size_);
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<std::vector<std::int64_t>> temp_basis(num_threads);
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         const int n_up   = static_cast<int>((total_electron - 2*n_up_down + total_2sz)/2);
         const int n_down = static_cast<int>((total_electron - 2*n_up_down - total_2sz)/2);
         const int n_vac  = system_size_ - total_electron + n_up_down;
         if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
            for (int s = 0; s < n_vac; ++s) {
               basis_list[s] = 0;
            }
            for (int s = 0; s < n_up; ++s) {
               basis_list[s + n_vac] = 1;
            }
            for (int s = 0; s < n_down; ++s) {
               basis_list[s + n_vac + n_up] = 2;
            }
            for (int s = 0; s < n_up_down; ++s) {
               basis_list[s + n_vac + n_up + n_down] = 3;
            }
            
            const std::int64_t size = utility::CalculateNumCombination(basis_list);
            std::vector<std::vector<int>> temp_basis_list(num_threads);
            
#pragma omp parallel num_threads (num_threads)
            {
               const int thread_num = omp_get_thread_num();
               const std::int64_t loop_begin = thread_num*size/num_threads;
               const std::int64_t loop_end   = (thread_num + 1)*size/num_threads;
               temp_basis_list[thread_num]   = basis_list;
               utility::CalculateNthPermutation(&temp_basis_list[thread_num], loop_begin);
               
               for (std::int64_t j = loop_begin; j < loop_end; ++j) {
                  std::int64_t basis_global = 0;
                  for (std::size_t k = 0; k < temp_basis_list[thread_num].size(); ++k) {
                     basis_global += temp_basis_list[thread_num][k]*site_constant[k];
                  }
                  temp_basis[thread_num].push_back(basis_global);
                  std::next_permutation(temp_basis_list[thread_num].begin(), temp_basis_list[thread_num].end());
               }
            }
         }
      }
      for (auto &&basis: temp_basis) {
         basis_ref.insert(basis_ref.end(), basis.begin(), basis.end());
         std::vector<std::int64_t>().swap(basis);
      }
      
#else
      
      basis_ref.reserve(dim_target);
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         const int n_up   = static_cast<int>((total_electron - 2*n_up_down + total_2sz)/2);
         const int n_down = static_cast<int>((total_electron - 2*n_up_down - total_2sz)/2);
         const int n_vac  = system_size_ - total_electron + n_up_down;
         if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
            for (int s = 0; s < n_vac; ++s) {
               basis_list[s] = 0;
            }
            for (int s = 0; s < n_up; ++s) {
               basis_list[s + n_vac] = 1;
            }
            for (int s = 0; s < n_down; ++s) {
               basis_list[s + n_vac + n_up] = 2;
            }
            for (int s = 0; s < n_up_down; ++s) {
               basis_list[s + n_vac + n_up + n_down] = 3;
            }
            
            do {
               std::int64_t basis_global = 0;
               for (std::size_t j = 0; j < basis_list.size(); ++j) {
                  basis_global += basis_list[j]*site_constant[j];
               }
               basis_ref.push_back(basis_global);
            } while (std::next_permutation(basis_list.begin(), basis_list.end()));
         }
      }
#endif
      
      if (static_cast<std::int64_t>(bases_.at({total_electron, total_2sz}).size()) != dim_target) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(basis_ref.begin(), basis_ref.end());
      
      bases_inv_[{total_electron, total_2sz}].clear();
      
      auto &basis_inv_ref = bases_inv_.at({total_electron, total_2sz});
      for (std::int64_t i = 0; i < dim_target; ++i) {
         basis_inv_ref[basis_ref[i]] = i;
      }
      
      if (basis_inv_ref.size() != basis_ref.size()) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         ss << "The same basis has been detected" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1 The matrix of an onsite operator.
   //! @param m_2 The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<std::pair<int, double>> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m2;
      for (std::int64_t i = 0; i < m_1.row_dim; ++i) {
         for (std::int64_t j = m_1.row[i]; j < m_1.row[i + 1]; ++j) {
            if (m_1.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_1.col[j])));
            }
         }
      }
      for (std::int64_t i = 0; i < m_2.row_dim; ++i) {
         for (std::int64_t j = m_2.row[i]; j < m_2.row[i + 1]; ++j) {
            if (m_2.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_2.col[j])));
            }
         }
      }
      std::vector<std::pair<int, double>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            const bool c1 = isValidQNumber(del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_);
            if (del_sec_m1 == del_sec_m2 && c1) {
               target_sector_set.push_back({del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_});
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
   std::vector<std::pair<std::pair<int, double>, std::pair<int, double>>> GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_ket, const CRS &m_3_ket) const {
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m2;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m3;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_1_bra.col[j])));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_ket.row_dim; ++i) {
         for (std::int64_t j = m_2_ket.row[i]; j < m_2_ket.row[i + 1]; ++j) {
            if (m_2_ket.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_2_ket.col[j])));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_3_ket.col[j])));
            }
         }
      }
      
      std::vector<std::pair<std::pair<int, double>, std::pair<int, double>>> target_sector_set;
      
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               const std::pair<int, double> del_sec_m2_m3 = {del_sec_m2.first + del_sec_m3.first, del_sec_m2.second + del_sec_m3.second};
               const bool c1 = isValidQNumber(del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_);
               const bool c2 = isValidQNumber(del_sec_m3.first + total_electron_, del_sec_m3.second + 0.5*total_2sz_);
               if (del_sec_m1 == del_sec_m2_m3 && c1 && c2) {
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
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1_bra The matrix of an onsite operator.
   //! @param m_2_bra The matrix of an onsite operator.
   //! @param m_3_ket The matrix of an onsite operator.
   //! @param m_4_ket The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<std::tuple<std::pair<int, double>, std::pair<int, double>, std::pair<int, double>>>
   GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_bra, const CRS &m_3_ket, const CRS &m_4_ket) const {
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m2;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m3;
      std::unordered_set<std::pair<int, double>, utility::PairHash> delta_sector_set_m4;
      
      for (std::int64_t i = 0; i < m_1_bra.row_dim; ++i) {
         for (std::int64_t j = m_1_bra.row[i]; j < m_1_bra.row[i + 1]; ++j) {
            if (m_1_bra.val[j] != 0.0) {
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_1_bra.col[j])));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2_bra.row_dim; ++i) {
         for (std::int64_t j = m_2_bra.row[i]; j < m_2_bra.row[i + 1]; ++j) {
            if (m_2_bra.val[j] != 0.0) {
               delta_sector_set_m2.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_2_bra.col[j])));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_3_ket.row_dim; ++i) {
         for (std::int64_t j = m_3_ket.row[i]; j < m_3_ket.row[i + 1]; ++j) {
            if (m_3_ket.val[j] != 0.0) {
               delta_sector_set_m3.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_3_ket.col[j])));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_4_ket.row_dim; ++i) {
         for (std::int64_t j = m_4_ket.row[i]; j < m_4_ket.row[i + 1]; ++j) {
            if (m_4_ket.val[j] != 0.0) {
               delta_sector_set_m4.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_4_ket.col[j])));
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
                  const bool c1 = isValidQNumber(del_sec_m1.first    + total_electron_, del_sec_m1.second    + 0.5*total_2sz_);
                  const bool c2 = isValidQNumber(del_sec_m1_m2.first + total_electron_, del_sec_m1_m2.second + 0.5*total_2sz_);
                  const bool c3 = isValidQNumber(del_sec_m4.first    + total_electron_, del_sec_m4.second    + 0.5*total_2sz_);
                  if (del_sec_m1_m2 == del_sec_m3_m4 && c1 && c2 && c3) {
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
      
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool isValidQNumber(const int system_size, const int total_electron, const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      const bool c1 = (0 <= total_electron && total_electron <= 2*system_size);
      const bool c2 = ((total_electron - total_2sz)%2 == 0);
      const bool c3 = (-total_electron <= total_2sz && total_2sz <= total_electron);
      if (c1 && c2 && c3) {
         return true;
      }
      else {
         return false;
      }
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param system_size The system size \f$ N\f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   static std::int64_t CalculateTargetDim(const int system_size, const int total_electron, const double total_sz) {
      if (!isValidQNumber(system_size, total_electron, total_sz)) {
         return 0;
      }
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size);
      const int max_n_up_down = static_cast<int>(total_electron/2);
      std::int64_t dim = 0;
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         const int n_up   = static_cast<int>((total_electron - 2*n_up_down + total_2sz)/2);
         const int n_down = static_cast<int>((total_electron - 2*n_up_down - total_2sz)/2);
         const int n_vac  = system_size - total_electron + n_up_down;
         if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
            // TODO: Detect Overflow
            dim += binom[system_size][n_up]*binom[system_size - n_up][n_down]*binom[system_size - n_up - n_down][n_up_down];
         }
      }
      return dim;
   }
   
   //! @brief Generate the annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorCUp() {
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      const int dim_onsite = 4;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; row++) {
         for (int col = 0; col < dim_onsite; col++) {
            if ((col == 1 && row == 0) || (col == 3 && row == 2)) {
               matrix.col.push_back(col);
               matrix.val.push_back(1.0);
            }
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      
      matrix.tag = sparse_matrix::CRSTag::FERMION;
      
      return matrix;
   }
   
   //! @brief Generate the annihilation operator for the electrons with the down spin \f$ \hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{c}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorCDown() {
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      const int dim_onsite = 4;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; row++) {
         for (int col = 0; col < dim_onsite; col++) {
            if (col == 2 && row == 0) {
               matrix.col.push_back(col);
               matrix.val.push_back(1.0);
            }
            else if (col == 3 && row == 1) {
               matrix.col.push_back(col);
               matrix.val.push_back(-1.0);
            }
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      
      matrix.tag = sparse_matrix::CRSTag::FERMION;

      return matrix;
   }
   
   //! @brief Generate the creation operator for the electrons with the up spin
   //! \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorCUpDagger() {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCUp());
   }
   
   //! @brief Generate the creation operator for the electrons with the down spin
   //! \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorCDownDagger() {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCDown());
   }
   
   //! @brief Generate the number operator for the electrons with the up spin
   //! \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorNCUp() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(), 1.0, CreateOnsiteOperatorCUp());
   }
   
   //! @brief Generate the number operator for the electrons with the down spin
   //! \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorNCDown() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(), 1.0, CreateOnsiteOperatorCDown());
   }
   
   //! @brief Generate the number operator for the electrons
   //! \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}\f$.
   static CRS CreateOnsiteOperatorNC() {
      return sparse_matrix::CalculateMatrixMatrixSum(1.0, CreateOnsiteOperatorNCUp(), 1.0, CreateOnsiteOperatorNCDown());
   }
   
   //! @brief Generate the spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   static CRS CreateOnsiteOperatorSx() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorSp(), 0.5, CreateOnsiteOperatorSm());
   }
   
   //! @brief Generate the spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   static CRS CreateOnsiteOperatoriSy() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorSp(), -0.5, CreateOnsiteOperatorSm());
   }
   
   //! @brief Generate the spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   static CRS CreateOnsiteOperatorSz() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorNCUp(), -0.5,CreateOnsiteOperatorNCDown());
   }
   
   //! @brief Generate the raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   static CRS CreateOnsiteOperatorSp() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(), 1.0, CreateOnsiteOperatorCDown());
   }
   
   //! @brief Generate the lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   static CRS CreateOnsiteOperatorSm() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(), 1.0, CreateOnsiteOperatorCUp());
   }
   
   //! @brief Calculate difference of the number of total electrons and the total sz
   //! from the rows and columns in the matrix representation of an onsite operator.
   //! @param row The row in the matrix representation of an onsite operator.
   //! @param col The column in the matrix representation of an onsite operator.
   //! @return The differences of the total electron and the total sz.
   static std::pair<int, double> CalculateQuntumNumberDifference(const int row, const int col) {
      if (row == col && 0 <= row && row < 4 && 0 <= col && col < 4) {
         return {+0, +0.0};
      }
      else if (row == 0 && col == 1) {
         return {-1, -0.5};
      }
      else if (row == 0 && col == 2) {
         return {-1, +0.5};
      }
      else if (row == 0 && col == 3) {
         return {-2, +0.0};
      }
      else if (row == 1 && col == 0) {
         return {+1, +0.5};
      }
      else if (row == 1 && col == 2) {
         return {+0, +1.0};
      }
      else if (row == 1 && col == 3) {
         return {-1, +0.5};
      }
      else if (row == 2 && col == 0) {
         return {+1, -0.5};
      }
      else if (row == 2 && col == 1) {
         return {+0, -1.0};
      }
      else if (row == 2 && col == 3) {
         return {-1, -0.5};
      }
      else if (row == 3 && col == 0) {
         return {+2, +0.0};
      }
      else if (row == 3 && col == 1) {
         return {+1, -0.5};
      }
      else if (row == 3 && col == 2) {
         return {+1, +0.5};
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The dimenstion of the matrix must be 4";
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Get the system size \f$ N\f$.
   //! @return The system size \f$ N\f$.
   inline int GetSystemSize() const { return system_size_; }
   
   //! @brief Get dimension of the local Hilbert space, 4.
   //! @return The dimension of the local Hilbert space, 4.
   inline int GetDimOnsite() const { return dim_onsite_; }
   
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline double GetTotalSz() const { return 0.5*total_2sz_; }
   
   //! @brief Get the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @return The total electrons.
   inline int GetTotalElectron() const { return total_electron_; }
   
   //! @brief Get the annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}_{\uparrow}\f$.
   inline const CRS &GetOnsiteOperatorCUp() const { return onsite_operator_c_up_;   }
   
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
   inline const CRS &GetOnsiteOperatorSx () const { return onsite_operator_sx_ ; }
   
   //! @brief Get the spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   inline const CRS &GetOnsiteOperatoriSy() const { return onsite_operator_isy_; }
   
   //! @brief Get the spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   inline const CRS &GetOnsiteOperatorSz () const { return onsite_operator_sz_ ; }
   
   //! @brief Get the raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   inline const CRS &GetOnsiteOperatorSp () const { return onsite_operator_sp_ ; }
   
   //! @brief Get the lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   inline const CRS &GetOnsiteOperatorSm () const { return onsite_operator_sm_ ; }
   
   //! @brief Get calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @return calculated_eigenvector_set_.
   inline const std::unordered_set<int> &GetCalculatedEigenvectorSet() const {
      return calculated_eigenvector_set_;
   }
   
   inline const std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, compnal::utility::PairHash> &GetBases() const {
      return bases_;
   }
   
   inline const std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, compnal::utility::PairHash> &GetBasesInv() const {
      return bases_inv_;
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron
   //! \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetBasis(const std::pair<int, double> &quantum_number) const {
      return bases_.at({quantum_number.first, utility::DoubleTheNumber(quantum_number.second)});
   }
   
   //! @brief Get inverse basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron
   //! \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @return Basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetBasisInv(const std::pair<int, double> &quantum_number) const {
      return bases_inv_.at({quantum_number.first, utility::DoubleTheNumber(quantum_number.second)});
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetTargetBasis() const {
      return bases_.at({total_electron_, total_2sz_});
   }
   
   //! @brief Get inverse basis of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Inverse basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetTargetBasisInv() const {
      return bases_inv_.at({total_electron_, total_2sz_});
   }
   
protected:
   //! @brief The annihilation operator for the electrons with the up spin \f$ \hat{c}_{\uparrow}\f$.
   CRS onsite_operator_c_up_;
   
   //! @brief The annihilation operator for the electrons with the down spin \f$ \hat{c}_{\downarrow}\f$.
   CRS onsite_operator_c_down_;
   
   //! @brief The creation operator for the electrons with the up spin.
   //! \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   CRS onsite_operator_c_up_dagger_;
   
   //! @brief The creation operator for the electrons with the down spin.
   //! \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   CRS onsite_operator_c_down_dagger_;
   
   //! @brief The number operator for the electrons with the up spin
   //! \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   CRS onsite_operator_nc_up_;
   
   //! @brief The number operator for the electrons with the down spin
   //! \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   CRS onsite_operator_nc_down_;
   
   //! @brief The number operator for the electrons
   //! \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   CRS onsite_operator_nc_;
   
   //! @brief The spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   CRS onsite_operator_sx_;
   
   //! @brief The spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   CRS onsite_operator_isy_;
   
   //! @brief The spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   CRS onsite_operator_sz_;
   
   //! @brief The raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   CRS onsite_operator_sp_;
   
   //! @brief The lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   CRS onsite_operator_sm_;
   
   //! @brief The system size
   int system_size_ = 0;
   
   //! @brief Twice the number of the total sz \f$ 2\langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   int total_2sz_ = 0;
   
   //! @brief The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   int total_electron_ = 0;
   
   //! @brief The dimension of the local Hilbert space, 4.
   const int dim_onsite_ = 4;
   
   //! @brief The calculated eigenvectors and eigenvalues.
   std::unordered_set<int> calculated_eigenvector_set_;
   
   //! @brief Bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, utility::PairHash> bases_;
   
   //! @brief Inverse bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, utility::PairHash> bases_inv_;
   
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_c_up_   = CreateOnsiteOperatorCUp();
      onsite_operator_c_down_ = CreateOnsiteOperatorCDown();
      onsite_operator_c_up_dagger_   = CreateOnsiteOperatorCUpDagger();
      onsite_operator_c_down_dagger_ = CreateOnsiteOperatorCDownDagger();
      onsite_operator_nc_up_   = CreateOnsiteOperatorNCUp();
      onsite_operator_nc_down_ = CreateOnsiteOperatorNCDown();
      onsite_operator_nc_ = CreateOnsiteOperatorNC();
      onsite_operator_sx_ = CreateOnsiteOperatorSx ();
      onsite_operator_isy_= CreateOnsiteOperatoriSy();
      onsite_operator_sz_ = CreateOnsiteOperatorSz ();
      onsite_operator_sp_ = CreateOnsiteOperatorSp ();
      onsite_operator_sm_ = CreateOnsiteOperatorSm ();
   }
   
};

}
}


#endif /* COMPNAL_MODEL_BASE_U1ELECTRON_1D_HPP_ */
