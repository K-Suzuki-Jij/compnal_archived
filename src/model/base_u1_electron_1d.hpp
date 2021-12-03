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

template<typename RealType>
class BaseU1Electron_1D {
   
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   using ValueType = RealType;
   
   BaseU1Electron_1D() {
      SetOnsiteOperator();
   }
   
   explicit BaseU1Electron_1D(const int system_size): BaseU1Electron_1D() {
      SetSystemSize(system_size);
   }
   
   BaseU1Electron_1D(const int system_size, const int total_electron): BaseU1Electron_1D(system_size) {
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
      if (!isValidQNumber({total_electron_, total_sz})) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "There is no target space specified by total_sz = " << total_sz << std::endl;
         throw std::runtime_error(ss.str());
      }
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetTotalElectron(const int total_electron) {
      if (!isValidQNumber({total_electron, 0.5*total_2sz_})) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "There is no target space specified by total_sz = " << total_2sz_ << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (total_electron_ != total_electron) {
         total_electron_ = total_electron;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetCalculatedEigenvectorSet(const std::int64_t level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   bool isValidQNumber(const std::pair<int, double> &quantum_number) const {
      const int total_electron = quantum_number.first;
      const int total_2sz      = utility::DoubleTheNumber(quantum_number.second);
      const bool c1 = (0 <= total_electron && total_electron <= 2*system_size_);
      const bool c2 = ((total_electron - total_2sz)%2 == 0);
      const bool c3 = (-total_electron <= total_2sz && total_2sz <= total_electron);
      if (c1 && c2 && c3) {
         return true;
      }
      else {
         return false;
      }
   }
   
   void PrintBasisOnsite() const {
      std::cout << "row " << 0 << ": |vac>" << std::endl;
      std::cout << "row " << 1 << ": |↑>"   << std::endl;
      std::cout << "row " << 2 << ": |↓>"   << std::endl;
      std::cout << "row " << 3 << ": |↑↓>"  << std::endl;
   }
   
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim({total_electron_, 0.5*total_2sz_});
   }
   
   std::int64_t CalculateTargetDim(const int total_electron, const double total_sz) const {
      return CalculateTargetDim({total_electron, total_sz});
   }
   
   std::int64_t CalculateTargetDim(const std::pair<int, double> &quantum_number) const {
      if (!isValidQNumber(quantum_number)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const int total_electron = quantum_number.first;
      const int total_2sz      = utility::DoubleTheNumber(quantum_number.second);
      const std::vector<std::vector<std::int64_t>> binom = utility::CalculateBinomialTable(system_size_);
      const int max_n_up_down = static_cast<int>(total_electron/2);
      std::int64_t dim = 0;
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         const int n_up   = static_cast<int>((total_electron - 2*n_up_down + total_2sz)/2);
         const int n_down = static_cast<int>((total_electron - 2*n_up_down - total_2sz)/2);
         const int n_vac  = system_size_ - total_electron + n_up_down;
         if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
            // TODO: Detect Overflow
            dim += binom[system_size_][n_up]*binom[system_size_ - n_up][n_down]*binom[system_size_ - n_up - n_down][n_up_down];
         }
      }
      return dim;
   }
   
   void GenerateBasis() {
      GenerateBasis({total_electron_, 0.5*total_2sz_});
   }
   
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
      
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
   }
   
   std::vector<std::pair<int, double>> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m1;
      std::unordered_set<std::pair<int, double>, utility::pair_hash> delta_sector_set_m2;
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
               delta_sector_set_m1.emplace(CalculateQuntumNumberDifference(static_cast<int>(i), static_cast<int>(m_2.col[j])));
            }
         }
      }
      std::vector<std::pair<int, double>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            if (del_sec_m1 == del_sec_m2) {
               target_sector_set.push_back(std::pair<int, double>{del_sec_m1.first + total_electron_, del_sec_m1.second + 0.5*total_2sz_});
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
      return matrix;
   }
   
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
      return matrix;
   }
   
   static CRS CreateOnsiteOperatorCUpDagger() {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCUp());
   }
   
   static CRS CreateOnsiteOperatorCDownDagger() {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCDown());
   }
   
   static CRS CreateOnsiteOperatorNCUp() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(), 1.0, CreateOnsiteOperatorCUp());
   }
   
   static CRS CreateOnsiteOperatorNCDown() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(), 1.0, CreateOnsiteOperatorCDown());
   }
   
   static CRS CreateOnsiteOperatorNC() {
      return sparse_matrix::CalculateMatrixMatrixSum(1.0, CreateOnsiteOperatorNCUp(), 1.0, CreateOnsiteOperatorNCUp());
   }
   
   static CRS CreateOnsiteOperatorSx() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorSp(), 0.5, CreateOnsiteOperatorSm());
   }
   
   static CRS CreateOnsiteOperatoriSy() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorSp(), -0.5, CreateOnsiteOperatorSm());
   }
   
   static CRS CreateOnsiteOperatorSz() {
      return sparse_matrix::CalculateMatrixMatrixSum(0.5, CreateOnsiteOperatorNCUp(), -0.5,CreateOnsiteOperatorNCDown());
   }
   
   static CRS CreateOnsiteOperatorSp() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCUpDagger(), 1.0, CreateOnsiteOperatorCDown());
   }
   
   static CRS CreateOnsiteOperatorSm() {
      return sparse_matrix::CalculateMatrixMatrixProduct(1.0, CreateOnsiteOperatorCDownDagger(), 1.0, CreateOnsiteOperatorCUp());
   }
   
   inline int    GetSystemSize()           const { return system_size_;            }
   inline int    GetDimOnsite()            const { return dim_onsite_;             }
   inline int    GetTotal2Sz()             const { return total_2sz_;              }
   inline double GetTotalSz()              const { return 0.5*total_2sz_;          }
   inline int    GetNumConservedQuantity() const { return num_conserved_quantity_; }
   inline int    GetTotalElectron()        const { return total_electron_; }
   
   inline const CRS &GetOnsiteOperatorCUp()         const { return onsite_operator_c_up_;   }
   inline const CRS &GetOnsiteOperatorCDown()       const { return onsite_operator_c_down_; }
   inline const CRS &GetOnsiteOperatorCUpDagger()   const { return onsite_operator_c_up_dagger_; }
   inline const CRS &GetOnsiteOperatorCDownDagger() const { return onsite_operator_c_down_dagger_; }
   inline const CRS &GetOnsiteOperatorNCUp()        const { return onsite_operator_nc_up_; }
   inline const CRS &GetOnsiteOperatorNCDown()      const { return onsite_operator_nc_down_; }
   inline const CRS &GetOnsiteOperatorNC()          const { return onsite_operator_nc_; }
   inline const CRS &GetOnsiteOperatorSx ()         const { return onsite_operator_sx_ ; }
   inline const CRS &GetOnsiteOperatoriSy()         const { return onsite_operator_isy_; }
   inline const CRS &GetOnsiteOperatorSz ()         const { return onsite_operator_sz_ ; }
   inline const CRS &GetOnsiteOperatorSp ()         const { return onsite_operator_sp_ ; }
   inline const CRS &GetOnsiteOperatorSm ()         const { return onsite_operator_sm_ ; }
   
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
   CRS onsite_operator_sx_;
   CRS onsite_operator_isy_;
   CRS onsite_operator_sz_;
   CRS onsite_operator_sp_;
   CRS onsite_operator_sm_;
   
   int system_size_    = 0;
   int total_2sz_      = 0;
   int total_electron_ = 0;
   
   const int dim_onsite_ = 4;
   
   const int num_conserved_quantity_ = 2;
   
   std::unordered_set<int> calculated_eigenvector_set_;
   
   std::unordered_map<std::pair<int, int>, std::vector<std::int64_t>, utility::pair_hash> bases_;
   std::unordered_map<std::pair<int, int>, std::unordered_map<std::int64_t, std::int64_t>, utility::pair_hash> bases_inv_;
   
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
   
   std::pair<int, double> CalculateQuntumNumberDifference(const int row, const int col) const {
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
   
   
};

}
}


#endif /* COMPNAL_MODEL_BASE_U1ELECTRON_1D_HPP_ */
