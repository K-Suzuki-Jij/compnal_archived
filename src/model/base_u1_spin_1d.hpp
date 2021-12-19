//
//  base.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/18.
//

#ifndef COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_
#define COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"

#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

//! @brief The base class for one-dimensional spin systems with the U(1) symmetry.
//! @tparam RealType The type of real values.
template<typename RealType>
class BaseU1Spin_1D {
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   //! @brief The type of real values.
   using ValueType = RealType;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of BaseU1Spin_1D class.
   BaseU1Spin_1D() {
      SetOnsiteOperator();
   }
   
   //! @brief Constructor of BaseU1Spin_1D class.
   //! @param system_size The system size \f$ N \f$.
   explicit BaseU1Spin_1D(const int system_size): BaseU1Spin_1D() {
      SetSystemSize(system_size);
   }
   
   //! @brief Constructor of BaseU1Spin_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   BaseU1Spin_1D(const int system_size, const double magnitude_spin): BaseU1Spin_1D(system_size) {
      SetMagnitudeSpin(magnitude_spin);
   }
   
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
   
   //! @brief Set the magnitude of the spin \f$ S \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   void SetMagnitudeSpin(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      if (magnitude_2spin <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Please set magnitude_2spin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_2spin_ != magnitude_2spin) {
         magnitude_2spin_ = magnitude_2spin;
         dim_onsite_      = magnitude_2spin + 1;
         SetOnsiteOperator();
         bases_.clear();
         bases_inv_.clear();
         calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set target Hilbert space specified by the total sz to be diagonalized.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$.
   void SetTotalSz(const double total_sz) {
      if (isValidQNumber(total_sz) == false) {
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
   
   //! @brief Set calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @param level Energy level.
   void SetCalculatedEigenvectorSet(const std::int64_t level) {
      calculated_eigenvector_set_.emplace(level);
   }
      
   //! @brief Check if there is a subspace specified by the input total sz.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const double total_sz) const {
      return isValidQNumber(system_size_, 0.5*magnitude_2spin_, total_sz);
   }
   
   //! @brief Print the onsite bases.
   void PrintBasisOnsite() const {
      const double magnitude_spin = magnitude_2spin_/2.0;
      for (int row = 0; row < dim_onsite_; ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin - row << ">" << std::endl;
      }
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
      return CalculateTargetDim(system_size_, 0.5*magnitude_2spin_, total_sz);
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   void GenerateBasis() {
      GenerateBasis(0.5*total_2sz_);
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   void GenerateBasis(const double total_sz) {
      if (!isValidQNumber(total_sz)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const auto start = std::chrono::system_clock::now();
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      
      if (bases_.count(total_2sz) != 0) {
         return;
      }
      
      std::cout << "Generating Basis..." << std::flush;
      
      const int shifted_2sz = (system_size_*magnitude_2spin_ - total_2sz)/2;
      const std::int64_t dim_target = CalculateTargetDim(total_sz);
      std::vector<std::vector<int>> partition_integers;
      utility::GenerateIntegerPartition(&partition_integers, shifted_2sz, magnitude_2spin_);
      
      std::vector<std::int64_t> site_constant(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant[site] = static_cast<std::int64_t>(std::pow(dim_onsite_, site));
      }
      
      if (bases_.count(total_2sz) == 0) {
         bases_[total_2sz] = std::vector<std::int64_t>();
      }
      
      std::vector<std::int64_t>().swap(bases_.at(total_2sz));
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<std::vector<std::int64_t>> temp_basis(num_threads);
      for (auto &&integer_list: partition_integers) {
         const bool condition1 = (0 < integer_list.size()) && (static_cast<int>(integer_list.size()) <= system_size_);
         const bool condition2 = (integer_list.size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            for (int j = static_cast<int>(integer_list.size()); j < system_size_; ++j) {
               integer_list.push_back(0);
            }
            
            const std::int64_t size = utility::CalculateNumCombination(integer_list);
            std::vector<std::vector<int>> temp_partition_integer(num_threads);
            
#pragma omp parallel num_threads (num_threads)
            {
               const int thread_num = omp_get_thread_num();
               const std::int64_t loop_begin = thread_num*size/num_threads;
               const std::int64_t loop_end   = (thread_num + 1)*size/num_threads;
               temp_partition_integer[thread_num] = integer_list;
               utility::CalculateNthPermutation(&temp_partition_integer[thread_num], loop_begin);
               
               for (std::int64_t j = loop_begin; j < loop_end; ++j) {
                  std::int64_t basis_global = 0;
                  const auto iter_begin = temp_partition_integer[thread_num].begin();
                  const auto iter_end   = temp_partition_integer[thread_num].end();
                  for (auto itr = iter_begin; itr != iter_end; ++itr) {
                     basis_global += *itr*site_constant[std::distance(iter_begin, itr)];
                  }
                  temp_basis[thread_num].push_back(basis_global);
                  std::next_permutation(temp_partition_integer[thread_num].begin(), temp_partition_integer[thread_num].end());
               }
            }
         }
      }
      for (auto &&basis: temp_basis) {
         bases_.at(total_2sz).insert(bases_.at(total_2sz).end(), basis.begin(), basis.end());
         std::vector<std::int64_t>().swap(basis);
      }
      
#else
      bases_.at(total_2sz).reserve(dim_target);
      
      for (auto &&integer_list: partition_integers) {
         const bool condition1 = (0 < integer_list.size()) && (static_cast<int>(integer_list.size()) <= system_size_);
         const bool condition2 = (integer_list.size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            
            for (std::int64_t j = integer_list.size(); j < system_size_; ++j) {
               integer_list.push_back(0);
            }
            
            std::sort(integer_list.begin(), integer_list.end());
            
            do {
               std::int64_t basis_global = 0;
               for (std::size_t j = 0; j < integer_list.size(); ++j) {
                  basis_global += integer_list[j]*site_constant[j];
               }
               bases_.at(total_2sz).push_back(basis_global);
            } while (std::next_permutation(integer_list.begin(), integer_list.end()));
         }
      }
      
#endif
      
      if (static_cast<std::int64_t>(bases_.at(total_2sz).size()) != dim_target) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(bases_.at(total_2sz).begin(), bases_.at(total_2sz).end());
      
      if (bases_inv_.count(total_2sz) == 0) {
         bases_inv_[total_2sz] = std::unordered_map<std::int64_t, std::int64_t>();
      }
      bases_inv_.at(total_2sz).clear();
      
      for (std::int64_t i = 0; i < dim_target; ++i) {
         bases_inv_.at(total_2sz)[bases_.at(total_2sz)[i]] = i;
      }
      
      const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
      std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1 The matrix of an onsite operator.
   //! @param m_2 The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<double> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      
      std::unordered_set<double> delta_sector_set_m1;
      std::unordered_set<double> delta_sector_set_m2;
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
      std::vector<double> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            const bool c1 = isValidQNumber(del_sec_m1 + 0.5*total_2sz_);
            if (del_sec_m1 == del_sec_m2 && c1) {
               target_sector_set.push_back(del_sec_m1 + 0.5*total_2sz_);
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
   std::vector<std::pair<double, double>> GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_ket, const CRS &m_3_ket) const {
      std::unordered_set<double> delta_sector_set_m1;
      std::unordered_set<double> delta_sector_set_m2;
      std::unordered_set<double> delta_sector_set_m3;
      
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
      
      std::vector<std::pair<double, double>> target_sector_set;
      
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               const bool c1 = isValidQNumber(del_sec_m1 + 0.5*total_2sz_);
               const bool c2 = isValidQNumber(del_sec_m3 + 0.5*total_2sz_);
               if (del_sec_m1 == del_sec_m2 + del_sec_m3 && c1 && c2) {
                  target_sector_set.push_back({
                     del_sec_m1 + 0.5*total_2sz_,
                     del_sec_m3 + 0.5*total_2sz_
                  });
               }
            }
         }
      }
      return target_sector_set;
   }
   
   //! @brief Calculate the quantum numbers of excited states that appear when calculating the correlation functions.
   //! @param m_1_bra The matrix of an onsite operator.
   //! @param m_2_bra The matrix of an onsite operator.
   //! @param m_3_ket The matrix of an onsite operator.
   //! @param m_4_ket The matrix of an onsite operator.
   //! @return The list of quantum numbers.
   std::vector<std::tuple<double, double, double>> GenerateTargetSector(const CRS &m_1_bra, const CRS &m_2_bra, const CRS &m_3_ket, const CRS &m_4_ket) const {
      std::unordered_set<double> delta_sector_set_m1;
      std::unordered_set<double> delta_sector_set_m2;
      std::unordered_set<double> delta_sector_set_m3;
      std::unordered_set<double> delta_sector_set_m4;
      
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
      
      std::vector<std::tuple<double, double, double>> target_sector_set;
      for (const auto &del_sec_m1: delta_sector_set_m1) {
         for (const auto &del_sec_m2: delta_sector_set_m2) {
            for (const auto &del_sec_m3: delta_sector_set_m3) {
               for (const auto &del_sec_m4: delta_sector_set_m4) {
                  const bool c1 = isValidQNumber(del_sec_m1 + 0.5*total_2sz_);
                  const bool c2 = isValidQNumber(del_sec_m1 + del_sec_m2 + 0.5*total_2sz_);
                  const bool c3 = isValidQNumber(del_sec_m4 + 0.5*total_2sz_);
                  if (del_sec_m1 + del_sec_m2 == del_sec_m3 + del_sec_m4 && c1 && c2 && c3) {
                     target_sector_set.push_back({
                        del_sec_m1 + 0.5*total_2sz_,
                        del_sec_m1 + del_sec_m2 + 0.5*total_2sz_,
                        del_sec_m4 + 0.5*total_2sz_
                     });
                  }
               }
            }
         }
      }
      return target_sector_set;
   }
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool isValidQNumber(const int system_size, const double magnitude_spin, const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const bool c1 = ((system_size*magnitude_2spin - total_2sz)%2 == 0);
      const bool c2 = (-system_size*magnitude_2spin <= total_2sz);
      const bool c3 = (total_2sz <= system_size*magnitude_2spin);
      if (c1 && c2 && c3) {
         return true;
      }
      else {
         return false;
      }
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   static std::int64_t CalculateTargetDim(const int system_size, const double magnitude_spin, const double total_sz) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      if (!isValidQNumber(system_size, magnitude_spin, total_sz)) {
         return 0;
      }
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      const int max_total_2sz = system_size*magnitude_2spin;
      std::vector<std::vector<std::int64_t>> dim(system_size, std::vector<std::int64_t>(max_total_2sz + 1));
      for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
         dim[0][(s + magnitude_2spin)/2] = 1;
      }
      for (int site = 1; site < system_size; site++) {
         for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
            for (int s_prev = -magnitude_2spin*site; s_prev <= magnitude_2spin*site; s_prev += 2) {
               const std::int64_t a = dim[site    ][(s + s_prev + magnitude_2spin*(site + 1))/2];
               const std::int64_t b = dim[site - 1][(s_prev + magnitude_2spin*site)/2];
               if (a >= INT64_MAX - b) {
                  throw std::runtime_error("Overflow detected for sumation using uint64_t");
               }
               dim[site][(s + s_prev + magnitude_2spin*(site + 1))/2] = a + b;
            }
         }
      }
      return dim[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator for the x-direction \f$ \hat{s}^{x}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   static CRS CreateOnsiteOperatorSx(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      int a = 0;
      int b = 1;
      
      matrix.val.push_back(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
      matrix.col.push_back(b);
      matrix.row[1] = matrix.col.size();
      
      for (int row = 1; row < dim_onsite - 1; ++row) {
         a = row;
         b = row - 1;
         matrix.val.push_back(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
         matrix.col.push_back(b);
         
         a = row;
         b = row + 1;
         matrix.val.push_back(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
         matrix.col.push_back(b);
         matrix.row[row + 1] = matrix.col.size();
      }
      
      a = dim_onsite - 1;
      b = dim_onsite - 2;
      
      matrix.val.push_back(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
      matrix.col.push_back(b);
      matrix.row[dim_onsite] = matrix.col.size();
      
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator for the y-direction \f$ i\hat{s}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   static CRS CreateOnsiteOperatoriSy(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      int a = 0;
      int b = 1;
      
      matrix.val.push_back(0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
      matrix.col.push_back(b);
      matrix.row[1] = matrix.col.size();
      
      for (int row = 1; row < dim_onsite - 1; ++row) {
         a = row;
         b = row - 1;
         matrix.val.push_back(-0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
         matrix.col.push_back(b);
         
         a = row;
         b = row + 1;
         matrix.val.push_back(0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
         matrix.col.push_back(b);
         
         matrix.row[row + 1] = matrix.col.size();
      }
      
      a = dim_onsite - 1;
      b = dim_onsite - 2;
      
      matrix.val.push_back(-0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
      matrix.col.push_back(b);
      matrix.row[dim_onsite] = matrix.col.size();
      
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator for the z-direction \f$ \hat{s}^{z}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   static CRS CreateOnsiteOperatorSz(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         const RealType val = magnitude_spin - row;
         if (val != 0.0) {
            matrix.val.push_back(val);
            matrix.col.push_back(row);
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator \f$ \hat{s}^{+}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   static CRS CreateOnsiteOperatorSp(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 1; row < dim_onsite; ++row) {
         matrix.val.push_back(std::sqrt((magnitude_spin + 1)*2*row - row*(row + 1)));
         matrix.col.push_back(row);
         matrix.row[row] = matrix.col.size();
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator \f$ \hat{s}^{-}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   static CRS CreateOnsiteOperatorSm(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 1; row < dim_onsite; ++row) {
         matrix.val.push_back(std::sqrt((magnitude_spin + 1)*2*row - row*(row + 1)));
         matrix.col.push_back(row - 1);
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   //! @brief Calculate difference of the total sz from the rows and columns in the matrix representation of an onsite operator.
   //! @param row The row in the matrix representation of an onsite operator.
   //! @param col The column in the matrix representation of an onsite operator.
   //! @return The differences of the total sz.
   inline static int CalculateQuntumNumberDifference(const int row, const int col) {
      return col - row;
   }

   //! @brief Get the system size \f$ N\f$.
   //! @return The system size \f$ N\f$.
   inline int GetSystemSize() const { return system_size_; }
   
   //! @brief Get dimension of the local Hilbert space, \f$ 2S+1\f$.
   //! @return The dimension of the local Hilbert space, \f$ 2S+1\f$.
   inline int GetDimOnsite() const { return dim_onsite_; }
      
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline double GetTotalSz() const { return 0.5*total_2sz_; }
   
   //! @brief Get the magnitude of the spin \f$ S\f$.
   //! @return The magnitude of the spin \f$ S\f$.
   inline double GetMagnitudeSpin() const { return 0.5*magnitude_2spin_; }
   
   //! @brief Get the spin-\f$ S\f$ operator for the x-direction \f$ \hat{s}^{x}\f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   inline const CRS &GetOnsiteOperatorSx () const { return onsite_operator_sx_; }
   
   //! @brief Get the spin-\f$ S\f$ operator for the y-direction \f$ i\hat{s}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   inline const CRS &GetOnsiteOperatoriSy() const { return onsite_operator_isy_; }
   
   //! @brief Get the spin-\f$ S\f$ operator for the z-direction \f$ \hat{s}^{z}\f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   inline const CRS &GetOnsiteOperatorSz () const { return onsite_operator_sz_; }
   
   //! @brief Get the spin-\f$ S\f$ raising operator \f$ \hat{s}^{+}\f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   inline const CRS &GetOnsiteOperatorSp () const { return onsite_operator_sp_; }
   
   //! @brief Get the spin-\f$ S\f$ lowering operator \f$ \hat{s}^{-}\f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   inline const CRS &GetOnsiteOperatorSm () const { return onsite_operator_sm_; }
   
   //! @brief Get calculated_eigenvector_set_, which represents the calculated eigenvectors and eigenvalues.
   //! @return calculated_eigenvector_set_.
   inline const std::unordered_set<int> &GetCalculatedEigenvectorSet() const {
      return calculated_eigenvector_set_;
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetBasis(const double total_sz) const {
      return bases_.at(utility::DoubleTheNumber(total_sz));
   }
   
   //! @brief Get inverse basis of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Inverse basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetBasisInv(const double total_sz) const {
      return bases_inv_.at(utility::DoubleTheNumber(total_sz));
   }
   
   //! @brief Get basis of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Basis.
   inline const std::vector<std::int64_t> &GetTargetBasis() const {
      return bases_.at(total_2sz_);
   }
   
   //! @brief Get inverse basis of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return Inverse basis.
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetTargetBasisInv() const {
      return bases_inv_.at(total_2sz_);
   }
   
protected:
   
   //! @brief The spin-\f$ S\f$ operator for the x-direction \f$ \hat{s}^{x}\f$.
   CRS onsite_operator_sx_;
   
   //! @brief The spin-\f$ S\f$ operator for the y-direction \f$ i\hat{s}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   CRS onsite_operator_isy_;
   
   //! @brief The spin-\f$ S\f$ operator for the z-direction \f$ \hat{s}^{z}\f$.
   CRS onsite_operator_sz_;
   
   //! @brief The spin-\f$ S\f$ raising operator \f$ \hat{s}^{+}\f$.
   CRS onsite_operator_sp_;
   
   //! @brief The the spin-\f$ S\f$ raising operator \f$ \hat{s}^{-}\f$.
   CRS onsite_operator_sm_;
   
   //! @brief The system size.
   int system_size_ = 0;
   
   //! @brief Twice the number of the total sz \f$ 2\langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   int total_2sz_ = 0;
   
   //! @brief The dimension of the local Hilbert space, \f$ 2S + 1\f$.
   int dim_onsite_ = 2;
   
   //! @brief The magnitude of the spin \f$ S\f$.
   int magnitude_2spin_ = 1;
      
   //! @brief The calculated eigenvectors and eigenvalues.
   std::unordered_set<int> calculated_eigenvector_set_;
   
   //! @brief Bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   std::unordered_map<int, std::vector<std::int64_t>> bases_;
   
   //! @brief Inverse bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>> bases_inv_;
   
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_sx_  = CreateOnsiteOperatorSx (0.5*magnitude_2spin_);
      onsite_operator_isy_ = CreateOnsiteOperatoriSy(0.5*magnitude_2spin_);
      onsite_operator_sz_  = CreateOnsiteOperatorSz (0.5*magnitude_2spin_);
      onsite_operator_sp_  = CreateOnsiteOperatorSp (0.5*magnitude_2spin_);
      onsite_operator_sm_  = CreateOnsiteOperatorSm (0.5*magnitude_2spin_);
   }
   
};



}
}


#endif /* COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_ */
