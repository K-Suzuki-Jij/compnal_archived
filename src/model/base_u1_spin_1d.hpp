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
//  Created by Kohei Suzuki on 2021/11/18.
//

#ifndef COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_
#define COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"
#include "../type.hpp"

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
   
   //! @brief Alias of quantum number (total sz) type.
   using QType = HalfInt;
   
   //! @brief Alias of quantum number hash.
   using QHash = utility::HalfIntHash;
   
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
   BaseU1Spin_1D(const int system_size, const HalfInt magnitude_spin): BaseU1Spin_1D(system_size) {
      SetMagnitudeSpin(magnitude_spin);
   }

   //! @brief Constructor of BaseU1Spin_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$.
   BaseU1Spin_1D(const int system_size,
                 const HalfInt magnitude_spin,
                 const HalfInt total_sz): BaseU1Spin_1D(system_size, magnitude_spin) {
      SetTotalSz(total_sz);
   }
   
   //! @brief Set system size.
   //! @param system_size The system size \f$ N \f$.
   void SetSystemSize(const int system_size) {
      if (system_size < 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "system_size must be a non-negative integer" << std::endl;
         ss << "system_size=" << system_size << "is not allowed" << std::endl;
         throw std::runtime_error(ss.str());
      }
      system_size_ = system_size;
   }
   
   //! @brief Set the magnitude of the spin \f$ S \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   void SetMagnitudeSpin(const HalfInt magnitude_spin) {
      if (magnitude_spin <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Please set magnitude_spin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_spin_ != magnitude_spin) {
         magnitude_spin_ = magnitude_spin;
         dim_onsite_     = 2*magnitude_spin + 1;
         SetOnsiteOperator();
      }
   }
   
   //! @brief Set target Hilbert space specified by the total sz to be diagonalized.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$.
   void SetTotalSz(const HalfInt total_sz) {
      total_sz_ = total_sz;
   }
         
   //! @brief Check if there is a subspace specified by the input total sz.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$
   //! @return ture if there exists corresponding subspace, otherwise false.
   bool isValidQNumber(const HalfInt total_sz) const {
      return isValidQNumber(system_size_, magnitude_spin_, total_sz);
   }
   
   //! @brief Calculate the number of electrons from the input onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The number of electrons.
   int CalculateNumElectron(const int basis_onsite) const {
      if (0 <= basis_onsite && basis_onsite < dim_onsite_) {
         return 0;
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__  << std::endl;
         ss << "Invalid onsite basis" << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Print the onsite bases.
   void PrintBasisOnsite() const {
      for (int row = 0; row < dim_onsite_; ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin_ - row << ">" << std::endl;
      }
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim(total_sz_);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const HalfInt total_sz) const {
      return CalculateTargetDim(system_size_, magnitude_spin_, total_sz);
   }
   
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$ and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   std::vector<std::int64_t> GenerateBasis(const HalfInt total_sz, const bool flag_display_info = false) const {
      if (!isValidQNumber(total_sz)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
   
      const auto start = std::chrono::system_clock::now();
   
      std::vector<std::int64_t> basis;
      
      if (flag_display_info) {
         std::cout << "Generating Basis..." << std::flush;
      }
      
      const int shifted_2sz = static_cast<int>(2*(system_size_*magnitude_spin_ - total_sz));
      const std::int64_t dim_target = CalculateTargetDim(total_sz);
      std::vector<std::vector<int>> partition_integers;
      utility::GenerateIntegerPartition(&partition_integers, shifted_2sz, magnitude_spin_);
      
      std::vector<std::int64_t> site_constant(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant[site] = static_cast<std::int64_t>(std::pow(dim_onsite_, site));
      }
      
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
      
      for (auto &&it: temp_basis) {
         basis.insert(basis.end(), it.begin(), it.end());
         std::vector<std::int64_t>().swap(it);
      }
      
#else
      basis.reserve(dim_target);
      
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
               basis.push_back(basis_global);
            } while (std::next_permutation(integer_list.begin(), integer_list.end()));
         }
      }
      
#endif
      
      if (static_cast<std::int64_t>(basis.size()) != dim_target) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(basis.begin(), basis.end());
      
      if (flag_display_info) {
         const auto   time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
         const double time_sec   = static_cast<double>(time_count)/sparse_matrix::TIME_UNIT_CONSTANT;
         std::cout << "\rElapsed time of generating basis:" << time_sec << "[sec]" << std::endl;
      }
      return basis;
   }
      
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool isValidQNumber(const int system_size, const HalfInt magnitude_spin, const HalfInt total_sz) {
      const int total_2sz       = 2*total_sz;
      const int magnitude_2spin = 2*magnitude_spin;
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
   static std::int64_t CalculateTargetDim(const int system_size, const HalfInt magnitude_spin, const HalfInt total_sz) {
      if (!isValidQNumber(system_size, magnitude_spin, total_sz)) {
         return 0;
      }
      if (system_size <= 0) {
         return 0;
      }
      const int magnitude_2spin = 2*magnitude_spin;
      const int total_2sz       = 2*total_sz;
      const int max_total_2sz   = system_size*magnitude_2spin;
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
   static CRS CreateOnsiteOperatorSx(const HalfInt magnitude_spin) {
      return 0.5*(CreateOnsiteOperatorSp(magnitude_spin) + CreateOnsiteOperatorSm(magnitude_spin));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator for the y-direction \f$ i\hat{s}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   static CRS CreateOnsiteOperatoriSy(const HalfInt magnitude_spin) {
      return 0.5*(CreateOnsiteOperatorSp(magnitude_spin) - CreateOnsiteOperatorSm(magnitude_spin));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator for the z-direction \f$ \hat{s}^{z}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   static CRS CreateOnsiteOperatorSz(const HalfInt magnitude_spin) {
      const int dim_onsite = 2*magnitude_spin + 1;
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
   static CRS CreateOnsiteOperatorSp(const HalfInt magnitude_spin) {
      const int dim_onsite = 2*magnitude_spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 1; row < dim_onsite; ++row) {
         matrix.val.push_back(std::sqrt((magnitude_spin + 1)*2.0*row - row*(row + 1)));
         matrix.col.push_back(row);
         matrix.row[row] = matrix.col.size();
      }
      matrix.row[dim_onsite] = matrix.col.size();
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator \f$ \hat{s}^{-}\f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   static CRS CreateOnsiteOperatorSm(const HalfInt magnitude_spin) {
      const int dim_onsite = 2*magnitude_spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 1; row < dim_onsite; ++row) {
         matrix.val.push_back(std::sqrt((magnitude_spin + 1)*2.0*row - row*(row + 1)));
         matrix.col.push_back(row - 1);
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   //! @brief Calculate difference of the total sz from the rows and columns in the matrix representation of an onsite operator.
   //! @param row The row in the matrix representation of an onsite operator.
   //! @param col The column in the matrix representation of an onsite operator.
   //! @return The differences of the total sz.
   template<typename IntegerType>
   inline HalfInt CalculateQNumber(const IntegerType row, const IntegerType col) {
      return static_cast<HalfInt>(col - row + total_sz_);
   }

   //! @brief Get the system size \f$ N\f$.
   //! @return The system size \f$ N\f$.
   inline int GetSystemSize() const { return system_size_; }
   
   //! @brief Get dimension of the local Hilbert space, \f$ 2S+1\f$.
   //! @return The dimension of the local Hilbert space, \f$ 2S+1\f$.
   inline int GetDimOnsite() const { return dim_onsite_; }
      
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline HalfInt GetTotalSz() const { return total_sz_; }
   
   inline QType GetQNumber() const {return total_sz_; }
   
   //! @brief Get the magnitude of the spin \f$ S\f$.
   //! @return The magnitude of the spin \f$ S\f$.
   inline HalfInt GetMagnitudeSpin() const { return magnitude_spin_; }
   
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
   HalfInt total_sz_ = 0_hi;
   
   //! @brief The dimension of the local Hilbert space, \f$ 2S + 1\f$.
   int dim_onsite_ = 2;
   
   //! @brief The magnitude of the spin \f$ S\f$.
   HalfInt magnitude_spin_ = 0.5_hi;
            
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_sx_  = CreateOnsiteOperatorSx (magnitude_spin_);
      onsite_operator_isy_ = CreateOnsiteOperatoriSy(magnitude_spin_);
      onsite_operator_sz_  = CreateOnsiteOperatorSz (magnitude_spin_);
      onsite_operator_sp_  = CreateOnsiteOperatorSp (magnitude_spin_);
      onsite_operator_sm_  = CreateOnsiteOperatorSm (magnitude_spin_);
   }
   
};



}
}


#endif /* COMPNAL_MODEL_BASE_U1SPIN_1D_HPP_ */
