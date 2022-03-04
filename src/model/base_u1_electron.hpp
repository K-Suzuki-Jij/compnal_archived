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
//  Created by Kohei Suzuki on 2021/11/27.
//

#ifndef COMPNAL_MODEL_BASE_U1ELECTRON_HPP_
#define COMPNAL_MODEL_BASE_U1ELECTRON_HPP_

#include "../blas/all.hpp"
#include "../utility/all.hpp"
#include "../type/all.hpp"

#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

//! @brief The base class for electron systems with the U(1) symmetry.
//! Conserved quantity is the total sz:\n
//! \f[ \langle\hat{S}^{z}_{\rm tot}\rangle=
//! \sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f] and the total electrons \n
//! \f[ \langle\hat{N}^{\rm tot}_{\rm e}\rangle=
//! \sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f].
//! @tparam RealType Type of real values.
template<typename RealType>
class BaseU1Electron {
   
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
   //! @brief Constructor of BaseU1Electron class.
   BaseU1Electron() {
      SetOnsiteOperator();
   }
   
   //! @brief Constructor of BaseU1Electron class.
   //! @param total_electron The number of the total electrons
   //! \f[ \langle \hat{N}_{e}\rangle=
   //! \sum^{N}_{i=1}\langle\hat{n}_{i}\rangle=
   //! \sum^{N}_{i=1}\sum_{\sigma=\uparrow,\downarrow}\langle\hat{c}^{\dagger}_{i,\sigma}\hat{c}_{i,\sigma}\rangle\f],
   //! with \f$ N\f$ being the system size.
   explicit BaseU1Electron(const int total_electron): BaseU1Electron() {
      SetTotalElectron(total_electron);
   }
   
   //! @brief Constructor of BaseU1Electron class.
   //! @param total_electron The number of the total electrons
   //! \f[ \langle \hat{N}_{e}\rangle=
   //! \sum^{N}_{i=1}\langle\hat{n}_{i}\rangle=
   //! \sum^{N}_{i=1}\sum_{\sigma=\uparrow,\downarrow}\langle\hat{c}^{\dagger}_{i,\sigma}\hat{c}_{i,\sigma}\rangle\f],
   //! with \f$ N\f$ being the system size.
   //! @param total_sz The total sz
   //! \f[ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f],
   //! with \f$ N\f$ being the system size.
   explicit BaseU1Electron(const int total_electron, const HalfInt total_sz): BaseU1Electron(total_electron) {
      SetTotalSz(total_sz);
   }
   
   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------
   //! @brief Set the target Hilbert space specified by the total sz.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle=\sum^{N}_{i=1}\langle\hat{S}^{z}_{i}\rangle \f$.
   void SetTotalSz(const HalfInt total_sz) {
      total_sz_ = total_sz;
   }
   
   //! @brief Set the number of total electrons.
   //! @param total_electron The number of total electrons is represented by the expectation value of the following operator:
   //! \f[ \hat{N}_{\rm e}=\sum^{N}_{i=1}\hat{n}_{i} \f]
   void SetTotalElectron(const int total_electron) {
      if (total_electron < 0) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "total_electron must be a non-negative integer" << std::endl;
         throw std::runtime_error(ss.str());
      }
      total_electron_ = total_electron;
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
         ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
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
      if (row == col && 0 <= row && row < 4 && 0 <= col && col < 4) {
         return {+0 + total_electron_, +0.0 + total_sz_};
      }
      else if (row == 0 && col == 1) {
         return {-1 + total_electron_, -0.5 + total_sz_};
      }
      else if (row == 0 && col == 2) {
         return {-1 + total_electron_, +0.5 + total_sz_};
      }
      else if (row == 0 && col == 3) {
         return {-2 + total_electron_, +0.0 + total_sz_};
      }
      else if (row == 1 && col == 0) {
         return {+1 + total_electron_, +0.5 + total_sz_};
      }
      else if (row == 1 && col == 2) {
         return {+0 + total_electron_, +1.0 + total_sz_};
      }
      else if (row == 1 && col == 3) {
         return {-1 + total_electron_, +0.5 + total_sz_};
      }
      else if (row == 2 && col == 0) {
         return {+1 + total_electron_, -0.5 + total_sz_};
      }
      else if (row == 2 && col == 1) {
         return {+0 + total_electron_, -1.0 + total_sz_};
      }
      else if (row == 2 && col == 3) {
         return {-1 + total_electron_, -0.5 + total_sz_};
      }
      else if (row == 3 && col == 0) {
         return {+2 + total_electron_, +0.0 + total_sz_};
      }
      else if (row == 3 && col == 1) {
         return {+1 + total_electron_, -0.5 + total_sz_};
      }
      else if (row == 3 && col == 2) {
         return {+1 + total_electron_, +0.5 + total_sz_};
      }
      else {
         std::stringstream ss;
         ss << "Error in " << __func__ << std::endl;
         ss << "The dimenstion of the matrix must be 4";
         throw std::runtime_error(ss.str());
      }
   }
   
   //! @brief Generate bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @param quantum_number The pair of the total electron \f$ \langle\hat{N}_{\rm e}\rangle \f$ and total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$
   //! @param flag_display_info If true, display the progress status. Set ture by default.
   //! @return Corresponding basis.
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
      const int total_2sz      = 2*total_sz;
      
      std::vector<std::int64_t> site_constant(system_size);
      for (int site = 0; site < system_size; ++site) {
         site_constant[site] = static_cast<std::int64_t>(std::pow(dim_onsite_, site));
      }
      
      const int max_n_up_down = static_cast<int>(total_electron/2);
      
      std::vector<std::vector<int>> partition_integers;
      for (int n_up_down = 0; n_up_down <= max_n_up_down; ++n_up_down) {
         const int n_up   = (total_electron - 2*n_up_down + total_2sz)/2;
         const int n_down = (total_electron - 2*n_up_down - total_2sz)/2;
         const int n_vac  = system_size - total_electron + n_up_down;
         if (0 <= n_up && 0 <= n_down && 0 <= n_vac) {
            std::vector<int> integer_list(system_size);
            for (int s = 0; s < n_vac; ++s) {
               integer_list[s] = 0;
            }
            for (int s = 0; s < n_up; ++s) {
               integer_list[s + n_vac] = 1;
            }
            for (int s = 0; s < n_down; ++s) {
               integer_list[s + n_vac + n_up] = 2;
            }
            for (int s = 0; s < n_up_down; ++s) {
               integer_list[s + n_vac + n_up + n_down] = 3;
            }
            partition_integers.push_back(integer_list);
         }
      }
      
      const std::int64_t dim_target = CalculateTargetDim(system_size, quantum_number.first, quantum_number.second);
      std::vector<std::int64_t> basis;
      basis.reserve(dim_target);
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<std::vector<std::int64_t>> temp_basis(num_threads);
      for (const auto &integer_list: partition_integers) {
         const std::int64_t size = utility::CalculateNumPermutation(integer_list);
#pragma omp parallel num_threads (num_threads)
         {
            const int thread_num = omp_get_thread_num();
            const std::int64_t loop_begin = thread_num*size/num_threads;
            const std::int64_t loop_end   = (thread_num + 1)*size/num_threads;
            std::vector<int> n_th_integer_list = utility::GenerateNthPermutation(integer_list, loop_begin);
            for (std::int64_t j = loop_begin; j < loop_end; ++j) {
               std::int64_t basis_global = 0;
               for (std::size_t k = 0; k < n_th_integer_list.size(); ++k) {
                  basis_global += n_th_integer_list[k]*site_constant[k];
               }
               temp_basis[thread_num].push_back(basis_global);
               std::next_permutation(n_th_integer_list.begin(), n_th_integer_list.end());
            }
         }
      }
      for (auto &&it: temp_basis) {
         basis.insert(basis.end(), it.begin(), it.end());
         std::vector<std::int64_t>().swap(it);
      }
#else
      for (std::size_t i = 0; i < partition_integers.size(); ++i) {
         auto &integer_list = partition_integers[i];
         std::sort(integer_list.begin(), integer_list.end());
         do {
            std::int64_t basis_global = 0;
            for (std::size_t j = 0; j < integer_list.size(); ++j) {
               basis_global += integer_list[j]*site_constant[j];
            }
            basis.push_back(basis_global);
         } while (std::next_permutation(integer_list.begin(), integer_list.end()));
      }
#endif
      basis.shrink_to_fit();
      
      if (static_cast<std::int64_t>(basis.size()) != dim_target) {
         std::stringstream ss;
         ss << "Unknown error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
         ss << "basis.size()=" << basis.size() << ", but dim_target=" << dim_target << std::endl;
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
      return ValidateQNumber(system_size, quantum_number.first, quantum_number.second);
   }
   
   //! @brief Calculate the dimension of the target Hilbert space specified by
   //! the system size \f$ N\f$, the number of the total electrons \f$ \langle\hat{N}_{\rm e}\rangle\f$, and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! @return The dimension of the target Hilbert space.
   std::int64_t CalculateTargetDim(const int system_size, const QType &quantum_number) const {
      return CalculateTargetDim(system_size, quantum_number.first, quantum_number.second);
   }
   
   
   //! @brief Check if there is a subspace specified by the input quantum numbers.
   //! @param system_size The system size \f$ N\f$.
   //! @param total_electron The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   //! @param total_sz The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return ture if there exists corresponding subspace, otherwise false.
   static bool ValidateQNumber(const int system_size, const int total_electron, const HalfInt total_sz) {
      if (system_size <= 0 || total_electron < 0) {
         return false;
      }
      const int total_2sz = 2*total_sz;
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
   static std::int64_t CalculateTargetDim(const int system_size, const int total_electron, const HalfInt total_sz) {
      if (!ValidateQNumber(system_size, total_electron, total_sz)) {
         return 0;
      }
      if (system_size <= 0) {
         return 0;
      }
      const int total_2sz = 2*total_sz;
      const std::vector<std::vector<std::int64_t>> binom = utility::GenerateBinomialTable(system_size);
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
      
      const RealType val = RealType{1.0};
      const int dim_onsite = 4;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; row++) {
         for (int col = 0; col < dim_onsite; col++) {
            if ((col == 1 && row == 0) || (col == 3 && row == 2)) {
               matrix.col.push_back(col);
               matrix.val.push_back(val);
            }
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      
      matrix.tag = type::CRSTag::FERMION;
      
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
      
      const RealType val = RealType{1.0};
      const int dim_onsite = 4;
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; row++) {
         for (int col = 0; col < dim_onsite; col++) {
            if (col == 2 && row == 0) {
               matrix.col.push_back(col);
               matrix.val.push_back(val);
            }
            else if (col == 3 && row == 1) {
               matrix.col.push_back(col);
               matrix.val.push_back(-val);
            }
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      
      matrix.tag = type::CRSTag::FERMION;
      
      return matrix;
   }
   
   //! @brief Generate the creation operator for the electrons with the up spin
   //! \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorCUpDagger() {
      return type::CalculateTransposedMatrix(CreateOnsiteOperatorCUp());
   }
   
   //! @brief Generate the creation operator for the electrons with the down spin
   //! \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorCDownDagger() {
      return type::CalculateTransposedMatrix(CreateOnsiteOperatorCDown());
   }
   
   //! @brief Generate the number operator for the electrons with the up spin
   //! \f$ \hat{n}_{\uparrow}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\uparrow}\f$.
   static CRS CreateOnsiteOperatorNCUp() {
      return CreateOnsiteOperatorCUpDagger()*CreateOnsiteOperatorCUp();
   }
   
   //! @brief Generate the number operator for the electrons with the down spin
   //! \f$ \hat{n}_{\downarrow}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}_{\downarrow}\f$.
   static CRS CreateOnsiteOperatorNCDown() {
      return CreateOnsiteOperatorCDownDagger()*CreateOnsiteOperatorCDown();
   }
   
   //! @brief Generate the number operator for the electrons
   //! \f$ \hat{n}=\hat{n}_{\uparrow} + \hat{n}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{n}\f$.
   static CRS CreateOnsiteOperatorNC() {
      return CreateOnsiteOperatorNCUp() + CreateOnsiteOperatorNCDown();
   }
   
   //! @brief Generate the spin operator for the x-direction for the electrons
   //! \f$ \hat{s}^{x}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} + \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{x}\f$.
   static CRS CreateOnsiteOperatorSx() {
      return RealType{0.5}*(CreateOnsiteOperatorSp() + CreateOnsiteOperatorSm());
   }
   
   //! @brief Generate the spin operator for the y-direction for the electrons
   //! \f$ i\hat{s}^{y}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @return The matrix of \f$ i\hat{s}^{y}\f$.
   static CRS CreateOnsiteOperatoriSy() {
      return RealType{0.5}*(CreateOnsiteOperatorSp() - CreateOnsiteOperatorSm());
   }
   
   //! @brief Generate the spin operator for the z-direction for the electrons
   //! \f$ \hat{s}^{z}=\frac{1}{2}(\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\uparrow} - \hat{c}^{\dagger}_{\downarrow}\hat{c}_{\downarrow})\f$.
   //! @return The matrix of \f$ \hat{s}^{z}\f$.
   static CRS CreateOnsiteOperatorSz() {
      return RealType{0.5}*(CreateOnsiteOperatorNCUp() - CreateOnsiteOperatorNCDown());
   }
   
   //! @brief Generate the raising operator for spin of the electrons
   //! \f$ \hat{s}^{+}=\hat{c}^{\dagger}_{\uparrow}\hat{c}_{\downarrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{+}\f$.
   static CRS CreateOnsiteOperatorSp() {
      return CreateOnsiteOperatorCUpDagger()*CreateOnsiteOperatorCDown();
   }
   
   //! @brief Generate the lowering operator for spin of the electrons
   //! \f$ \hat{s}^{-}=\hat{c}^{\dagger}_{\downarrow}\hat{c}_{\uparrow}\f$.
   //! @return The matrix of \f$ \hat{s}^{-}\f$.
   static CRS CreateOnsiteOperatorSm() {
      return CreateOnsiteOperatorCDownDagger()*CreateOnsiteOperatorCUp();
   }
      
   //! @brief Get dimension of the local Hilbert space, 4.
   //! @return The dimension of the local Hilbert space, 4.
   inline int GetDimOnsite() const { return dim_onsite_; }
   
   //! @brief Get the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   //! @return The total sz.
   inline HalfInt GetTotalSz() const { return total_sz_; }
   
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

   
protected:
   //! @brief The total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle\f$.
   HalfInt total_sz_ = 0;
   
   //! @brief The total electron \f$ \langle\hat{N}_{\rm e}\rangle\f$.
   int total_electron_ = 0;
   
   //! @brief The dimension of the local Hilbert space, 4.
   const int dim_onsite_ = 4;
   
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
   
   //------------------------------------------------------------------
   //----------------------Private Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Set onsite operators.
   void SetOnsiteOperator() {
      onsite_operator_c_up_          = CreateOnsiteOperatorCUp();
      onsite_operator_c_down_        = CreateOnsiteOperatorCDown();
      onsite_operator_c_up_dagger_   = CreateOnsiteOperatorCUpDagger();
      onsite_operator_c_down_dagger_ = CreateOnsiteOperatorCDownDagger();
      onsite_operator_nc_up_         = CreateOnsiteOperatorNCUp();
      onsite_operator_nc_down_       = CreateOnsiteOperatorNCDown();
      onsite_operator_nc_            = CreateOnsiteOperatorNC();
      onsite_operator_sx_            = CreateOnsiteOperatorSx();
      onsite_operator_isy_           = CreateOnsiteOperatoriSy();
      onsite_operator_sz_            = CreateOnsiteOperatorSz();
      onsite_operator_sp_            = CreateOnsiteOperatorSp();
      onsite_operator_sm_            = CreateOnsiteOperatorSm();
   }
   
};

} //namespace model
} //namespace compnal


#endif /* COMPNAL_MODEL_BASE_U1ELECTRON_HPP_ */
