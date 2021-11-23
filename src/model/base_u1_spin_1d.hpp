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

template<typename RealType>
class BaseU1Spin_1D {
   
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   using ValueType = RealType;
   
   BaseU1Spin_1D() {
      SetOnsiteOperator();
   }
   
   explicit BaseU1Spin_1D(const int system_size): BaseU1Spin_1D() {
      SetSystemSize(system_size);
   }
   
   BaseU1Spin_1D(const int system_size,
                 const double magnitude_spin
                 ): BaseU1Spin_1D(system_size) {
      SetMagnitudeSpin(magnitude_spin);
   }
   
   BaseU1Spin_1D(const int system_size,
                 const utility::BoundaryCondition bc
                 ): BaseU1Spin_1D(system_size) {
      SetBoundaryCondition(bc);
   }
   
   BaseU1Spin_1D(const int system_size,
                 const double magnitude_spin,
                 const utility::BoundaryCondition bc
                 ): BaseU1Spin_1D(system_size, magnitude_spin) {
      SetBoundaryCondition(bc);
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
   
   void SetTotalSz(const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      if (isValidQNumber(total_2sz)) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "There is no target space specified by total_sz = " << total_sz << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetBoundaryCondition(const utility::BoundaryCondition bc) {
      boundary_condition_ = bc;
   }
   
   void SetCalculatedEigenvectorSet(const std::int64_t level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   void RemoveCalculatedEigenvectorSet(const std::int64_t level) {
      calculated_eigenvector_set_.erase(level);
   }
   
   bool isValidQNumber(const int total_2sz) const {
      const bool c1 = ((system_size_*magnitude_2spin_ - total_2sz)%2 == 0);
      const bool c2 = (-system_size_*magnitude_2spin_ < total_2sz);
      const bool c3 = (total_2sz < +system_size_*magnitude_2spin_);
      if (c1 && c2 && c3) {
         return true;
      }
      else {
         return false;
      }
   }
   
   void PrintBasisOnsite() const {
      const double magnitude_spin = magnitude_2spin_/2.0;
      for (int row = 0; row < dim_onsite_; ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin - row << ">" << std::endl;
      }
   }
   
   std::int64_t CalculateTargetDim() const {
      return CalculateTargetDim(0.5*total_2sz_);
   }
   
   std::int64_t CalculateTargetDim(const double total_sz) const {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      if (isValidQNumber(total_2sz) == false) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const int max_total_2sz = system_size_*magnitude_2spin_;
      std::vector<std::vector<std::int64_t>> dim(system_size_, std::vector<std::int64_t>(max_total_2sz + 1));
      for (int s = -magnitude_2spin_; s <= magnitude_2spin_; s += 2) {
         dim[0][(s + magnitude_2spin_)/2] = 1;
      }
      for (int site = 1; site < system_size_; site++) {
         for (int s = -magnitude_2spin_; s <= magnitude_2spin_; s += 2) {
            for (int s_prev = -magnitude_2spin_*site; s_prev <= magnitude_2spin_*site; s_prev += 2) {
               const std::int64_t a = dim[site    ][(s + s_prev + magnitude_2spin_*(site + 1))/2];
               const std::int64_t b = dim[site - 1][(s_prev + magnitude_2spin_*site)/2];
               if (a >= INT64_MAX - b) {
                  throw std::runtime_error("Overflow detected for sumation using uint64_t");
               }
               dim[site][(s + s_prev + magnitude_2spin_*(site + 1))/2] = a + b;
            }
         }
      }
      return dim[system_size_ - 1][(total_2sz + max_total_2sz)/2];
   }
   
   void GenerateBasis() {
      GenerateBasis(0.5*total_2sz_);
   }
   
   void GenerateBasis(const double total_sz) {
      const auto start = std::chrono::system_clock::now();
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      
      if (isValidQNumber(total_2sz) == false) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (bases_.count(total_2sz) != 0) {
         return;
      }
      
      std::cout << "Generating Basis..." << std::flush;
      
      const int shifted_2sz = (system_size_*magnitude_2spin_ - total_2sz)/2;
      const std::int64_t dim_target = CalculateTargetDim(total_2sz);
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
   
   std::unordered_set<int> GenerateTargetSector(const CRS &m_1, const CRS &m_2) const {
      std::unordered_set<int> level_set;
      std::unordered_set<int> level_set_m1;
      std::unordered_set<int> level_set_m2;
      for (std::int64_t i = 0; i < m_1.row_dim; ++i) {
         for (std::int64_t j = m_1.row[i]; j < m_1.row[i + 1]; ++j) {
            if (m_1.val[j] != 0.0) {
               level_set.emplace(static_cast<int>((m_1.col[j] - i) + total_2sz_));
               level_set_m1.emplace(static_cast<int>((m_1.col[j] - i) + total_2sz_));
            }
         }
      }
      
      for (std::int64_t i = 0; i < m_2.row_dim; ++i) {
         for (std::int64_t j = m_2.row[i]; j < m_2.row[i + 1]; ++j) {
            if (m_2.val[j] != 0.0) {
               level_set.emplace(static_cast<int>((m_2.col[j] - i) + total_2sz_));
               level_set_m2.emplace(static_cast<int>((m_2.col[j] - i) + total_2sz_));
            }
         }
      }
      
      std::unordered_set<int> level_set_intersection;
      
      for (const auto &level: level_set) {
         if (level_set_m1.count(level) != 0 && level_set_m2.count(level)) {
            level_set_intersection.emplace(level);
         }
      }
      
      return level_set_intersection;
   }
   
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
   
   inline utility::BoundaryCondition GetBoundaryCondition() const { return boundary_condition_; }
   inline int GetSystemSize()           const { return system_size_;            }
   inline int GetDimOnsite()            const { return dim_onsite_;             }
   inline int GetMagnitude2Spin()       const { return magnitude_2spin_;        }
   inline int GetTotal2Sz()             const { return total_2sz_;              }
   inline double GetTotalSz()           const { return 0.5*total_2sz_;          }
   inline int GetNumConservedQuantity() const { return num_conserved_quantity_; }
   inline double GetMagnitudeSpin()     const { return 0.5*magnitude_2spin_;    }
   
   inline const CRS &GetOnsiteOperatorSx () const { return onsite_operator_sx_ ; }
   inline const CRS &GetOnsiteOperatoriSy() const { return onsite_operator_isy_; }
   inline const CRS &GetOnsiteOperatorSz () const { return onsite_operator_sz_ ; }
   inline const CRS &GetOnsiteOperatorSp () const { return onsite_operator_sp_ ; }
   inline const CRS &GetOnsiteOperatorSm () const { return onsite_operator_sm_ ; }
   
   inline const std::unordered_set<int> &GetCalculatedEigenvectorSet() const {
      return calculated_eigenvector_set_;
   }
   inline const std::vector<std::int64_t> &GetBasis(const double total_sz) const {
      return bases_.at(utility::DoubleTheNumber(total_sz));
   }
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetBasisInv(const double total_sz) const {
      return bases_inv_.at(utility::DoubleTheNumber(total_sz));
   }
   inline const std::vector<std::int64_t> &GetTargetBasis() const {
      return bases_.at(total_2sz_);
   }
   inline const std::unordered_map<std::int64_t, std::int64_t> &GetTargetBasisInv() const {
      return bases_inv_.at(total_2sz_);
   }
   
protected:
   CRS onsite_operator_sx_;
   CRS onsite_operator_isy_;
   CRS onsite_operator_sz_;
   CRS onsite_operator_sp_;
   CRS onsite_operator_sm_;
   
   utility::BoundaryCondition boundary_condition_ = utility::BoundaryCondition::OBC;
   
   int system_size_     = 0;
   int total_2sz_       = 0;
   int dim_onsite_      = 2;
   int magnitude_2spin_ = 1;
   
   const int num_conserved_quantity_ = 1;
   
   std::unordered_set<int> calculated_eigenvector_set_;
   
   std::unordered_map<int, std::vector<std::int64_t>> bases_;
   std::unordered_map<int, std::unordered_map<std::int64_t, std::int64_t>> bases_inv_;
   
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
