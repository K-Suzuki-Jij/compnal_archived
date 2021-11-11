//
//  xxz_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/06.
//

#ifndef COMPNAL_MODEL_XXZ_1D_HPP_
#define COMPNAL_MODEL_XXZ_1D_HPP_

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
class XXZ_1D {
   
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   using ValueType = RealType;
   
   XXZ_1D() {
      SetOnsiteOperator();
   };
   
   explicit XXZ_1D(const int system_size): XXZ_1D() {
      SetSystemSize(system_size);
   }
   
   XXZ_1D(const int system_size, const double magnitude_spin): XXZ_1D(system_size) {
      SetMagnitudeSpin(magnitude_spin);
   }
   
   XXZ_1D(const int system_size, const utility::BoundaryCondition bc): XXZ_1D(system_size) {
      SetBoundaryCondition(bc);
   }
   
   XXZ_1D(const int system_size, const double magnitude_spin, const utility::BoundaryCondition bc): XXZ_1D(system_size, magnitude_spin) {
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
         flag_recalc_basis_ = true;
         flag_recalc_ham_   = true;
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
         flag_recalc_basis_ = true;
         flag_recalc_ham_   = true;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetTotalSz(const double total_sz) {
      const int total_2sz = utility::DoubleTheNumber(total_sz);
      if (total_2sz < -system_size_*magnitude_2spin_ || system_size_*magnitude_2spin_ < total_2sz) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "There is no target space specified by total_2sz = " << total_2sz << std::endl;
         ss << "Please set as follows:" << std::endl;
         ss << -system_size_*magnitude_2spin_ << " <= total_2sz <= " << system_size_*magnitude_2spin_;
         throw std::runtime_error(ss.str());
      }
      if (total_2sz_ != total_2sz) {
         total_2sz_ = total_2sz;
         flag_recalc_basis_ = true;
         flag_recalc_ham_   = true;
         calculated_eigenvector_set_.clear();
      }
   }
   
   void SetBoundaryCondition(const utility::BoundaryCondition bc) {
      boundary_condition_ = bc;
      flag_recalc_ham_ = true;
   }
   
   void SetJz(const std::vector<RealType> &J_z) {
      if (J_z_ != J_z) {
         J_z_ = J_z;
         flag_recalc_ham_ = true;
      }
   }
   
   void SetJz(const RealType J_z) {
      if (J_z_.size() == 0) {
         J_z_.push_back(J_z);
         flag_recalc_ham_ = true;
      }
      else if (J_z_[0] != J_z) {
         J_z_[0] = J_z;
         flag_recalc_ham_ = true;
      }
   }
   
   void SetJxy(const std::vector<RealType> &J_xy) {
      if (J_xy_ != J_xy) {
         J_xy_ = J_xy;
         flag_recalc_ham_ = true;
      }
   }
   
   void SetJxy(const RealType J_xy) {
      if (J_xy_.size() == 0) {
         J_xy_.push_back(J_xy);
         flag_recalc_ham_ = true;
      }
      else if (J_xy_[0] != J_xy) {
         J_xy_[0] = J_xy;
         flag_recalc_ham_ = true;
      }
   }
   
   void SetHz(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(magnitude_2spin_, h_z_, D_z_);
         flag_recalc_ham_ = true;
      }
   }
   
   void SetDz(const RealType D_z) {
      if (D_z_ != D_z) {
         D_z_ = D_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(magnitude_2spin_, h_z_, D_z_);
         flag_recalc_ham_ = true;
      }
   }
   
   void SetFlagRecalcBasis(const bool flag) {
      flag_recalc_basis_ = flag;
   }
   
   void SetFlagRecalcHam(const bool flag) {
      flag_recalc_ham_ = flag;
   }
   
   void SetCalculatedEigenvectorSet(const std::size_t level) {
      calculated_eigenvector_set_.emplace(level);
   }
   
   void RemoveCalculatedEigenvectorSet(const std::size_t level) {
      calculated_eigenvector_set_.erase(level);
   }
   
   void PrintInfo() const {
      std::string bc = "None";
      if (boundary_condition_ == utility::BoundaryCondition::OBC) {
         bc = "OBC";
      }
      else if (boundary_condition_ == utility::BoundaryCondition::PBC) {
         bc = "PBC";
      }
      else if (boundary_condition_ == utility::BoundaryCondition::SSD) {
         bc = "SSD";
      }
      std::cout << "Print Heisenberg Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << bc                      << std::endl;
      std::cout << "system_size            = " << system_size_            << std::endl;
      std::cout << "magnitute_2spin        = " << magnitude_2spin_        << std::endl;
      std::cout << "total_2sz              = " << total_2sz_              << std::endl;
      std::cout << "dim_target             = " << CalculateTargetDim()    << std::endl;
      std::cout << "dim_onsite             = " << dim_onsite_             << std::endl;
      std::cout << "num_conserved_quantity = " << num_conserved_quantity_ << std::endl;
      
      std::cout << "Print Heisenberg Interaction" << std::endl;
      std::cout << "Sz-Sz Interaction: J_z =" << std::endl;
      for (std::size_t i = 0; i < J_z_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_z_.at(i) << std::endl;
      }
      std::cout << "Sx-Sx, Sy-Sy Interactions: J_xy =" << std::endl;
      for (std::size_t i = 0; i < J_xy_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_xy_.at(i) << std::endl;
      }
      std::cout << "External Magnetic Fields for the z-direction: h_z =" << h_z_ << std::endl;
      std::cout << "Uniaxial Anisotropy for the z-direction: D_z =" << D_z_ << std::endl;
   }
   
   void PrintBasisOnsite() const {
      const double magnitude_spin = magnitude_2spin_/2.0;
      for (int row = 0; row < dim_onsite_; ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin - row << ">" << std::endl;
      }
   }
   
   std::size_t CalculateTargetDim() const {
      const int max_total_2sz = system_size_*magnitude_2spin_;
      std::vector<std::vector<std::size_t>> dim(system_size_, std::vector<std::size_t>(max_total_2sz + 1));
      for (int s = -magnitude_2spin_; s <= magnitude_2spin_; s += 2) {
         dim[0][(s + magnitude_2spin_)/2] = 1;
      }
      for (int site = 1; site < system_size_; site++) {
         for (int s = -magnitude_2spin_; s <= magnitude_2spin_; s += 2) {
            for (int s_prev = -magnitude_2spin_*site; s_prev <= magnitude_2spin_*site; s_prev += 2) {
               const std::size_t a = dim[site    ][(s + s_prev + magnitude_2spin_*(site + 1))/2];
               const std::size_t b = dim[site - 1][(s_prev + magnitude_2spin_*site)/2];
               if (a >= UINT64_MAX - b) {
                  throw std::runtime_error("Overflow detected for sumation using uint64_t");
               }
               dim[site][(s + s_prev + magnitude_2spin_*(site + 1))/2] = a + b;
            }
         }
      }
      return dim[system_size_ - 1][(total_2sz_ + max_total_2sz)/2];
   }
   
   void GenerateBasis(std::vector<std::size_t> *basis, std::unordered_map<std::size_t, std::size_t> *basis_inv) const {
      if ((system_size_*magnitude_2spin_ - total_2sz_)%2 == 1) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      const int shifted_2sz = (system_size_*magnitude_2spin_ - total_2sz_)/2;
      const std::size_t dim_target = CalculateTargetDim();
      std::vector<std::vector<int>> partition_integers;
      utility::GenerateIntegerPartition(&partition_integers, shifted_2sz, magnitude_2spin_);
      
      std::vector<std::size_t> site_constant(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant[site] = static_cast<std::size_t>(std::pow(dim_onsite_, site));
      }
      
      std::vector<std::size_t>().swap(*basis);
      
#ifdef _OPENMP
      const int num_threads = omp_get_max_threads();
      std::vector<std::vector<std::size_t>> temp_basis(num_threads);
      
      for (std::size_t i = 0; i < partition_integers.size(); ++i) {
         const bool condition1 = (0 < partition_integers[i].size()) && (static_cast<int>(partition_integers[i].size()) <= system_size_);
         const bool condition2 = (partition_integers[i].size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            
            for (int j = static_cast<int>(partition_integers[i].size()); j < system_size_; ++j) {
               partition_integers[i].push_back(0);
            }
            
            const std::size_t size = utility::CalculateNumCombination(partition_integers[i]);
            std::vector<std::vector<int>> temp_partition_integer(num_threads);
            
#pragma omp parallel num_threads (num_threads)
            {
               const int thread_num = omp_get_thread_num();
               const std::size_t loop_begin = thread_num*size/num_threads;
               const std::size_t loop_end   = (thread_num + 1)*size/num_threads;
               temp_partition_integer[thread_num] = partition_integers[i];
               utility::CalculateNthPermutation(&temp_partition_integer[thread_num], loop_begin);
               
               for (std::size_t j = loop_begin; j < loop_end; ++j) {
                  std::size_t basis_global = 0;
                  for (std::size_t k = 0; k < temp_partition_integer[thread_num].size(); ++k) {
                     basis_global += temp_partition_integer[thread_num][k]*site_constant[k];
                  }
                  temp_basis[thread_num].push_back(basis_global);
                  std::next_permutation(temp_partition_integer[thread_num].begin(), temp_partition_integer[thread_num].end());
               }
            }
         }
      }
      
      for (std::size_t i = 0; i < temp_basis.size(); ++i) {
         basis->insert(basis->end(), temp_basis[i].begin(), temp_basis[i].end());
         std::vector<std::size_t>().swap(temp_basis[i]);
      }
      
#else
      basis->reserve(dim_target);
      
      for (std::size_t i = 0; i < partition_integers.size(); ++i) {
         const bool condition1 = (0 < partition_integers[i].size()) && (partition_integers[i].size() <= system_size_);
         const bool condition2 = (partition_integers[i].size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            
            for (std::size_t j = partition_integers[i].size(); j < system_size_; ++j) {
               partition_integers[i].push_back(0);
            }
            
            std::sort(partition_integers[i].begin(), partition_integers[i].end());
            
            do {
               std::size_t basis_global = 0;
               for (std::size_t j = 0; j < partition_integers[i].size(); ++j) {
                  basis_global += partition_integers[i][j]*site_constant[j];
               }
               basis->push_back(basis_global);
            } while (std::next_permutation(partition_integers[i].begin(), partition_integers[i].end()));
         }
      }
#endif
      
      if (basis->size() != dim_target) {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      std::sort(basis->begin(), basis->end());
      basis_inv->clear();
      
      for (std::size_t i = 0; i < basis->size(); ++i) {
         (*basis_inv)[(*basis)[i]] = i;
      }
   }
   
   inline utility::BoundaryCondition GetBoundaryCondition() const { return boundary_condition_;     }
   inline int GetSystemSize()           const { return system_size_;            }
   inline int GetDimOnsite()            const { return dim_onsite_;             }
   inline int GetMagnitude2Spin()       const { return magnitude_2spin_;        }
   inline int GetTotal2Sz()             const { return total_2sz_;              }
   inline int GetNumConservedQuantity() const { return num_conserved_quantity_; }
   inline double GetMagnitudeSpin()     const { return 0.5*magnitude_2spin_;    }
   inline double GetTotalSz()           const { return 0.5*total_2sz_;          }
   
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }
   inline const CRS &GetOnsiteOperatorSx () const { return onsite_operator_sx_ ; }
   inline const CRS &GetOnsiteOperatoriSy() const { return onsite_operator_isy_; }
   inline const CRS &GetOnsiteOperatorSz () const { return onsite_operator_sz_ ; }
   inline const CRS &GetOnsiteOperatorSp () const { return onsite_operator_sp_ ; }
   inline const CRS &GetOnsiteOperatorSm () const { return onsite_operator_sm_ ; }
   
   inline const std::vector<RealType> &GetJz () const { return J_z_ ; }
   inline const std::vector<RealType> &GetJxy() const { return J_xy_; }
   
   inline RealType GetJz (const std::size_t index) const { return J_z_ .at(index); }
   inline RealType GetJxy(const std::size_t index) const { return J_xy_.at(index); }
   
   inline RealType GetHz() const { return h_z_; }
   inline RealType GetDz() const { return D_z_; }
   
   inline bool GetFlagRecalcBasis() const { return flag_recalc_basis_; }
   inline bool GetFlagRecalcHam  () const { return flag_recalc_ham_; }
   inline const std::unordered_set<std::size_t> &GetCalculatedEigenvectorSet() const { return calculated_eigenvector_set_; }

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
   
   static CRS CreateOnsiteOperatorSzSz(const double magnitude_spin) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         const RealType val = magnitude_spin - row;
         if (val != 0.0) {
            matrix.val.push_back(val*val);
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
   
   static CRS CreateOnsiteOperatorHam(const double magnitude_spin, const RealType h_z = 0.0, const RealType D_z = 0.0) {
      const int magnitude_2spin = utility::DoubleTheNumber(magnitude_spin);
      const int dim_onsite      = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         const RealType val = h_z*(magnitude_spin - row) + D_z*D_z*(magnitude_spin - row);
         if (val != 0.0) {
            matrix.val.push_back(val);
            matrix.col.push_back(row);
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   
private:
   CRS onsite_operator_ham_;
   CRS onsite_operator_sx_;
   CRS onsite_operator_isy_;
   CRS onsite_operator_sz_;
   CRS onsite_operator_sp_;
   CRS onsite_operator_sm_;
   
   utility::BoundaryCondition boundary_condition_ = utility::BoundaryCondition::OBC;
   
   int system_size_     = 1;
   int total_2sz_       = 0;
   int dim_onsite_      = 2;
   int magnitude_2spin_ = 1;
   
   std::vector<RealType> J_z_  = {1.0};
   std::vector<RealType> J_xy_ = {1.0};
   RealType h_z_ = 0.0;
   RealType D_z_ = 0.0;
   
   const int num_conserved_quantity_ = 1;
   
   bool flag_recalc_basis_ = true;
   bool flag_recalc_ham_   = true;
   std::unordered_set<std::size_t> calculated_eigenvector_set_;
   
   void SetOnsiteOperator() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*magnitude_2spin_);
      onsite_operator_sx_  = CreateOnsiteOperatorSx (0.5*magnitude_2spin_);
      onsite_operator_isy_ = CreateOnsiteOperatoriSy(0.5*magnitude_2spin_);
      onsite_operator_sz_  = CreateOnsiteOperatorSz (0.5*magnitude_2spin_);
      onsite_operator_sp_  = CreateOnsiteOperatorSp (0.5*magnitude_2spin_);
      onsite_operator_sm_  = CreateOnsiteOperatorSm (0.5*magnitude_2spin_);
   }
   
};


}
}


#endif /* COMPNAL_MODEL_XXZ_1D_HPP_ */
