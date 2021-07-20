//
//  heisenberg_model_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef heisenberg_model_1d_hpp
#define heisenberg_model_1d_hpp

#include "model_utility.hpp"
#include "sparse_matrix.hpp"

#include <cmath>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace compnal {
namespace model {

template<typename RealType>
class Heisenberg1D {

   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   using ValueType = RealType;
   
   explicit Heisenberg1D(const int system_size) {
      SetSystemSize(system_size);
      SetOperator();
   }
   
   Heisenberg1D(const int system_size, const double magnitude_spin) {
      SetSystemSize(system_size);
      magnitude_2spin_ = DoubleTheNumber(magnitude_spin);
      dim_onsite_ = magnitude_2spin_ + 1;
      SetOperator();
   }
   
   Heisenberg1D(const int system_size, const BoundaryCondition bc) {
      SetSystemSize(system_size);
      SetBoundaryCondition(bc);
      SetOperator();
   }
   
   Heisenberg1D(const int system_size, const double magnitude_spin, const BoundaryCondition bc) {
      SetSystemSize(system_size);
      SetBoundaryCondition(bc);
      magnitude_2spin_ = DoubleTheNumber(magnitude_spin);
      dim_onsite_      = magnitude_2spin_ + 1;
      SetOperator();
   }
      
   inline BoundaryCondition GetBoundaryCondition()    const { return boundary_condition_;     }
   inline int               GetSystemSize()           const { return system_size_;            }
   inline int               GetDimOnsite()            const { return dim_onsite_;             }
   inline double            GetMagnitudeSpin()        const { return magnitude_2spin_/2.0;    }
   inline double            GetTotalSz()              const { return total_2sz_/2.0;          }
   inline int               GetNumConservedQuantity() const { return num_conserved_quantity_; }
   inline bool              GetFlagRecalcBasis()      const { return flag_recalc_basis;       }
   
   const CRS &GetOperatorHam() const { return ham_; }
   const CRS &GetOperatorSx () const { return sx_ ; }
   const CRS &GetOperatoriSy() const { return isy_; }
   const CRS &GetOperatorSz () const { return sz_ ; }
   const CRS &GetOperatorSp () const { return sp_ ; }
   const CRS &GetOperatorSm () const { return sm_ ; }
   
   const std::vector<RealType> &GetJz()  const { return J_z_ ; }
   const std::vector<RealType> &GetJxy() const { return J_xy_; }
   RealType GetJz (const std::size_t index) const { return J_z_[index];  }
   RealType GetJxy(const std::size_t index) const { return J_xy_[index]; }
   RealType GetHz () const { return h_z_; }
   RealType GetDz () const { return D_z_; }
   
   void SetMagnitudeSpin(const double magnitude_spin) {
      const int magnitude_2spin = DoubleTheNumber(magnitude_spin);
      if (magnitude_2spin <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__      << std::endl;
         ss << "Please set magnitude_2spin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_2spin_ != magnitude_2spin) {
         magnitude_2spin_ = magnitude_2spin;
         dim_onsite_      = magnitude_2spin + 1;
         SetOperator();
         flag_recalc_basis = true;
      }
   }
   
   void SetTotalSz(const double total_sz) {
      const int total_2sz = DoubleTheNumber(total_sz);
      if (system_size_ <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "Please set system_size > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (magnitude_2spin_ <= 0) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__      << std::endl;
         ss << "Please set magnitude_2spin > 0" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (total_2sz < -system_size_*magnitude_2spin_ || system_size_*magnitude_2spin_ < total_2sz) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__  << std::endl;
         ss << "There is no target space specified by total_2sz = " << total_2sz << std::endl;
         ss << "Please set as follows:" << std::endl;
         ss << -system_size_*magnitude_2spin_ << " <= total_2sz <= " << system_size_*magnitude_2spin_;
         throw std::runtime_error(ss.str());
      }
      total_2sz_ = total_2sz;
      flag_recalc_basis = true;
   }
   
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         std::stringstream ss;
         ss << "system_size must be more than 0" << std::endl;
         ss << "system_size=" << system_size << "is not allowed" << std::endl;
         throw std::runtime_error(ss.str());
      }
      system_size_ = system_size;
      flag_recalc_basis = true;
   }
   
   void SetBoundaryCondition(const BoundaryCondition bc) {
      boundary_condition_ = bc;
   }
   
   void SetJz(const std::vector<RealType> &J_z) {
      J_z_  = J_z;
   }
   void SetJz(const RealType J_z) {
      J_z_  = std::vector<RealType>{J_z};
   }
   
   template<typename... Args>
   void SetJz(Args... args) {
      J_z_ = std::vector<RealType>{args...};
   }
   
   void SetJxy(const std::vector<RealType> &J_xy) {
      J_xy_  = J_xy;
   }
   void SetJxy(const RealType J_xy) {
      J_xy_  = std::vector<RealType>{J_xy};
   }
   
   template<typename... Args>
   void SetJxy(Args... args) {
      J_xy_ = std::vector<RealType>{args...};
   }
   
   void SetHz(const RealType h_z) {
      h_z_ = h_z;
      ham_ = CreateOperatorHam(static_cast<double>(magnitude_2spin_/2.0));
   }
   
   void SetDz(const RealType D_z) {
      D_z_ = D_z;
      ham_ = CreateOperatorHam(static_cast<double>(magnitude_2spin_/2.0));
   }
   
   void PrintInfo() const {
      std::string bc = "None";
      if (boundary_condition_ == BoundaryCondition::OBC) {
         bc = "OBC";
      }
      else if (boundary_condition_ == BoundaryCondition::PBC) {
         bc = "PBC";
      }
      else if (boundary_condition_ == BoundaryCondition::SSD) {
         bc = "SSD";
      }
      std::cout << "Print Heisenberg Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << bc                      << std::endl;
      std::cout << "system_size            = " << system_size_            << std::endl;
      std::cout << "magnitute_spin         = " << magnitude_2spin_/2.0    << std::endl;
      std::cout << "total_ sz              = " << total_2sz_/2.0          << std::endl;
      std::cout << "dim_target             = " << GetDim()                << std::endl;
      std::cout << "dim_onsite             = " << dim_onsite_             << std::endl;
      std::cout << "num_conserved_quantity = " << num_conserved_quantity_ << std::endl;
      
      std::cout << "Print Heisenberg Interaction" << std::endl;
      std::cout << "Sz-Sz Interaction: J_z =" << std::endl;
      for (int64_t i = 0; i < J_z_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_z_.at(i) << std::endl;
      }
      std::cout << "Sx-Sx, Sy-Sy Interactions: J_xy =" << std::endl;
      for (int64_t i = 0; i < J_xy_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_xy_.at(i) << std::endl;
      }
      std::cout << "External Magnetic Fields for the z-direction: h_z =" << h_z_ << std::endl;
      std::cout << "Uniaxial Anisotropy for the z-direction: D_z =" << D_z_ << std::endl;
   }
   
   void PrintBasisOnsite() const {
      const double magnitude_spin = magnitude_2spin_/2.0;
      for (int64_t row = 0; row < sz_.GetRowDim(); ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin - row << ">" << std::endl;
      }
   }
   
   void PrintBasisOnsite(int basis, const bool flag_new_line = true) const {
      if (basis < 0 || dim_onsite_ < basis) {
         return;
      }
      const double magnitude_spin = magnitude_2spin_/2.0;
      std::cout << "|Sz=" << std::showpos << magnitude_spin - basis << std::noshowpos << ">";
      if (flag_new_line) {
         std::cout << std::endl;
      }
   }
      
   int64_t GetDim() const {
      const int system_size     = GetSystemSize();
      const int total_2sz       = total_2sz_;
      const int magnitude_2spin = magnitude_2spin_;
      const int max_total_2sz   = system_size*magnitude_2spin;
      
      std::vector<std::vector<int64_t>> dim;
      dim.assign(system_size, std::vector<int64_t>(max_total_2sz + 1));
      
      for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
         dim[0][(s + magnitude_2spin)/2] = 1;
      }
      
      for (int site = 1; site < system_size; site++) {
         for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
            for (int s_prev = -magnitude_2spin*site; s_prev <= magnitude_2spin*site; s_prev += 2) {
               int64_t a = dim[site][(s + s_prev + magnitude_2spin*(site + 1))/2];
               int64_t b = dim[site - 1][(s_prev + magnitude_2spin*site)/2];
               CheckSumOverflowINT64(a, b);
               dim[site][(s + s_prev + magnitude_2spin*(site + 1))/2] = a + b;
            }
         }
      }
      return dim[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   void GenerateBasis(std::vector<int64_t> *basis, std::unordered_map<int64_t, int64_t> *basis_inv) {
      if ((system_size_*magnitude_2spin_ - total_2sz_)%2 == 1) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid parameters (system_size or magnitude_spin or total_sz)" << std::endl;
         throw std::runtime_error(ss.str());
      }
      const int shifted_2sz = (system_size_*magnitude_2spin_ - total_2sz_)/2;
      const int64_t dim_target = GetDim();
      std::vector<std::vector<int>> partition_integers;
      model::GenerateIntegerPartition(&partition_integers, shifted_2sz, magnitude_2spin_);
      
      std::vector<int64_t> site_constant(system_size_);
      for (int site = 0; site < system_size_; ++site) {
         site_constant[site] = CalculatePower(dim_onsite_, site);
      }
      
      std::vector<int64_t>().swap(*basis);
      
#ifdef _OPENMP
      const int num_threads = omp_get_num_threads();
      std::vector<std::vector<int64_t>> bases(num_threads);
      
      for (int64_t i = 0; i < partition_integers.size(); ++i) {
         const bool condition1 = (0 < partition_integers[i].size()) && (partition_integers[i].size() <= system_size_);
         const bool condition2 = (partition_integers[i].size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            
            for (int64_t j = partition_integers[i].size(); j < system_size_; ++j) {
               partition_integers[i].push_back(0);
            }
            
            const int64_t size = CalculateNumCombination(partition_integers[i]);
            
#pragma omp parallel num_threads (num_threads)
            {
               const int thread_num = omp_get_thread_num();
               const int64_t loop_begin = thread_num*size/num_threads;
               const int64_t loop_end   = (thread_num + 1)*size/num_threads;
               
               std::vector<int> temp_partition_integer = partition_integers[i];
               NthPermutation(temp_partition_integer, loop_begin);
               
               for (int64_t j = loop_begin; j < loop_end; ++j) {
                  int64_t basis_onsite = 0;
                  for (int64_t k = 0; k < partition_integers[i].size(); ++k) {
                     basis_onsite += temp_partition_integer[k]*site_constant[k];
                  }
                  bases[thread_num].push_back(basis_onsite);
                  std::next_permutation(temp_partition_integer.begin(), temp_partition_integer.end());
               }
            }
         }
      }
            
      for (int64_t i = 0; i < bases.size(); ++i) {
         basis->insert(basis->end(), bases[i].begin(), bases[i].end());
         std::vector<int64_t>().swap(basis[i]);
      }
      
#else
      basis->reserve(dim_target);
         
      for (int64_t i = 0; i < partition_integers.size(); ++i) {
         const bool condition1 = (0 < partition_integers[i].size()) && (partition_integers[i].size() <= system_size_);
         const bool condition2 = (partition_integers[i].size() == 0) && (shifted_2sz  == 0);
         if (condition1 || condition2) {
            
            for (int64_t j = partition_integers[i].size(); j < system_size_; ++j) {
               partition_integers[i].push_back(0);
            }
            
            std::sort(partition_integers[i].begin(), partition_integers[i].end());
            
            do {
               int64_t basis_onsite = 0;
               for (int64_t j = 0; j < partition_integers[i].size(); ++j) {
                  basis_onsite += partition_integers[i][j]*site_constant[j];
               }
               basis->push_back(basis_onsite);
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
      
      for (int64_t i = 0; i < basis->size(); ++i) {
         (*basis_inv)[(*basis)[i]] = i;
      }

      flag_recalc_basis = false;
   }
   
   static CRS CreateOperatorSx(double magnitude_spin) {
      const int dim_onsite = DoubleTheNumber(magnitude_spin) + 1;
      CRS matrix(dim_onsite, dim_onsite);
      RealType a = 0;
      RealType b = 1;
      
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
   
   static CRS CreateOperatoriSy(double magnitude_spin) {
      const int dim_onsite = DoubleTheNumber(magnitude_spin) + 1;
      CRS matrix(dim_onsite, dim_onsite);
      RealType a = 0;
      RealType b = 1;
      
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
   
   static CRS CreateOperatorSz(double magnitude_spin) {
      const int dim_onsite = DoubleTheNumber(magnitude_spin) + 1;
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         RealType val = magnitude_spin - row;
         if (val != 0.0) {
            matrix.val.push_back(val);
            matrix.col.push_back(row);
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   static CRS CreateOperatorSzSz(double magnitude_spin) {
      const int dim_onsite = DoubleTheNumber(magnitude_spin) + 1;
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         RealType val = magnitude_spin - row;
         if (val != 0.0) {
            matrix.val.push_back(val*val);
            matrix.col.push_back(row);
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }
   
   static CRS CreateOperatorSp(double magnitude_spin) {
      return CreateOperatorSx(magnitude_spin) + CreateOperatoriSy(magnitude_spin);
   }
   
   static CRS CreateOperatorSm(double magnitude_spin) {
      return CreateOperatorSx(magnitude_spin) - CreateOperatoriSy(magnitude_spin);
   }
   
private:
   CRS ham_, sx_, isy_, sz_, sp_, sm_;
   
   BoundaryCondition boundary_condition_ = BoundaryCondition::OBC;
   
   int system_size_     = 0;
   int total_2sz_       = 0;
   int dim_onsite_      = 2;
   int magnitude_2spin_ = 1;
   
   const int num_conserved_quantity_ = 1;
   
   std::vector<RealType> J_z_  = {1.0};
   std::vector<RealType> J_xy_ = {1.0};
   RealType h_z_ = 0.0;
   RealType D_z_ = 0.0;
   
   bool flag_recalc_basis = true;
   
   void SetOperator() {
      double magnitude_spin = static_cast<double>(magnitude_2spin_/2.0);
      ham_ = CreateOperatorHam(magnitude_spin);
      sx_  = CreateOperatorSx (magnitude_spin);
      isy_ = CreateOperatoriSy(magnitude_spin);
      sz_  = CreateOperatorSz (magnitude_spin);
      sp_  = CreateOperatorSp (magnitude_spin);
      sm_  = CreateOperatorSm (magnitude_spin);
   }
   
   CRS CreateOperatorHam(double magnitude_spin) const {
      return D_z_*CreateOperatorSzSz(magnitude_spin) + h_z_*CreateOperatorSz(magnitude_spin);
   }
      
};



} // namespace lattice
} // namespace compnal

#endif /* heisenberg_model_1d_hpp */
