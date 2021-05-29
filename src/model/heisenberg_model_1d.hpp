//
//  heisenberg_model.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef heisenberg_model_hpp
#define heisenberg_model_hpp

#include "model_utility.hpp"
#include "sparse_matrix.hpp"

#include <cmath>
#include <sstream>

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
      SetDim();
   }
   
   Heisenberg1D(const int system_size, const double magnitude_spin) {
      SetSystemSize(system_size);
      magnitude_2spin_ = DoubleTheNumber(magnitude_spin);
      dim_onsite_ = magnitude_2spin_ + 1;
      SetOperator();
      SetDim();
   }
   
   Heisenberg1D(const int system_size, const BoundaryCondition bc) {
      SetSystemSize(system_size);
      SetBoundaryCondition(bc);
      SetOperator();
      SetDim();
   }
   
   Heisenberg1D(const int system_size, const double magnitude_spin, const BoundaryCondition bc) {
      SetSystemSize(system_size);
      SetBoundaryCondition(bc);
      magnitude_2spin_ = DoubleTheNumber(magnitude_spin);
      dim_onsite_      = magnitude_2spin_ + 1;
      SetOperator();
      SetDim();
   }
      
   inline BoundaryCondition GetBoundaryCondition()    const { return boundary_condition_;     }
   inline int               GetSystemSize()           const { return system_size_;            }
   inline int               GetDimOnsite()            const { return dim_onsite_;             }
   inline double            GetMagnitudeSpin()        const { return magnitude_2spin_/2.0;    }
   inline double            GetTotalSz()              const { return total_2sz_/2.0;          }
   inline const int         GetNumConservedQuantity() const { return num_conserved_quantity_; }
   
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
      }
   }
   
   void SetTotal2Sz(const double total_sz) {
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
   }
   
   void SetSystemSize(const int system_size) {
      if (system_size <= 0) {
         std::stringstream ss;
         ss << "system_size must be more than 0" << std::endl;
         ss << "system_size=" << system_size << "is not allowed" << std::endl;
         throw std::runtime_error(ss.str());
      }
      system_size_ = system_size;
   }
   
   void SetBoundaryCondition(const BoundaryCondition bc) {
      boundary_condition_ = bc;
   }
   
   void SetJz (const std::vector<RealType> &J_z) {
      J_z_  = J_z;
   }
   void SetJz (const RealType J_z) {
      J_z_  = std::vector<RealType>{J_z};
   }
   
   template<typename... Args>
   void SetJz (Args... args) {
      J_z_ = std::vector<RealType>{args...};
   }
   
   void SetJxy (const std::vector<RealType> &J_xy) {
      J_xy_  = J_xy;
   }
   void SetJxy (const RealType J_xy) {
      J_xy_  = std::vector<RealType>{J_xy};
   }
   
   template<typename... Args>
   void SetJxy (Args... args) {
      J_xy_ = std::vector<RealType>{args...};
   }
   
   void SetHz(const RealType h_z) {
      h_z_ = h_z;
      ham_ = CreateHam();
   }
   
   void SetDz(const RealType D_z) {
      D_z_ = D_z;
      ham_ = CreateHam();
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
      std::cout << "num_conserved_quantity = " << num_conserved_quantity_ << std::endl;
      std::cout << "dim_onsite             = " << dim_onsite_             << std::endl;
      
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
   
   void PrintBasis() {
      const double magnitude_spin = magnitude_2spin_/2.0;
      for (int64_t row = 0; row < sz_.GetRowDim(); ++row) {
         std::cout << "row " << row << ": |Sz=" << magnitude_spin - row << ">" << std::endl;
      }
   }
   
   int DoubleTheNumber(double s) const {
      s = 2*s;
      if (std::floor(s) != s) {
         throw std::runtime_error("Invalid value of magnitude_spin or total_sz");
      }
      return static_cast<int>(s);
   }
   
   int64_t GetDim() const {
      int system_size   = GetSystemSize();
      int total_2sz     = total_2sz_;
      int max_total_2sz = system_size*magnitude_2spin_;
      if (dim_[system_size - 1][(total_2sz + max_total_2sz)/2] < 0) {
         throw std::runtime_error("Maybe too large Hilbert space. Overflow detected.");
      }
      return dim_[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   int64_t GetDim(double target_sz) const {
      int system_size   = GetSystemSize();
      int total_2sz     = DoubleTheNumber(target_sz);
      int max_total_2sz = system_size*magnitude_2spin_;
      if (dim_[system_size - 1][(total_2sz + max_total_2sz)/2] < 0) {
         throw std::runtime_error("Maybe too large Hilbert space. Overflow detected.");
      }
      if (std::abs(total_2sz) > max_total_2sz) {
         return 0;
      }
      return dim_[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   const CRS &GetOperatorHam() const { return ham_; }
   const CRS &GetOperatorSx () const { return sx_ ; }
   const CRS &GetOperatoriSy() const { return isy_; }
   const CRS &GetOperatorSz () const { return sz_ ; }
   const CRS &GetOperatorSp () const { return sp_ ; }
   const CRS &GetOperatorSm () const { return sm_ ; }
   
   const std::vector<RealType> &GetJz()  const { return J_z_ ; }
   const std::vector<RealType> &GetJxy() const { return J_xy_; }
   RealType GetHz() const { return h_z_; }
   RealType GetDz() const { return D_z_; }

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
   
   std::vector<std::vector<int64_t>> dim_;
   
   void SetOperator() {
      ham_ = CreateHam();
      sx_  = CreateSx ();
      isy_ = CreateiSy();
      sz_  = CreateSz ();
      sp_  = CreateSp ();
      sm_  = CreateSm ();
   }
   
   CRS CreateHam() const {
      return D_z_*CreateSzSz() + h_z_*CreateSz();
   }
   
   CRS CreateSp() const {
      return CreateSx() + CreateiSy();
   }
   
   CRS CreateSm() const {
      return CreateSx() - CreateiSy();
   }
   
   CRS CreateSzSz() const {
      return CreateSz()*CreateSz();
   }
   
   CRS CreateSx() const {
      int    dim_onsite     = GetDimOnsite();
      double magnitude_spin = GetMagnitudeSpin();
      RealType a = 0;
      RealType b = 1;
      
      CRS matrix(dim_onsite, dim_onsite);

      matrix.PushVal(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
      matrix.PushCol(b);
      matrix.Row(1) = matrix.GetSizeCol();
      
      for (int row = 1; row < dim_onsite - 1; ++row) {
         a = row;
         b = row - 1;
         matrix.PushVal(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
         matrix.PushCol(b);
         
         a = row;
         b = row + 1;
         matrix.PushVal(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
         matrix.PushCol(b);
         matrix.Row(row + 1) = matrix.GetSizeCol();
      }
      
      a = dim_onsite - 1;
      b = dim_onsite - 2;
      
      matrix.PushVal(0.5*std::sqrt((magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1)) );
      matrix.PushCol(b);
      matrix.Row(dim_onsite) = matrix.GetSizeCol();
      

      return matrix;
   }
   
   CRS CreateiSy() const {
      int   dim_onsite      = GetDimOnsite();
      double magnitude_spin = GetMagnitudeSpin();
      RealType a = 0;
      RealType b = 1;
      
      CRS matrix(dim_onsite, dim_onsite);
      
      matrix.PushVal(0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
      matrix.PushCol(b);
      matrix.Row(1) = matrix.GetSizeCol();
      
      for (int row = 1; row < dim_onsite - 1; ++row) {
         a = row;
         b = row - 1;
         matrix.PushVal(-0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
         matrix.PushCol(b);
         
         a = row;
         b = row + 1;
         matrix.PushVal(0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
         matrix.PushCol(b);
         
         matrix.Row(row + 1) = matrix.GetSizeCol();
      }
      
      a = dim_onsite - 1;
      b = dim_onsite - 2;
      
      matrix.PushVal(-0.5*std::sqrt( (magnitude_spin + 1)*(a + b + 1) - (a + 1)*(b + 1) ) );
      matrix.PushCol(b);
      matrix.Row(dim_onsite) = matrix.GetSizeCol();
      
      return matrix;
   }

   CRS CreateSz() const {
      int   dim_onsite      = GetDimOnsite();
      double magnitude_spin = GetMagnitudeSpin();

      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         RealType val = magnitude_spin - row;
         if (val != 0.0) {
            matrix.PushVal(val);
            matrix.PushCol(row);
         }
         matrix.Row(row + 1) = matrix.GetSizeCol();
      }
      
      return matrix;
   }
   
   void SetDim() {
      int system_size = GetSystemSize();
      int magnitude_2spin = magnitude_2spin_;
      int max_total_2sz = system_size*magnitude_2spin;
      dim_.assign(system_size, std::vector<int64_t>(max_total_2sz + 1));
      
      for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
         dim_[0][(s + magnitude_2spin)/2] = 1;
      }
      
      for (int site = 1; site < system_size; site++) {
         for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
            for (int s_prev = -magnitude_2spin*site; s_prev <= magnitude_2spin*site; s_prev += 2) {
               dim_[site][(s + s_prev + magnitude_2spin*(site + 1))/2] += dim_[site - 1][(s_prev + magnitude_2spin*site)/2];
            }
         }
      }
   }
      
};



} // namespace lattice
} // namespace compnal

#endif /* heisenberg_model_hpp */
