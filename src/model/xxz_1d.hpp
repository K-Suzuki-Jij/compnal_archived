//
//  xxz_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/18.
//

#ifndef COMPNAL_MODEL_XXZ_1D_HPP_
#define COMPNAL_MODEL_XXZ_1D_HPP_

#include "./base_u1_spin_1d.hpp"

namespace compnal {
namespace model {

template<typename RealType>
class XXZ_1D: public BaseU1Spin_1D<RealType> {
   
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   XXZ_1D(): BaseU1Spin_1D<RealType>() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2spin_);
   }
   
   explicit XXZ_1D(const int system_size): BaseU1Spin_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2spin_);
   }
   
   XXZ_1D(const int system_size, const double magnitude_spin):
   BaseU1Spin_1D<RealType>(system_size, magnitude_spin) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2spin_);
   }
   
   XXZ_1D(const int system_size, const utility::BoundaryCondition bc):
   BaseU1Spin_1D<RealType>(system_size, bc) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2spin_);
   }
   
   XXZ_1D(const int system_size, const double magnitude_spin, const utility::BoundaryCondition bc):
   BaseU1Spin_1D<RealType>(system_size, magnitude_spin, bc) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2spin_);
   }
   
   void SetJz(const std::vector<RealType> &J_z) {
      if (J_z_ != J_z) {
         J_z_ = J_z;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetJz(const RealType J_z) {
      if (J_z_.size() == 0) {
         J_z_.push_back(J_z);
         this->calculated_eigenvector_set_.clear();
      }
      else if (J_z_[0] != J_z) {
         J_z_[0] = J_z;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetJxy(const std::vector<RealType> &J_xy) {
      if (J_xy_ != J_xy) {
         J_xy_ = J_xy;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetJxy(const RealType J_xy) {
      if (J_xy_.size() == 0) {
         J_xy_.push_back(J_xy);
         this->calculated_eigenvector_set_.clear();
      }
      else if (J_xy_[0] != J_xy) {
         J_xy_[0] = J_xy;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetHz(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(this->magnitude_2spin_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetDz(const RealType D_z) {
      if (D_z_ != D_z) {
         D_z_ = D_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(this->magnitude_2spin_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void PrintInfo() const {
      std::string bc = "None";
      if (this->boundary_condition_ == utility::BoundaryCondition::OBC) {
         bc = "OBC";
      }
      else if (this->boundary_condition_ == utility::BoundaryCondition::PBC) {
         bc = "PBC";
      }
      else if (this->boundary_condition_ == utility::BoundaryCondition::SSD) {
         bc = "SSD";
      }
      std::cout << "Print Heisenberg Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << this->bc                      << std::endl;
      std::cout << "system_size            = " << this->system_size_            << std::endl;
      std::cout << "magnitute_2spin        = " << this->magnitude_2spin_        << std::endl;
      std::cout << "total_2sz              = " << this->total_2sz_              << std::endl;
      std::cout << "dim_target             = " << this->CalculateTargetDim()    << std::endl;
      std::cout << "dim_onsite             = " << this->dim_onsite_             << std::endl;
      std::cout << "num_conserved_quantity = " << this->num_conserved_quantity_ << std::endl;
      
      std::cout << "Print Heisenberg Interaction" << std::endl;
      std::cout << "Sz-Sz Interaction: J_z =" << std::endl;
      for (std::int64_t i = 0; i < J_z_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_z_.at(i) << std::endl;
      }
      std::cout << "Sx-Sx, Sy-Sy Interactions: J_xy =" << std::endl;
      for (std::int64_t i = 0; i < J_xy_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << J_xy_.at(i) << std::endl;
      }
      std::cout << "External Magnetic Fields for the z-direction: h_z =" << h_z_ << std::endl;
      std::cout << "Uniaxial Anisotropy for the z-direction: D_z =" << D_z_ << std::endl;
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
   
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }
   
   inline const std::vector<RealType> &GetJz () const { return J_z_ ; }
   inline const std::vector<RealType> &GetJxy() const { return J_xy_; }
   
   inline RealType GetJz (const std::int64_t index) const { return J_z_ .at(index); }
   inline RealType GetJxy(const std::int64_t index) const { return J_xy_.at(index); }
   
   inline RealType GetHz() const { return h_z_; }
   inline RealType GetDz() const { return D_z_; }
   
   
private:
   CRS onsite_operator_ham_;
   
   std::vector<RealType> J_z_  = {1.0};
   std::vector<RealType> J_xy_ = {1.0};
   RealType h_z_ = 0.0;
   RealType D_z_ = 0.0;
   
   
};



}
}


#endif /* COMPNAL_MODEL_XXZ_1D_HPP_ */
