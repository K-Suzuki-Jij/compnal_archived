//
//  hubbard_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/04.
//

#ifndef COMPNAL_MODEL_HUBBARD_1D_HPP_
#define COMPNAL_MODEL_HUBBARD_1D_HPP_

#include "./base_u1_electron_1d.hpp"

namespace compnal {
namespace model {

template<typename RealType>
class Hubbard_1D: public BaseU1Electron_1D<RealType> {
   
   using CRS = sparse_matrix::CRS<RealType>;
   

public:
   Hubbard_1D(): BaseU1Electron_1D<RealType>() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   explicit Hubbard_1D(const int system_size): BaseU1Electron_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   Hubbard_1D(const int system_size, const int total_electron): BaseU1Electron_1D<RealType>(system_size, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   Hubbard_1D(const int system_size, const utility::BoundaryCondition bc): BaseU1Electron_1D<RealType>(system_size) {
      SetBoundaryCondition(bc);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   Hubbard_1D(const int system_size, const int total_electron, const utility::BoundaryCondition bc): BaseU1Electron_1D<RealType>(system_size, total_electron) {
      SetBoundaryCondition(bc);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   void SetBoundaryCondition(const utility::BoundaryCondition bc) {
      boundary_condition_ = bc;
   }
   
   void SetHopping(const std::vector<RealType> &t) {
      if (t_ != t) {
         t_ = t;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetHopping(const RealType t) {
      if (t_.size() == 0) {
         t_.push_back(t);
         this->calculated_eigenvector_set_.clear();
      }
      else if (t_[0] != t) {
         t_[0] = t;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetIntersiteCoulomb(const std::vector<RealType> &V) {
      if (V_ != V) {
         V_ = V;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetIntersiteCoulomb(const RealType &V) {
      if (V_.size() == 0) {
         V_.push_back(V);
         this->calculated_eigenvector_set_.clear();
      }
      else if (V_[0] != V) {
         V_[0] = V;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetOnsiteCoulomb(const RealType U) {
      if (U_ != U) {
         U_ = U;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_, U_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   void SetMagneticField(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_, U_);
         this->calculated_eigenvector_set_.clear();
      }
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
      std::cout << "Print Hubbard Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << this->boundary_condition_     << std::endl;
      std::cout << "system_size            = " << this->system_size_            << std::endl;
      std::cout << "total_sz               = " << 0.5*this->total_2sz_          << std::endl;
      std::cout << "dim_target             = " << this->CalculateTargetDim()    << std::endl;
      std::cout << "dim_onsite             = " << this->dim_onsite_             << std::endl;
      std::cout << "num_conserved_quantity = " << this->num_conserved_quantity_ << std::endl;
      
      std::cout << "Print Interaction" << std::endl;
      std::cout << "Electron Hopping: t =" << std::endl;
      for (std::size_t i = 0; i < t_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << t_.at(i) << std::endl;
      }
      std::cout << "Onsite Coulomb: U =" << U_ << std::endl;
      std::cout << "Magnetic Field for the z-direction: h_z =" << h_z_ << std::endl;
      std::cout << "Intersite Coulomb: V =" << std::endl;
      for (std::size_t i = 0; i < V_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << V_.at(i) << std::endl;
      }
   }
   
   static CRS CreateOnsiteOperatorHam(const RealType h_z = 0.0, const RealType U = 0.0) {
      const CRS m_s_z       = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorSz();
      const CRS m_n_up      = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorNCUp();
      const CRS m_n_down    = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorNCDown();
      const CRS m_n_up_down = sparse_matrix::CalculateMatrixMatrixProduct(1.0, m_n_up, 1.0, m_n_down);
      return sparse_matrix::CalculateMatrixMatrixSum(U, m_n_up_down, h_z, m_s_z);
   }
   
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }
   
   inline const std::vector<RealType> &GetHopping         () const { return t_ ; }
   inline const std::vector<RealType> &GetIntersiteCoulomb() const { return V_; }
   
   inline RealType GetHopping         (const std::int64_t index) const { return t_.at(index); }
   inline RealType GetIntersiteCoulomb(const std::int64_t index) const { return V_.at(index); }
   
   inline RealType GetOnsiteCoulomb() const { return U_; }
   inline RealType GetMagneticField() const { return h_z_; }
   
   inline utility::BoundaryCondition GetBoundaryCondition() const { return boundary_condition_; }

private:
   CRS onsite_operator_ham_;

   utility::BoundaryCondition boundary_condition_ = utility::BoundaryCondition::OBC;

   std::vector<RealType> t_ = {1.0};
   std::vector<RealType> V_ = {};
   RealType U_   = 1.0;
   RealType h_z_ = 0.0;
   
};

}
}

#endif /* COMPNAL_MODEL_HUBBARD_1D_HPP_ */
