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

//! @brief The class for the one-dimensional Hubbard model.
//! The Hamiltonian reads
//! \f[ \hat{H}=\sum_{d}t_{d}\sum^{N}_{i=1}\sum_{\sigma=\uparrow,\downarrow}
//! \left(\hat{c}^{\dagger}_{i,\sigma}\hat{c}_{i+d,\sigma}+{\rm h.c.}\right)+
//! U\sum^{N}_{i=1}\hat{n}_{i,\uparrow}\hat{n}_{i,\downarrow}+
//! \sum_{d}V_{d}\sum^{N}_{i=1}\hat{n}_{i}\hat{n}_{i+d}+
//! h_z\sum^{N}_{i=1}s^{z}_{i} + \mu\sum^{N}_{i=1}\hat{n}_{i} \f]
//! @tparam RealType The type of real values.
template<typename RealType>
class Hubbard_1D: public BaseU1Electron_1D<RealType> {
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of Hubbard_1D class.
   Hubbard_1D(): BaseU1Electron_1D<RealType>() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   //! @brief Constructor of Hubbard_1D class.
   //! @param system_size The system size \f$ N \f$.
   explicit Hubbard_1D(const int system_size): BaseU1Electron_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   //! @brief Constructor of Hubbard_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   Hubbard_1D(const int system_size, const int total_electron): BaseU1Electron_1D<RealType>(system_size, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   //! @brief Constructor of Hubbard_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param boundary_condition Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   Hubbard_1D(const int system_size, const utility::BoundaryCondition boundary_condition): BaseU1Electron_1D<RealType>(system_size) {
      SetBoundaryCondition(boundary_condition);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   //! @brief Constructor of Hubbard_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   //! @param boundary_condition Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   Hubbard_1D(const int system_size, const int total_electron,
              const utility::BoundaryCondition boundary_condition): BaseU1Electron_1D<RealType>(system_size, total_electron) {
      SetBoundaryCondition(boundary_condition);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_);
   }
   
   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------
   //! @brief Set the boundary condition.
   //! @param boundary_condition Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   void SetBoundaryCondition(const utility::BoundaryCondition boundary_condition) {
      boundary_condition_ = boundary_condition;
   }
   
   //! @brief Set hopping energy \f$ t \f$.
   //! @param t The hopping energy.
   void SetHopping(const std::vector<RealType> &t) {
      if (t_ != t) {
         t_ = t;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set hopping energy \f$ t \f$.
   //! @param t The hopping energy.
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
   
   //! @brief Set the intersite density-density interactions.
   //! @param V The intersite density-density interactions.
   void SetIntersiteCoulomb(const std::vector<RealType> &V) {
      if (V_ != V) {
         V_ = V;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the intersite density-density interactions.
   //! @param V The intersite density-density interactions.
   void SetIntersiteCoulomb(const RealType V) {
      if (V_.size() == 0) {
         V_.push_back(V);
         this->calculated_eigenvector_set_.clear();
      }
      else if (V_[0] != V) {
         V_[0] = V;
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the onsite density interactions.
   //! @param U The onsite density interactions.
   void SetOnsiteCoulomb(const RealType U) {
      if (U_ != U) {
         U_ = U;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_, U_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the magnetic fields for the z-direction.
   //! @param h_z The magnetic fields for the z-direction.
   void SetMagneticField(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_, U_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the chemical potential.
   //! @param mu The chemical potential.
   void SetChemicalPotential(const RealType mu) {
      if (mu_ != mu) {
         mu_ = mu;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(h_z_, U_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Print information about this class.
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
   
   //! @brief Create the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{s}^{z}+U\hat{n}_{\uparrow}\hat{n}_{\downarrow}+\mu(\hat{n}_{\uparrow}+\hat{n}_{\downarrow})\f]
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   static CRS CreateOnsiteOperatorHam(const RealType h_z = 0.0, const RealType U = 0.0, const RealType mu = 0.0) {
      const CRS m_s_z    = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorSz();
      const CRS m_n_up   = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorNCUp();
      const CRS m_n_down = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorNCDown();
      const CRS m_mu     = BaseU1Electron_1D<RealType>::CreateOnsiteOperatorNC();
      return h_z*m_s_z + U*m_n_up*m_n_down + mu*m_mu;
   }
   
   //! @brief Get the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{s}^{z}+U\hat{n}_{\uparrow}\hat{n}_{\downarrow}+\mu(\hat{n}_{\uparrow}+\hat{n}_{\downarrow})\f]
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }
   
   //! @brief Get hopping energy \f$ t \f$.
   //! @return The hopping energy \f$ t \f$.
   inline const std::vector<RealType> &GetHopping() const { return t_ ; }
   
   //! @brief Get intersite density-density interactions. \f$ V \f$.
   //! @return The intersite density-density interactions. \f$ V \f$.
   inline const std::vector<RealType> &GetIntersiteCoulomb() const { return V_; }
   
   //! @brief Get hopping energy \f$ t_{d} \f$ at the distance \f$ d \f$.
   //! @param index The distance \f$ d \f$.
   //! @return The hopping energy \f$ t_{d} \f$ at the distance \f$ d \f$.
   inline RealType GetHopping(const std::int64_t index) const { return t_.at(index); }
   
   //! @brief Get intersite density-density interactions. \f$ V_{d} \f$ at the distance \f$ d \f$.
   //! @return The intersite density-density interactions. \f$ V_{d} \f$ at the distance \f$ d \f$.
   inline RealType GetIntersiteCoulomb(const std::int64_t index) const { return V_.at(index); }
   
   //! @brief Get onsite density interactions. \f$ U \f$.
   //! @return The onsite density interactions. \f$ U \f$.
   inline RealType GetOnsiteCoulomb() const { return U_; }
   
   //! @brief Get the magnetic fields along the z-direction. \f$ h_z \f$.
   //! @return The magnetic fields along the z-direction. \f$ h_z \f$.
   inline RealType GetMagneticField() const { return h_z_; }
   
   //! @brief Get the chemical potential \f$ mu \f$.
   //! @return The chemical potential \f$ mu \f$.
   inline RealType GetChemicalPotential() const { return mu_; }
   
   //! @brief Get the boundary condition.
   //! @return The boundary condition.
   inline utility::BoundaryCondition GetBoundaryCondition() const { return boundary_condition_; }

private:
   //! @brief the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{s}^{z}+U\hat{n}_{\uparrow}\hat{n}_{\downarrow}+\mu(\hat{n}_{\uparrow}+\hat{n}_{\downarrow})\f]
   CRS onsite_operator_ham_;

   //! @brief Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   utility::BoundaryCondition boundary_condition_ = utility::BoundaryCondition::OBC;

   //! @brief Hopping energy \f$ t \f$.
   std::vector<RealType> t_ = {1.0};
   
   //! @brief The intersite density-density interactions \f$ V \f$.
   std::vector<RealType> V_ = {};
   
   //! @brief The onsite density interactions \f$ U \f$.
   RealType U_   = 1.0;
   
   //! @brief The magnetic fields along the z-direction \f$ h_z \f$.
   RealType h_z_ = 0.0;
   
   //! @brief The chemical potential \f$ mu \f$.
   RealType mu_  = 0.0;
   
};

}
}

#endif /* COMPNAL_MODEL_HUBBARD_1D_HPP_ */
