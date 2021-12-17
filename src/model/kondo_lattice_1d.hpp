//
//  kondo.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/17.
//

#ifndef COMPNAL_MODEL_KONDO_LATTICE_1D_HPP_
#define COMPNAL_MODEL_KONDO_LATTICE_1D_HPP_

#include "./base_u1_spin_electron_1d.hpp"

namespace compnal {
namespace model {

//! @brief The class for the one-dimensional Kondo lattice model with the magnitude of the spin \f$ S\f$.
//! The Hamiltonian reads
//! \f[ \hat{H}=\sum_{d}t_{d}\sum^{N-1}_{i=1}\left(\hat{c}^{\dagger}_{i,\sigma}\hat{c}_{i+d,\sigma}+{\rm h.c.}\right)+
//! \sum^{N}_{i=1}\left(J_{xy}\hat{s}^{x}_{i}\hat{S}^{x}_{i}+J_{xy}\hat{s}^{y}_{i}\hat{S}^{y}_{i}+J_{z}\hat{s}^{z}_{i}\hat{S}^{z}_{i}\right) +
//! h_z\sum^{N}_{i=1}\left(\hat{s}^{x}_{i} + \hat{S}^{x}_{i}\right) + D_z\sum^{N}_{i=1}\left(\hat{S}^{z}_{i}\right)^2 \f]
//! @tparam RealType The type of real values.
template<typename RealType>
class KondoLattice_1D: public BaseU1SpinElectron_1D<RealType> {

   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of KondoLattice_1D.
   KondoLattice_1D(): BaseU1SpinElectron_1D<RealType>() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   KondoLattice_1D(const int system_size): BaseU1SpinElectron_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   KondoLattice_1D(const int system_size,
                   const double magnitude_lspin):
   BaseU1SpinElectron_1D<RealType>(system_size, magnitude_lspin) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   KondoLattice_1D(const int system_size,
                   const int total_electron):
   BaseU1SpinElectron_1D<RealType>(system_size, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   KondoLattice_1D(const int system_size,
                   const double magnitude_lspin,
                   const int total_electron):
   BaseU1SpinElectron_1D<RealType>(system_size, magnitude_lspin, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param boundary_condition Boundary condition.
   KondoLattice_1D(const int system_size,
                   const utility::BoundaryCondition boundary_condition):
   BaseU1SpinElectron_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
      SetBoundaryCondition(boundary_condition);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param boundary_condition Boundary condition.
   KondoLattice_1D(const int system_size,
                   const double magnitude_lspin,
                   const utility::BoundaryCondition boundary_condition):
   BaseU1SpinElectron_1D<RealType>(system_size, magnitude_lspin) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
      SetBoundaryCondition(boundary_condition);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   //! @param boundary_condition Boundary condition.
   KondoLattice_1D(const int system_size,
                   const int total_electron,
                   const utility::BoundaryCondition boundary_condition):
   BaseU1SpinElectron_1D<RealType>(system_size, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
      SetBoundaryCondition(boundary_condition);
   }
   
   //! @brief Constructor of KondoLattice_1D.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param total_electron The number of the total electrons
   //! \f$ \langle \hat{N}_{e}\rangle =\sum^{N}_{i=1}\langle\hat{n}_{i}\rangle\f$.
   //! @param boundary_condition Boundary condition.
   KondoLattice_1D(const int system_size,
                   const double magnitude_lspin,
                   const int total_electron,
                   const utility::BoundaryCondition boundary_condition):
   BaseU1SpinElectron_1D<RealType>(system_size, magnitude_lspin, total_electron) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
      SetBoundaryCondition(boundary_condition);
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
   
   //! @brief Set the magnetic fields for the z-direction.
   //! @param h_z The magnetic fields for the z-direction.
   void SetMagneticField(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin, J_z_, J_xy_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @param D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   void SetDz(const RealType D_z) {
      if (D_z_ != D_z) {
         D_z_ = D_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the Kondo exchange coupling along the z-direction \f$J_z \f$.
   //! @param J_z The Kondo exchange coupling along the z-direction \f$J_z \f$.
   void SetJz(const RealType J_z) {
      if (J_z_ != J_z) {
         J_z_ = J_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }
   
   //! @brief Set the Kondo exchange coupling along the x, y-direction \f$J_{xy} \f$.
   //! @param J_xy The Kondo exchange coupling along the x, y-direction \f$J_{xy} \f$.
   void SetJxy(const RealType J_xy) {
      if (J_xy_ != J_xy) {
         J_xy_ = J_xy;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5*this->magnitude_2lspin_, J_z_, J_xy_, h_z_, D_z_);
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
      std::cout << "Print Kondo Lattice Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << this->boundary_condition_     << std::endl;
      std::cout << "system_size            = " << this->system_size_            << std::endl;
      std::cout << "magnitute_2lspin       = " << this->magnitude_2lspin_       << std::endl;
      std::cout << "total_2sz              = " << this->total_2sz_              << std::endl;
      std::cout << "dim_target             = " << this->CalculateTargetDim()    << std::endl;
      std::cout << "dim_onsite             = " << this->dim_onsite_             << std::endl;
      std::cout << "num_conserved_quantity = " << this->num_conserved_quantity_ << std::endl;
      
      std::cout << "Print Interactions" << std::endl;
      std::cout << "Electron Hopping: t =" << std::endl;
      for (std::size_t i = 0; i < t_.size(); ++i) {
         std::cout << i + 1 << "-th neighber: " << t_.at(i) << std::endl;
      }
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
   
   //! @brief Create the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=
   //! J_{xy}\left(\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}\right)+J_{z}\hat{s}^{z}\hat{S}^{z} +
   //! h_z\left(\hat{s}^{x} + \hat{S}^{x}\right) +
   //! D_z\left(\hat{S}^{z}\right)^2
   //! \f]
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param J_z The Kondo exchange coupling along the z-direction \f$J_z \f$.
   //! @param J_xy The Kondo exchange coupling along the x, y-direction \f$J_{xy} \f$.
   //! @param h_z The magnetic fields along the z-direction \f$ h_z \f$.
   //! @param D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   static CRS CreateOnsiteOperatorHam(const double magnitude_lspin, const RealType J_z, const RealType J_xy, const RealType h_z, const RealType D_z) {
      const CRS spc = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSpC(magnitude_lspin);
      const CRS smc = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSmC(magnitude_lspin);
      const CRS szc = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSzC(magnitude_lspin);
      const CRS spl = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSpL(magnitude_lspin);
      const CRS sml = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSmL(magnitude_lspin);
      const CRS szl = BaseU1SpinElectron_1D<RealType>::CreateOnsiteOperatorSzL(magnitude_lspin);
      return J_z*szc*szl + J_xy*0.5*(spc*sml + smc*spl) + h_z*(szc + szl) * D_z*szl*szl;
   }
   
   //! @brief Create the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=
   //! J\hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}} +
   //! h_z\left(\hat{s}^{x} + \hat{S}^{x}\right) +
   //! D_z\left(\hat{S}^{z}\right)^2
   //! \f]
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param J The Kondo exchange coupling.
   //! @param h_z The magnetic fields along the z-direction \f$ h_z \f$.
   //! @param D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   static CRS CreateOnsiteOperatorHam(const double magnitude_lspin, const RealType J, const RealType h_z, const RealType D_z) {
      return CreateOnsiteOperatorHam(magnitude_lspin, J, J, h_z, D_z);
   }

   //! @brief Get the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=
   //! J_{xy}\left(\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}\right)+J_{z}\hat{s}^{z}\hat{S}^{z} +
   //! h_z\left(\hat{s}^{x} + \hat{S}^{x}\right) +
   //! D_z\left(\hat{S}^{z}\right)^2
   //! \f]
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }
   
   //! @brief Get hopping energy \f$ t \f$.
   //! @return The hopping energy \f$ t \f$.
   inline const std::vector<RealType> &GetHopping() const { return t_ ; }
   
   //! @brief Get hopping energy \f$ t_{d} \f$ at the distance \f$ d \f$.
   //! @param index The distance \f$ d \f$.
   //! @return The hopping energy \f$ t_{d} \f$ at the distance \f$ d \f$.
   inline RealType GetHopping(const std::int64_t index) const { return t_.at(index); }
   
   //! @brief Get the spin-spin interaction along the z-direction \f$ J_{z} \f$.
   //! @return J_z The spin-spin interaction along the z-direction \f$ J_{z} \f$.
   inline RealType &GetJz() const { return J_z_ ; }
   
   //! @brief Get the spin-spin interaction along the x, y-direction \f$ J_{xy} \f$.
   //! @return J_xy The spin-spin interaction along the x, y-direction \f$ J_{xy} \f$.
   inline RealType &GetJxy() const { return J_xy_; }
      
   //! @brief Get the magnetic field for the z-direction \f$ h_z\f$.
   //! @return h_z The magnetic field for the z-direction \f$ h_z\f$.
   inline RealType GetMagneticField() const { return h_z_; }
   
   //! @brief Get the uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @return D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   inline RealType GetDz() const { return D_z_; }
   
   //! @brief Get the boundary condition.
   //! @return The boundary condition.
   inline utility::BoundaryCondition GetBoundaryCondition() const { return boundary_condition_; }
   
private:
   //! @brief The onsite Hamiltonian.
   //! \f$ \hat{H}_{\rm onsite}=
   //! J_{xy}\left(\hat{s}^{x}\hat{S}^{x}+\hat{s}^{y}\hat{S}^{y}\right)+J_{z}\hat{s}^{z}\hat{S}^{z} +
   //! h_z\left(\hat{s}^{x} + \hat{S}^{x}\right) +
   //! D_z\left(\hat{S}^{z}\right)^2
   //! \f$
   CRS onsite_operator_ham_;
   
   //! @brief Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   utility::BoundaryCondition boundary_condition_ = utility::BoundaryCondition::OBC;
   
   //! @brief Hopping energy \f$ t \f$.
   std::vector<RealType> t_ = {1.0};
   
   //! @brief The Kondo exchange coupling along the z-direction \f$J_z \f$.
   RealType J_z_ = 1.0;
   
   //! @brief The Kondo exchange coupling along the x, y-direction \f$J_{xy} \f$.
   RealType J_xy_ = 1.0;
   
   //! @brief The magnetic fields along the z-direction \f$ h_z \f$.
   RealType h_z_;
   
   //! @brief The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   RealType D_z_;
   
};

}
}


#endif /* COMPNAL_MODEL_KONDO_LATTICE_1D_HPP_ */
