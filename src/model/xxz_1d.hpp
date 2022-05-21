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

#ifndef COMPNAL_MODEL_XXZ_1D_HPP_
#define COMPNAL_MODEL_XXZ_1D_HPP_

#include "./base_u1_spin_1d.hpp"
#include "./utility.hpp"

namespace compnal {
namespace model {

//! @brief The class for the one-dimensional XXZ model with the magnitude of the spin \f$ S\f$.
//! The Hamiltonian reads
//! \f[ \hat{H}=\sum_{d}J_{d}\sum^{N}_{i=1}\left(J^{xy}_{d}\hat{S}^{x}_{i}\hat{S}^{x}_{i+d}+
//! J^{xy}_{d}\hat{S}^{y}_{i}\hat{S}^{y}_{i+d}+
//! J^{z}_{d}\hat{S}^{z}_{i}\hat{S}^{z}_{i+d}\right) \f]
//! @tparam RealType The type of real values.
template <typename RealType>
class XXZ_1D : public BaseU1Spin_1D<RealType> {
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = type::CRS<RealType>;

  public:
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of XXZ_1D class.
   XXZ_1D() : BaseU1Spin_1D<RealType>() {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5 * this->magnitude_2spin_, h_z_, D_z_);
   }

   //! @brief Constructor of XXZ_1D class.
   //! @param system_size The system size \f$ N \f$.
   explicit XXZ_1D(const int system_size) : BaseU1Spin_1D<RealType>(system_size) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5 * this->magnitude_2spin_, h_z_, D_z_);
   }

   //! @brief Constructor of XXZ_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   XXZ_1D(const int system_size, const double magnitude_spin) : BaseU1Spin_1D<RealType>(system_size, magnitude_spin) {
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5 * this->magnitude_2spin_, h_z_, D_z_);
   }

   //! @brief Constructor of XXZ_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param boundary_condition Boundary condition.
   XXZ_1D(const int system_size, const BoundaryCondition boundary_condition) : BaseU1Spin_1D<RealType>(system_size) {
      SetBoundaryCondition(boundary_condition);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5 * this->magnitude_2spin_, h_z_, D_z_);
   }

   //! @brief Constructor of XXZ_1D class.
   //! @param system_size The system size \f$ N \f$.
   //! @param magnitude_spin The magnitude of the spin \f$ S \f$.
   //! @param boundary_condition Boundary condition.
   XXZ_1D(const int system_size, const double magnitude_spin, const BoundaryCondition boundary_condition)
       : BaseU1Spin_1D<RealType>(system_size, magnitude_spin) {
      SetBoundaryCondition(boundary_condition);
      onsite_operator_ham_ = CreateOnsiteOperatorHam(0.5 * this->magnitude_2spin_, h_z_, D_z_);
   }

   //------------------------------------------------------------------
   //----------------------Public Member functions---------------------
   //------------------------------------------------------------------
   //! @brief Set the boundary condition.
   //! @param boundary_condition Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   void SetBoundaryCondition(const BoundaryCondition boundary_condition) { boundary_condition_ = boundary_condition; }

   //! @brief Set the spin-spin interaction along the z-direction \f$ J^{z} \f$.
   //! @param J_z The spin-spin interaction along the z-direction \f$ J^{z} \f$.
   void SetJz(const std::vector<RealType> &J_z) {
      if (J_z_ != J_z) {
         J_z_ = J_z;
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Set the spin-spin interaction along the z-direction \f$ J^{z}_{0} \f$.
   //! @param J_z The spin-spin interaction along the z-direction \f$ J^{z}_{0} \f$.
   void SetJz(const RealType J_z) {
      if (J_z_.size() == 0) {
         J_z_.push_back(J_z);
         this->calculated_eigenvector_set_.clear();
      } else if (J_z_[0] != J_z) {
         J_z_[0] = J_z;
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Set the spin-spin interaction along the x, y-direction \f$ J^{xy} \f$.
   //! @param J_xy The spin-spin interaction along the x, y-direction \f$ J^{xy} \f$.
   void SetJxy(const std::vector<RealType> &J_xy) {
      if (J_xy_ != J_xy) {
         J_xy_ = J_xy;
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Set the spin-spin interaction along the x, y-direction \f$ J^{xy}_{0} \f$.
   //! @param J_xy The spin-spin interaction along the x, y-direction \f$ J^{xy}_{0} \f$.
   void SetJxy(const RealType J_xy) {
      if (J_xy_.size() == 0) {
         J_xy_.push_back(J_xy);
         this->calculated_eigenvector_set_.clear();
      } else if (J_xy_[0] != J_xy) {
         J_xy_[0] = J_xy;
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Set the magnetic field for the z-direction \f$ h_z\f$.
   //! @param h_z The magnetic field for the z-direction \f$ h_z\f$.
   void SetHz(const RealType h_z) {
      if (h_z_ != h_z) {
         h_z_ = h_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(this->magnitude_2spin_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Set the uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @param D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   void SetDz(const RealType D_z) {
      if (D_z_ != D_z) {
         D_z_ = D_z;
         onsite_operator_ham_ = CreateOnsiteOperatorHam(this->magnitude_2spin_, h_z_, D_z_);
         this->calculated_eigenvector_set_.clear();
      }
   }

   //! @brief Print information about this class.
   void PrintInfo() const {
      std::string bc = "None";
      if (boundary_condition_ == BoundaryCondition::OBC) {
         bc = "OBC";
      } else if (boundary_condition_ == BoundaryCondition::PBC) {
         bc = "PBC";
      } else if (boundary_condition_ == BoundaryCondition::SSD) {
         bc = "SSD";
      }
      std::cout << "Print Heisenberg Model Infomation:" << std::endl;
      std::cout << "boundary_condition     = " << this->boundary_condition_ << std::endl;
      std::cout << "system_size            = " << this->system_size_ << std::endl;
      std::cout << "magnitute_2spin        = " << this->magnitude_2spin_ << std::endl;
      std::cout << "total_2sz              = " << this->total_2sz_ << std::endl;
      std::cout << "dim_target             = " << this->CalculateTargetDim() << std::endl;
      std::cout << "dim_onsite             = " << this->dim_onsite_ << std::endl;

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

   //! @brief Create the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{S}^{z}+D_z\left(\hat{S}^{z}\right)^{2}\f]
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   static CRS CreateOnsiteOperatorHam(const double magnitude_spin, const RealType h_z = 0.0, const RealType D_z = 0.0) {
      const int magnitude_2spin = utility::DoubleHalfInteger(magnitude_spin);
      const int dim_onsite = magnitude_2spin + 1;
      CRS matrix(dim_onsite, dim_onsite);

      for (int row = 0; row < dim_onsite; ++row) {
         const RealType val = h_z * (magnitude_spin - row) + D_z * D_z * (magnitude_spin - row);
         if (val != 0.0) {
            matrix.val.push_back(val);
            matrix.col.push_back(row);
         }
         matrix.row[row + 1] = matrix.col.size();
      }
      return matrix;
   }

   //! @brief Get the onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{S}^{z}+D_z\left(\hat{S}^{z}\right)^{2}\f]
   //! @return The matrix of \f$ \hat{H}_{\rm onsite}\f$.
   inline const CRS &GetOnsiteOperatorHam() const { return onsite_operator_ham_; }

   //! @brief Get the spin-spin interaction along the z-direction \f$ J^{z} \f$.
   //! @return J_z The spin-spin interaction along the z-direction \f$ J^{z} \f$.
   inline const std::vector<RealType> &GetJz() const { return J_z_; }

   //! @brief Get the spin-spin interaction along the x, y-direction \f$ J^{xy} \f$.
   //! @return J_xy The spin-spin interaction along the x, y-direction \f$ J^{xy} \f$.
   inline const std::vector<RealType> &GetJxy() const { return J_xy_; }

   //! @brief Get the spin-spin interaction along the z-direction \f$ J^{z}_{d} \f$.
   //! @param index The distance \f$ d\f$.
   //! @return J_z The spin-spin interaction along the z-direction \f$ J^{z}_{d} \f$.
   inline RealType GetJz(const std::int64_t index) const { return J_z_.at(index); }

   //! @brief Get the spin-spin interaction along the x, y-direction \f$ J^{xy}_{d} \f$.
   //! @param index The distance \f$ d\f$.
   //! @return J_xy The spin-spin interaction along the x, y-direction \f$ J^{xy}_{d} \f$.
   inline RealType GetJxy(const std::int64_t index) const { return J_xy_.at(index); }

   //! @brief Get the magnetic field for the z-direction \f$ h_z\f$.
   //! @return h_z The magnetic field for the z-direction \f$ h_z\f$.
   inline RealType GetHz() const { return h_z_; }

   //! @brief Get the uniaxial anisotropy to the z-direction \f$ D_z\f$.
   //! @return D_z The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   inline RealType GetDz() const { return D_z_; }

   //! @brief Get the boundary condition.
   //! @return The boundary condition.
   inline BoundaryCondition GetBoundaryCondition() const { return boundary_condition_; }

  private:
   //! @brief The onsite Hamiltonian.
   //! \f[ \hat{H}_{\rm onsite}=h_z\hat{S}^{z}+D_z\left(\hat{S}^{z}\right)^{2}\f]
   CRS onsite_operator_ham_;

   //! @brief Boundary condition.
   //! Open boundary condition (OBC), periodic boundary condition (PBC), or sine square deformation (SSD).
   BoundaryCondition boundary_condition_ = BoundaryCondition::OBC;

   //! @brief The spin-spin interaction along the z-direction \f$ J^{z} \f$.
   std::vector<RealType> J_z_ = {1.0};

   //! @brief The spin-spin interaction along the x, y-direction \f$ J^{xy} \f$.
   std::vector<RealType> J_xy_ = {1.0};

   //! @brief The magnetic field for the z-direction \f$ h_z\f$.
   RealType h_z_ = 0.0;

   //! @brief The uniaxial anisotropy to the z-direction \f$ D_z\f$.
   RealType D_z_ = 0.0;
};

}  // namespace model
}  // namespace compnal

#endif /* COMPNAL_MODEL_XXZ_1D_HPP_ */
