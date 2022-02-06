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

#ifndef COMPNAL_MODEL_U1_ELECTRON_1D_HPP_
#define COMPNAL_MODEL_U1_ELECTRON_1D_HPP_

#include "../sparse_matrix/all.hpp"
#include "../utility/all.hpp"

#include <unordered_map>
#include <vector>

namespace compnal {
namespace model {

//! @brief The general model class for one-dimensional systems.
//! @tparam BaseClass The base model classes.
template<class BaseClass>
class GeneralModel_1D: public BaseClass {
   
   //! @brief The type of real values.
   using RealType = typename BaseClass::ValueType;
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   using BaseClass::BaseClass;
   
   //! @brief Add an onsite potential term to the system \f$ c\times \hat{O}_{i} \f$.
   //! @param m The onsite operator \f$ \hat{O}_{i} \f$.
   //! @param site The site \f$ i\f$.
   //! @param val The coefficient of the potential term \f$ c \f$.
   void AddOnsitePotential(const CRS &m, const int site, const RealType val = 1.0) {
      if (val == 0.0) {
         return;
      }
      // TODO: Check of m is valid with BaseClass operator
      if (m.row_dim != m.col_dim) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "The input matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (m.row_dim != this->dim_onsite_) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "The input matrix is invalid" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site < 0 || this->system_size_ <= site) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Invalid site: " << site << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (onsite_operator_list_.count(site) == 0) {
         onsite_operator_list_[site] = m.MultiplyByScalar(val);
      }
      else {
         onsite_operator_list_.at(site) = CalculateMatrixMatrixSum(val, m, 1.0, onsite_operator_list_.at(site));
      }
   }
      
   //! @brief Add an interaction term to the system \f$ c\times\hat{O}_{i}\hat{P}_{j} \f$.
   //! @param m_1 The onsite operator \f$ \hat{O}_{i} \f$.
   //! @param site_1 The site \f$ i\f$.
   //! @param m_2 The onsite operator \f$ \hat{P}_{j} \f$.
   //! @param site_2 The site \f$ j\f$.
   //! @param val The coefficient of the interaction term \f$ c \f$.
   void AddInteraction(const CRS &m_1, const int site_1, const CRS &m_2, const int site_2, const RealType val = 1.0) {
      if (val == 0.0) {
         return;
      }
      if (m_1.row_dim != m_1.col_dim || m_2.row_dim != m_2.col_dim) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "The input matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (m_1.row_dim != this->dim_onsite_) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "The input matrix is invalid" << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site_1 < 0 || this->system_size_ <= site_1) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Invalid site: " << site_1 << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site_2 < 0 || this->system_size_ <= site_2) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
         ss << "Invalid site: " << site_2 << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site_1 == site_2) {
         AddOnsitePotential(m_1*m_2, site_1, val);
         return;
      }
      intersite_operator_list_[site_1][site_2].push_back({m_1.MultiplyByScalar(val), m_2});
   }
   
   //! @brief Get onsite potential term list added to the system.
   //! @return Onsite potential term list added to the system.
   inline const std::unordered_map<int, CRS> &GetOnsiteOperatorList() const { return onsite_operator_list_; }
   
   //! @brief Get an onsite potential term added to the system \f$ \hat{O}_{i}\f$.
   //! @param site The site \f$ i\f$.
   //! @return The onsite potential term added to the system \f$ \hat{O}_{i}\f$.
   inline const CRS &GetOnsiteOperator(const int site) const { return onsite_operator_list_.at(site); }
   
   //! @brief Get interaction term list added to the system.
   //! @return Interaction term list added to the system.
   inline const std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<CRS, CRS>>>> &GetIntersiteOperatorList() const { return intersite_operator_list_; }
   
   //! @brief Get an interaction term added to the system \f$ \hat{O}_{i}\hat{P}_{j}\f$.
   //! @param site_1 The site \f$ i\f$.
   //! @param site_2 The site \f$ j\f$.
   //! @return The interaction term added to the system \f$ \hat{O}_{i}\f$ and \f$ \hat{P}_{j}\f$.
   inline const std::vector<std::pair<CRS, CRS>> &GetIntersiteOperator(const int site_1, const int site_2) const { return intersite_operator_list_.at(site_1).at(site_2); }
   
private:
   //! @brief Onsite potential term list added to the system.
   std::unordered_map<int, CRS> onsite_operator_list_;
   
   //! @brief Interaction term list added to the system.
   std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<CRS, CRS>>>> intersite_operator_list_;
   
};

}
}





#endif /* COMPNAL_MODEL_U1_ELECTRON_1D_HPP_ */
