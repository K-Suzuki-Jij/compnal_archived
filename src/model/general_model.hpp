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

#ifndef COMPNAL_MODEL_GENERAL_MODEL_HPP_
#define COMPNAL_MODEL_GENERAL_MODEL_HPP_

#include "../blas/all.hpp"
#include "../utility/all.hpp"

#include <variant>
#include <unordered_set>

namespace compnal {
namespace model {

//! @brief The general model class.
//! @tparam BaseClass The base model classes.
template<class BaseClass>
class GeneralModel: public BaseClass {
   
   //------------------------------------------------------------------
   //------------------------Private Type Alias------------------------
   //------------------------------------------------------------------
   //! @brief Alias of ValueType.
   using RealType = typename BaseClass::ValueType;
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = blas::CRS<RealType>;
   
public:
   //------------------------------------------------------------------
   //------------------------Public Type Alias-------------------------
   //------------------------------------------------------------------
   //! @brief Alias of integer type.
   using IntegerType = std::int64_t;
   
   //! @brief Alias of string type.
   using StringType = std::string;
   
   //! @brief Alias of variant type as mixture of IntegerType and StringType.
   using VariantType = std::variant<IntegerType, StringType>;
   
   //! @brief Alias of hash struct for IndexType.
   using IndexHash = utility::VariantHash<std::int64_t>;
   
   //! @brief Alias of index type.
   using IndexType = std::variant<IntegerType, StringType, std::vector<VariantType>>;
   
   //! @brief Type for onsite operator lists.
   using OnsiteListType = std::unordered_map<IndexType, CRS, IndexHash>;
   
   //! @brief Type for lists of onsite operator pairs.
   using IntersiteListType = std::unordered_map<IndexType, std::unordered_map<IndexType, std::vector<std::pair<CRS, CRS>>, IndexHash>, IndexHash>;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   using BaseClass::BaseClass;
      
   //------------------------------------------------------------------
   //----------------------Public Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Add an onsite potential term \f$ c\hat{O}_{i} \f$ to the system.
   //! @param site The site index \f$ i \f$.
   //! @param m The onsite operator \f$ \hat{O}_{i} \f$.
   //! @param val The coefficient of the potential term \f$ c \f$.
   void AddPotential(const IndexType &site, const CRS &m, const RealType val = 1.0) {
      
      if (val == 0.0) {
         return;
      }

      if (m.row_dim != m.col_dim) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __FUNCTION__ << " in "<< __FILE__ << std::endl;
         ss << "The input matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
      
      if (m.row_dim != this->GetDimOnsite()) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __FUNCTION__ << " in "<< __FILE__ << std::endl;
         ss << "The input matrix is invalid." << std::endl;
         ss << "The dimension of input matrix is " << m.row_dim << std::endl;
         ss << "But the dimensiotn of this model is " << this->GetDimOnsite() << std::endl;
         m.PrintInfo();
         throw std::runtime_error(ss.str());
      }

      if (potential_list_.count(site) == 0) {
         potential_list_[site] = val*m;
      }
      else {
         potential_list_.at(site) = potential_list_.at(site) + val*m;
      }
      index_list_.emplace(site);
   }
   
   //! @brief Add an interaction term \f$ c\hat{O}_{i}\hat{P}_{j} \f$ to the system.
   //! @param site_1 The site index \f$ i\f$.
   //! @param m_1 The onsite operator \f$ \hat{O}_{i} \f$.
   //! @param site_2 The site index \f$ j\f$.
   //! @param m_2 The onsite operator \f$ \hat{P}_{j} \f$.
   //! @param val The coefficient of the interaction term \f$ c \f$.
   void AddInteraction(const IndexType &site_1, const CRS &m_1, const IndexType &site_2, const CRS &m_2, const RealType val = 1.0) {
      if (val == 0.0) {
         return;
      }
      if (m_1.row_dim != m_1.col_dim || m_2.row_dim != m_2.col_dim) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __FUNCTION__ << " in "<< __FILE__ << std::endl;
         ss << "The input matrix is not a square matrix." << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (m_1.row_dim != this->GetDimOnsite()) {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __FUNCTION__ << " in "<< __FILE__ << std::endl;
         ss << "The input matrix is invalid" << std::endl;
         throw std::runtime_error(ss.str());
      }

      if (site_1 == site_2) {
         AddPotential(site_1, m_1*m_2, val);
         return;
      }
      interaction_list_[site_1][site_2].push_back({val*m_1, m_2});
      index_list_.emplace(site_1);
      index_list_.emplace(site_2);
   }
   
   //------------------------------------------------------------------
   //----------------------Access Member variables---------------------
   //------------------------------------------------------------------
   //! @brief Get the system size.
   //! @return The system size.
   int GetSystemSize() const {
      return static_cast<int>(index_list_.size());
   }
   
   //! @brief Get the site index list.
   //! @return The site index list.
   const std::unordered_set<IndexType, IndexHash> &GetIndexSet() const {
      return index_list_;
   }
   
   //! @brief Get onsite potential term list added to the system.
   //! @return Onsite potential term list added to the system.
   inline const OnsiteListType &GetPotentialList() const { return potential_list_; }
   
   //! @brief Get an onsite potential term \f$ c\hat{O}_{i}\f$ added to the system.
   //! @param site The site index \f$ i\f$.
   //! @return The onsite potential term \f$ c\hat{O}_{i}\f$ added to the system.
   inline const CRS &GetPotential(const IndexType &site) const { return potential_list_.at(site); }
   
   //! @brief Get interaction term list added to the system.
   //! @return Interaction term list added to the system.
   inline const IntersiteListType &GetInteractionList() const { return interaction_list_; }
   
   //! @brief Get an interaction term \f$ c\hat{O}_{i}\hat{P}_{j}\f$ added to the system.
   //! @param site_1 The site index \f$ i\f$.
   //! @param site_2 The site index \f$ j\f$.
   //! @return The interaction term added to the system \f$ c\hat{O}_{i}\f$ and \f$ \hat{P}_{j}\f$.
   inline const std::vector<std::pair<CRS, CRS>> &GetInteraction(const IndexType &site_1, const IndexType &site_2) const {
      return interaction_list_.at(site_1).at(site_2);
   }
   
private:
   //------------------------------------------------------------------
   //---------------------Private Member Variables---------------------
   //------------------------------------------------------------------
   //! @brief Onsite potential term list added to the system.
   OnsiteListType potential_list_;
   
   //! @brief Interaction term list added to the system.
   IntersiteListType interaction_list_;
   
   //! @brief The site index list.
   std::unordered_set<IndexType, IndexHash> index_list_;
   
};

} //namespace model
} //namespace compnal





#endif /* COMPNAL_MODEL_GENERAL_MODEL_HPP_ */
