//
//  u1_electron_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/27.
//

#ifndef COMPNAL_MODEL_U1_ELECTRON_1D_HPP_
#define COMPNAL_MODEL_U1_ELECTRON_1D_HPP_

#include "../sparse_matrix/all.hpp"

#include <unordered_map>
#include <vector>

namespace compnal {
namespace model {

template<class BaseClass>
class GeneralModel_1D: public BaseClass {
   
   using RealType = typename BaseClass::ValueType;
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   using BaseClass::BaseClass;
   
   void AddOnsitePotential(const RealType val, const CRS &m, const int site) {
      if (val == 0.0) {
         return;
      }
      if (site < 0 || this->system_size_ <= site) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
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
   
   void AddOnsitePotential(const CRS &m, const int site) {
      AddOnsitePotential(1.0, m, site);
   }
   
   void AddInteraction(const RealType val, const CRS &m_1, const int site_1, const CRS &m_2, const int site_2) {
      if (val == 0.0) {
         return;
      }
      if (site_1 < 0 || this->system_size_ <= site_1) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid site: " << site_1 << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site_2 < 0 || this->system_size_ <= site_2) {
         std::stringstream ss;
         ss << "Error in " << __FUNCTION__ << std::endl;
         ss << "Invalid site: " << site_2 << std::endl;
         throw std::runtime_error(ss.str());
      }
      if (site_1 == site_2) {
         AddOnsitePotential(val, sparse_matrix::CalculateMatrixMatrixProduct(val, m_1, 1.0, m_2), site_1);
         return;
      }
      
      intersite_operator_list_[site_1][site_2].push_back({m_1.MultiplyByScalar(val), m_2});
   }
   
   void AddInteraction(const CRS &m_1, const int site_1, const CRS &m_2, const int site_2) {
      AddInteraction(1.0, m_1, site_1, m_2, site_2);
   }
   
   inline const std::unordered_map<int, CRS> &GetOnsiteOperatorList() const { return onsite_operator_list_; }
   inline const CRS &GetOnsiteOperator(const int site) const { return onsite_operator_list_.at(site); }
   inline const std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<CRS, CRS>>>> &GetIntersiteOperatorList() const { return intersite_operator_list_; }
   inline const std::vector<std::pair<CRS, CRS>> &GetIntersiteOperator(const int site_1, const int site_2) const { return intersite_operator_list_.at(site_1).at(site_2); }
   
private:
   std::unordered_map<int, CRS> onsite_operator_list_;
   std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<CRS, CRS>>>> intersite_operator_list_;
   
};

}
}





#endif /* COMPNAL_MODEL_U1_ELECTRON_1D_HPP_ */
