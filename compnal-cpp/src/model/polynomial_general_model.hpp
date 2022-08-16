//
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
//  polynomial_general_model.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_MODEL_POLYNOMIAL_GENERAL_MODEL_HPP_
#define COMPNAL_MODEL_POLYNOMIAL_GENERAL_MODEL_HPP_

#include "../utility/hash.hpp"
#include "../utility/type.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace compnal {
namespace model {

template<typename RealType>
class PolynomialGeneralModel {
  
public:
   using IndexType = utility::IndexType;
   using IndexHash = utility::IndexHash;
   using VectorHash = utility::IndexVectorHash;
   using InteractionType = std::unordered_map<std::vector<IndexType>, RealType, VectorHash>;
   
   PolynomialGeneralModel() {}
   
   void AddInteraction(std::vector<IndexType> index_list, const RealType value) {
      std::sort(index_list.begin(), index_list.end());
      interaction_[index_list] += value;
      index_set_.insert(index_list.begin(), index_list.end());
      flag_need_reset_index_map_ = true;
   }
   
   std::vector<IndexType> GenerateIndexList() const {
      std::vector<IndexType> index_list(index_set_.begin(), index_set_.end());
      std::sort(index_list.begin(), index_list.end());
      return index_list;
   }

   const InteractionType &GetInteraction() const {
      return interaction_;
   }
   
   const std::unordered_set<IndexType, IndexHash> &GetIndexSet() const {
      return index_set_;
   }
   
   std::size_t GetSystemSize() const {
      return index_set_.size();
   }
   
   const std::unordered_map<IndexType, std::int64_t> GetIndexMap() const {
      if (flag_need_reset_index_map_) {
         std::vector<IndexType> index_list = GenerateIndexList();
         std::int64_t count = 0;
         index_map_.clear();
         for (std::size_t i = 0; i < index_list.size(); ++i) {
            index_map_[index_list[i]] = count;
            count++;
         }
         flag_need_reset_index_map_ = false;
      }
      return index_map_;
   }

private:
   std::unordered_set<IndexType, IndexHash> index_set_;
   std::unordered_map<IndexType, std::int64_t, IndexHash> index_map_;
   InteractionType interaction_;
   bool flag_need_reset_index_map_ = true;
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_GENERAL_MODEL_HPP_ */
