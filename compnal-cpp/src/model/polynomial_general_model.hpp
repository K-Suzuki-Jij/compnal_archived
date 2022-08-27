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
   using IndexType = utility::AnyIndexType;
   using IndexHash = utility::AnyIndexHash;
   using VectorHash = utility::AnyIndexVectorHash;
   using InteractionType = std::unordered_map<std::vector<IndexType>, RealType, VectorHash>;
   
   PolynomialGeneralModel(const InteractionType &interaction) {
      for (const auto &it: interaction) {
         std::vector<IndexType> index_list = it.first;
         std::sort(index_list.begin(), index_list.end());
         index_set_.insert(index_list.begin(), index_list.end());
         interaction_[index_list] += it.second;
         if (degree_ < index_list.size()) {
            degree_ = static_cast<std::int32_t>(index_list.size());
         }
      }
      std::vector<IndexType> index_list = GenerateIndexList();
      std::int64_t count = 0;
      for (std::size_t i = 0; i < index_list.size(); ++i) {
         index_map_[index_list[i]] = count;
         count++;
      }
   }
      
   std::vector<IndexType> GenerateIndexList() const {
      std::vector<IndexType> index_list(index_set_.begin(), index_set_.end());
      std::sort(index_list.begin(), index_list.end());
      return index_list;
   }
   
   std::pair<std::vector<std::vector<IndexType>>, std::vector<RealType>> GenerateInteractionAsPair() const {
      std::vector<std::vector<IndexType>> key_list;
      std::vector<RealType> value_list;
      key_list.reserve(interaction_.size());
      value_list.reserve(interaction_.size());
      for (const auto &it: interaction_) {
         key_list.push_back(it.first);
         value_list.push_back(it.second);
      }
      return std::pair<std::vector<std::vector<IndexType>>, std::vector<RealType>>{key_list, value_list};
   }
   
   std::int32_t GetSystemSize() const {
      return static_cast<std::int32_t>(index_set_.size());
   }
   
   const std::unordered_set<IndexType, IndexHash> &GetIndexSet() const {
      return index_set_;
   }
   
   const std::unordered_map<IndexType, std::int64_t, IndexHash> &GetIndexMap() const {
      return index_map_;
   }
   
   const InteractionType &GetInteraction() const {
      return interaction_;
   }
   
   std::int32_t GetDegree() const {
      return degree_;
   }
   
private:
   int32_t degree_ = 0;
   std::unordered_set<IndexType, IndexHash> index_set_;
   std::unordered_map<IndexType, std::int64_t, IndexHash> index_map_;
   InteractionType interaction_;

};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_GENERAL_MODEL_HPP_ */
