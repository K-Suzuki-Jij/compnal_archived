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
   using OPType = utility::SpinType;
   using IndexType = utility::AnyIndexType;
   using IndexHash = utility::AnyIndexHash;
   using VectorHash = utility::AnyIndexVectorHash;
   using InteractionType = std::unordered_map<std::vector<IndexType>, RealType, VectorHash>;
   
   PolynomialGeneralModel(const InteractionType &interaction) {
      std::unordered_set<IndexType, IndexHash> index_set;
      for (const auto &it: interaction) {
         index_set.insert(it.first.begin(), it.first.end());
      }
      index_list_.reserve(index_set.size());
      for (const auto &index: index_set) {
         index_list_.push_back(index);
      }
      std::sort(index_list_.begin(), index_list_.end());
      
      std::int32_t count = 0;
      for (std::size_t i = 0; i < index_list_.size(); ++i) {
         index_map_[index_list_[i]] = count;
         count++;
      }
      
      key_list_.reserve(interaction.size());
      value_list_.reserve(interaction.size());
      for (const auto &it: interaction) {
         std::vector<std::int32_t> keys(it.first.size());
         for (std::size_t i = 0; i < it.first.size(); ++i) {
            keys[i] = index_map_.at(it.first[i]);
         }
         std::sort(keys.begin(), keys.end());
         key_list_.push_back(keys);
         value_list_.push_back(it.second);
         if (degree_ < keys.size()) {
            degree_ = static_cast<std::int32_t>(keys.size());
         }
      }
      
      adjacency_list_.resize(index_list_.size());
      for (std::size_t i = 0; i < key_list_.size(); ++i) {
         for (const auto &index: key_list_[i]) {
            adjacency_list_[index].push_back(i);
         }
      }
      
      // Save memory
      for (std::size_t i = 0; i < index_list_.size(); ++i) {
         adjacency_list_[i].shrink_to_fit();
         std::sort(adjacency_list_[i].begin(), adjacency_list_[i].end());
      }
   }
   
   std::pair<std::vector<std::vector<IndexType>>, std::vector<RealType>> GenerateInteractionAsPair() const {
      std::vector<std::vector<IndexType>> key_list;
      std::vector<RealType> value_list;
      key_list.reserve(key_list_.size());
      value_list.reserve(key_list_.size());
      for (std::size_t i = 0; i < key_list_.size(); ++i) {
         std::vector<IndexType> keys(key_list_[i].size());
         for (std::size_t j = 0; j < key_list_[i].size(); ++j) {
            keys[j] = index_list_[key_list_[i][j]];
         }
         key_list.push_back(keys);
         value_list.push_back(value_list_[i]);
      }
      return std::pair<std::vector<std::vector<IndexType>>, std::vector<RealType>>{key_list, value_list};
   }
   
   std::int32_t GetSystemSize() const {
      return static_cast<std::int32_t>(index_list_.size());
   }
   
   const std::vector<std::vector<std::int32_t>> &GetKeyList() const {
      return key_list_;
   }
   
   const std::vector<RealType> &GetValueList() const {
      return value_list_;
   }
   
   std::int32_t GetDegree() const {
      return degree_;
   }
   
   const std::vector<IndexType> &GetIndexList() const {
      return index_list_;
   }
   
   const std::unordered_map<IndexType, std::int32_t, IndexHash> &GetIndexMap() const {
      return index_map_;
   }
   
   const std::vector<std::vector<std::size_t>> &GetAdjacencyList() const {
      return adjacency_list_;
   }
   
private:
   int32_t degree_ = 0;
   std::unordered_map<IndexType, std::int32_t, IndexHash> index_map_;
   std::vector<std::vector<std::int32_t>> key_list_;
   std::vector<RealType> value_list_;
   std::vector<IndexType> index_list_;
   std::vector<std::vector<std::size_t>> adjacency_list_;
};


} // namespace model
} // namespace compnal


#endif /* COMPNAL_MODEL_POLYNOMIAL_GENERAL_MODEL_HPP_ */
