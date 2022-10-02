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
//  quadratic_general_model.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_MODEL_QUADRATIC_GENERAL_MODEL_HPP_
#define COMPNAL_MODEL_QUADRATIC_GENERAL_MODEL_HPP_

#include "../utility/hash.hpp"
#include "../utility/type.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace compnal {
namespace model {

template<typename RealType>
class QuadraticGeneralModel {
  
public:
   using OPType = utility::SpinType;
   using IndexType = utility::AnyIndexType;
   using IndexHash = utility::AnyIndexHash;
   using PairHash  = utility::AnyIndexPairHash;
   using LinearType = std::unordered_map<IndexType, RealType, IndexHash>;
   using QuadraticType = std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash>;
      
   QuadraticGeneralModel(const LinearType &linear,
                         const QuadraticType &quadratic) {
      
      std::unordered_set<IndexType, IndexHash> index_set;
      for (const auto &it: linear) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            index_set.emplace(it.first);
         }
      }
      for (const auto &it: quadratic) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            index_set.emplace(it.first.first);
            index_set.emplace(it.first.second);
         }
      }
      
      index_list_ = std::vector<IndexType>(index_set.begin(), index_set.end());
      std::sort(index_list_.begin(), index_list_.end());
      
      std::int32_t count = 0;
      for (std::size_t i = 0; i < index_list_.size(); ++i) {
         index_map_[index_list_[i]] = count;
         count++;
      }
      
      linear_.resize(index_list_.size());
      for (const auto &it: linear) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            linear_[index_map_.at(it.first)] = it.second;
            degree_ = 1;
         }
      }
      
      std::unordered_map<std::pair<std::int32_t, std::int32_t>, RealType, PairHash> new_quadratic;
      for (const auto &it: quadratic) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            const std::int32_t key_1 = index_map_.at(it.first.first);
            const std::int32_t key_2 = index_map_.at(it.first.second);
            if (key_1 < key_2) {
               new_quadratic[{key_1, key_2}] += it.second;
               degree_ = 2;
            }
            else if (key_1 > key_2) {
               new_quadratic[{key_2, key_1}] += it.second;
               degree_ = 2;
            }
            else {
               constant_ += it.second;
            }
         }
      }
      
      std::vector<std::vector<std::pair<std::int32_t, RealType>>> adjacency_list(index_list_.size());
      for (const auto &it: new_quadratic) {
         adjacency_list[it.first.first].push_back({it.first.second, it.second});
         adjacency_list[it.first.second].push_back({it.first.first, it.second});
      }
      
      for (std::size_t i = 0; i < adjacency_list.size(); ++i) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end(), [](const auto &a, const auto &b) {
            return a.first < b.first;
         });
      }
      
      row_ptr_.resize(adjacency_list.size() + 1);
      row_ptr_[0] = 0;
      for (std::size_t i = 0; i < adjacency_list.size(); ++i) {
         row_ptr_[i + 1] = static_cast<std::int64_t>(adjacency_list[i].size());
         for (const auto it: adjacency_list[i]) {
            col_ptr_.push_back(it.first);
            val_ptr_.push_back(it.second);
         }
      }
      
      for (std::size_t i = 0; i < row_ptr_.size(); ++i) {
         row_ptr_[i + 1] += row_ptr_[i];
      }
      
   }
   
   const std::vector<IndexType> &GetIndexList() const {
      return index_list_;
   }
   
   RealType GetConstant() const {
      return constant_;
   }
      
   const std::vector<RealType> &GetLinear() const {
      return linear_;
   }
   
   std::int32_t GetSystemSize() const {
      return static_cast<std::int32_t>(index_list_.size());
   }

   const std::unordered_map<IndexType, std::int32_t, IndexHash> &GetIndexMap() const {
      return index_map_;
   }
   
   std::int32_t GetDegree() const {
      return degree_;
   }
   
   const std::vector<std::int64_t> &GetRowPtr() const {
      return row_ptr_;
   }
   
   const std::vector<std::int32_t> &GetColPtr() const {
      return col_ptr_;
   }
   
   const std::vector<RealType> &GetValPtr() const {
      return val_ptr_;
   }
      
private:
   int32_t degree_ = 0;
   std::unordered_map<IndexType, std::int32_t, IndexHash> index_map_;
   RealType constant_ = 0;
   std::vector<RealType> linear_;
   std::vector<IndexType> index_list_;
   
   std::vector<std::int64_t> row_ptr_;
   std::vector<std::int32_t> col_ptr_;
   std::vector<RealType> val_ptr_;
   
};


} // namespace model
} // namespace compnal

#endif /* COMPNAL_MODEL_QUADRATIC_GENERAL_MODEL_HPP_ */
