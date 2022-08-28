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

namespace compnal {
namespace model {

template<typename RealType>
class QuadraticGeneralModel {
  
public:
   using IndexType = utility::AnyIndexType;
   using IndexHash = utility::AnyIndexHash;
   using PairHash  = utility::AnyIndexPairHash;
   
   QuadraticGeneralModel(const std::unordered_map<IndexType, RealType, IndexHash> &linear,
                         const std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash> &quadratic) {
      for (const auto &it: linear) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            index_set_.emplace(it.first);
            linear_[it.first] += it.second;
         }
      }
      for (const auto &it: quadratic) {
         if (std::abs(it.second) > std::numeric_limits<RealType>::epsilon()) {
            index_set_.emplace(it.first.first);
            index_set_.emplace(it.first.second);
            if (it.first.first < it.first.second) {
               quadratic_[{it.first.first, it.first.second}] += it.second;
            }
            else {
               quadratic_[{it.first.second, it.first.first}] += it.second;
            }
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

   
private:
   int32_t degree_ = 0;
   std::unordered_set<IndexType, IndexHash> index_set_;
   std::unordered_map<IndexType, std::int64_t, IndexHash> index_map_;
   RealType constant_ = 0;
   std::unordered_map<IndexType, RealType, IndexHash> linear_;
   std::unordered_map<std::pair<IndexType, IndexType>, RealType, PairHash> quadratic_;

};


} // namespace model
} // namespace compnal

#endif /* COMPNAL_MODEL_QUADRATIC_GENERAL_MODEL_HPP_ */
