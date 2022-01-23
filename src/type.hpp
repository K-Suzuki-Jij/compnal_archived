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
//  Created by Kohei Suzuki on 2022/01/23.
//

#ifndef COMPNAL_TYPE_HPP_
#define COMPNAL_TYPE_HPP_

#include <cstdint>
#include <utility>
#include <unordered_map>

namespace compnal {

//! @brief Half integer type.
using HalfInt = double;

//! @brief Long integer type.
using LInt = std::int64_t;

//! @brief Integer pair.
using IntPair = std::pair<int, int>;

//! @brief Long integer pair.
using LIntPair = std::pair<LInt, LInt>;

//! @brief std::vector<int> and int pair.
using VecIntPair = std::pair<std::vector<int>, int>;

//! @brief Alias of std::unordered_map.
//! @tparam Key Type of key.
//! @tparam T Type of value.
//! @tparam Hash Type of hash function.
template<class Key, class T, class Hash = std::hash<Key>>
using Map = std::unordered_map<Key, T, Hash>;

//! @brief Alias of hash class for std::pair.
struct PairHash {
   template<class T1, class T2>
   std::size_t operator() (const std::pair<T1, T2>& p) const {
      std::size_t lhs = std::hash<T1>()(p.first), rhs = std::hash<T2>()(p.second);
      return lhs^(rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2));
   }
};

//! @brief Alias of hash class for std::vector.
struct VecHash {
   template<class T>
   std::size_t operator() (const std::vector<T> &V) const {
      std::size_t hash = V.size();
      for (auto &i : V) {
         hash ^= std::hash<T>()(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      return hash;
   }
};

//! @brief Alias of hash class for std::pair including std::vector.
struct VecIntHash {
   template<class T1, class T2>
   std::size_t operator() (const std::pair<std::vector<T1>, T2>& p) const {
      std::size_t hash = p.first.size();
      for (auto &i: p.first) {
         hash ^= std::hash<T1>()(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      std::size_t lhs = std::hash<T2>()(p.second);
      hash ^= lhs^(0x9e3779b9 + (lhs << 6));
      return hash;
   }
};


} //namespace compnal


#endif /* COMPNAL_TYPE_HPP_ */
