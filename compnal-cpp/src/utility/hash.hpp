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
//  hash.hpp
//  compnal
//
//  Created by kohei on 2022/08/10.
//  
//

#ifndef COMPNAL_UTILITY_HASH_HPP_
#define COMPNAL_UTILITY_HASH_HPP_

#include <variant>

namespace compnal {
namespace utility {

//! @brief Hash struct of IndexType used in model::GeneralModel
//! @tparam IntegerType Integer type.
template <typename IntegerType>
struct VariantHash {
   static_assert(std::is_integral<IntegerType>::value, "Template parameter, IntegerType must be integer type");

   using VariantVecType = std::vector<std::variant<IntegerType, std::string>>;

   template <class... Types>
   std::size_t operator()(const std::variant<Types...> &v) const {
      if (std::holds_alternative<IntegerType>(v)) {
         return std::hash<IntegerType>()(std::get<IntegerType>(v));
      }
      else if (std::holds_alternative<std::string>(v)) {
         return std::hash<std::string>()(std::get<std::string>(v));
      }
      else if (std::holds_alternative<VariantVecType>(v)) {
         const auto &variant_vec = std::get<VariantVecType>(v);
         std::size_t hash = variant_vec.size();
         for (const auto &i : variant_vec) {
            if (std::holds_alternative<IntegerType>(i)) {
               hash ^= std::hash<IntegerType>()(std::get<IntegerType>(i)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            else if (std::holds_alternative<std::string>(i)) {
               hash ^= std::hash<std::string>()(std::get<std::string>(i)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            else {
               throw std::runtime_error("Invalid template parameters");
            }
         }
         return hash;
      }
      else {
         throw std::runtime_error("Invalid template parameters");
      }
   }
};






} // namespace utility
} // namespace compnal


#endif /* COMPNAL_UTILITY_HASH_HPP_ */
