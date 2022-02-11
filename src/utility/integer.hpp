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
//  Created by Kohei Suzuki on 2021/11/07.
//

#ifndef COMPNAL_UTILITY_INTEGER_HPP_
#define COMPNAL_UTILITY_INTEGER_HPP_

#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>

namespace compnal {
namespace utility {


template<typename IntegerType>
std::vector<std::vector<IntegerType>> GenerateIntegerPartition(const IntegerType partitioned_number, IntegerType max_number) {
   static_assert(std::is_integral<IntegerType>::value, "Template parameter IntegerType must be integer type");
   
   if (partitioned_number <= 0 || max_number <= 0) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "Invalid input parameters" << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   auto generate_next = [](const std::vector<IntegerType> &vec) -> std::vector<IntegerType> {
      std::int64_t size = static_cast<std::int64_t>(vec.size());
      std::vector<IntegerType> out = vec;
      for (std::int64_t i = size - 1; i >= 0; --i) {
         if (vec[i] > 1) {
            if (i + 1 < size) {
               if (vec[i] - 1 >= vec[i + 1] + 1) {
                  out[i]--;
                  out[i + 1]++;
                  return out;
               }
               else {
                  out[i]--;
                  out.push_back(1);
                  return out;
               }
            }
            else {
               out[i]--;
               out.push_back(1);
               return out;
            }
         }
      }
      return std::vector<IntegerType>();
   };
   
   std::vector<std::vector<IntegerType>> out;
   
   if (max_number >= partitioned_number) {
      max_number = partitioned_number;
      out.push_back({max_number});
   }
   else {
      std::vector<IntegerType> temp = {max_number};
      IntegerType rem = partitioned_number - max_number;
      while (rem != 0) {
         if (rem <= max_number) {
            temp.push_back(rem);
            rem = 0;
         }
         else {
            temp.push_back(max_number);
            rem -= max_number;
         }
      }
      out.push_back(temp);
   }
   
   while (true) {
      auto next_partition = generate_next(out.back());
      if (next_partition.size() > 0) {
         out.push_back(next_partition);
      }
      else {
         break;
      }
   }
   
   return out;
   
}

std::int64_t CalculateBinomialCoefficient(int n, const int k) {
   if (n < 0 || k < 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   std::int64_t r = 1;
   
   for (int d = 1; d <= k; d++) {
       r *= n--;
       r /= d;
   }
   
   return r;
   
}

std::vector<std::vector<std::int64_t>> CalculateBinomialTable(int N) {
   std::vector<std::vector<std::int64_t>> vec(N + 1, std::vector<std::int64_t>(N + 1));
   for (int i = 0; i <= N; ++i) {
      for (int j = 0; j <= i; j++) {
         if (j == 0 || j == i) {
            vec[i][j] = 1;
         }
         else {
            vec[i][j] = vec[i - 1][j - 1] + vec[i - 1][j];
         }
      }
   }
   return vec;
}

std::int64_t CalculateNumCombination(const std::vector<int> &list) {
   
   std::unordered_map<int, int> u_map;
   
   for (const auto &it: list) {
      u_map[it]++;
   }
   
   std::int64_t result = 1;
   int size_list = static_cast<int>(list.size());
   
   for (const auto &it: u_map) {
      result *= CalculateBinomialCoefficient(size_list, it.second);
      size_list -= it.second;
   }
   
   return result;
   
}

void CalculateNthPermutation(std::vector<int> *vec, std::int64_t target_num) {
   
   std::map<int, int> u_map;
   
   for (const auto &it: *vec) {
      u_map[it]++;
   }
   
   int size_vec = static_cast<int>(vec->size());
   
   for (int i = 0; i < size_vec; ++i) {
      std::int64_t temp1 = 0;
      for (auto &&it: u_map) {
         if (it.second > 0) {
            it.second -= 1;
            std::int64_t temp2 = 1;
            int size = size_vec - (i + 1);
            for (const auto it2: u_map) {
               temp2 *= CalculateBinomialCoefficient(size, it2.second);
               size  -= it2.second;
            }
            temp1 += temp2;
            if (temp1 > target_num) {
               (*vec)[i] = it.first;
               temp1 -= temp2;
               break;
            }
            it.second += 1;
         }
      }
      target_num -= temp1;
   }
   
   if (target_num != 0) {
      std::stringstream ss;
      ss << "Error in " << __FUNCTION__ << " at " << __LINE__ << std::endl;
      ss << "Can't find corresponding permutation, " << target_num << std::endl;
      throw std::runtime_error(ss.str());
   }
   
}

}
}

#endif /* COMPNAL_UTILITY_INTEGER_HPP_ */
