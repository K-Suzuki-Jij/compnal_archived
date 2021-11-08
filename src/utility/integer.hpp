//
//  integer.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/07.
//

#ifndef COMPNAL_UTILITY_INTEGER_HPP_
#define COMPNAL_UTILITY_INTEGER_HPP_

#include <sstream>
#include <vector>
#include <unordered_map>

namespace compnal {
namespace utility {

void GenerateIntegerPartition(std::vector<std::vector<int>> *results, int n, int max_partition_num) {
   
   std::vector<int> temp_vec;
   
   if (results->size() > 0 && n >= max_partition_num) {
      temp_vec = results->back();
   }
   
   if (results->size() == 0 && n == 0) {
      results->push_back(std::vector<int>());
   }
   else if (n == 1) {
      if (results->size() == 0) {
         results->push_back(std::vector<int>());
      }
      results->back().push_back(1);
   }
   else if (max_partition_num == 1) {
      if (results->size() == 0) {
         results->push_back(std::vector<int>());
      }
      for (int i = 0; i < n; ++i) {
         results->back().push_back(1);
      }
   }
   else {
      if (n >= max_partition_num) {
         if (results->size() == 0) {
            results->push_back(std::vector<int>());
         }
         results->back().push_back(max_partition_num);
         GenerateIntegerPartition(results, n - max_partition_num, max_partition_num);
         results->push_back(std::vector<int>());
         for (const auto &v: temp_vec) {
            results->back().push_back(v);
         }
      }
      GenerateIntegerPartition(results, n, max_partition_num - 1);
   }
}


std::size_t CalculateBinomialCoefficient(int n, const int k) {
   if (n < 0 || k < 0) {
      std::stringstream ss;
      ss << "Error in " << __func__ << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   std::size_t r = 1;
   
   for (int d = 1; d <= k; d++) {
       r *= n--;
       r /= d;
   }
   
   return r;
   
}

std::size_t CalculateNumCombination(const std::vector<int> &list) {
   
   std::unordered_map<int, int> u_map;
   
   for (const auto &it: list) {
      u_map[it]++;
   }
   
   std::size_t result = 1;
   int size_list = static_cast<int>(list.size());
   
   for (const auto &it: u_map) {
      result *= CalculateBinomialCoefficient(size_list, it.second);
      size_list -= it.second;
   }
   
   return result;
   
}

void CalculateNthPermutation(std::vector<int> *vec, std::size_t target_num) {
   
   std::unordered_map<int, int> u_map;
   
   for (const auto &it: *vec) {
      u_map[it]++;
   }
   
   int size_vec = static_cast<int>(vec->size());
   
   for (int i = 0; i < size_vec; ++i) {
      std::size_t temp1 = 0;
      for (auto &&it: u_map) {
         if (it.second > 0) {
            it.second -= 1;
            std::size_t temp2 = 1;
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
      ss << "Error in " << __FUNCTION__ << std::endl;
      ss << "Can't find corresponding permutation, " << target_num << std::endl;
      throw std::runtime_error(ss.str());
   }
   
}


}
}

#endif /* COMPNAL_UTILITY_INTEGER_HPP_ */
