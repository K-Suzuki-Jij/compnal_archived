//
//  model_utility.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef model_utility_hpp
#define model_utility_hpp

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <unordered_map>

namespace compnal {
namespace model {

enum ElectronState {
   
   VACUUM  = 0,
   UP      = 1,
   DOWN    = 2,
   UP_DOWN = 3,

   EVEN     = 1,
   ODD      = 2,
   EVEN_ODD = 3

};

enum BoundaryCondition {
  
   OBC = 0,
   PBC = 1,
   SSD = 2
   
};

int DoubleTheNumber(double s) {
   s = 2*s;
   if (std::floor(s) != s) {
      throw std::runtime_error("Invalid value of magnitude_spin or total_sz");
   }
   return static_cast<int>(s);
};

void CheckSumOverflowINT64(int64_t a, int64_t b) {
   if (a < 0 || b < 0) {
      std::stringstream ss;
      ss << "Error in " << __FUNCTION__ << std::endl;
      ss << "Invalid parameters" << std::endl;
      throw std::runtime_error(ss.str());
   }
   if (a >= INT64_MAX - b) {
      throw std::runtime_error("Overflow detected for sumation using int64_t");
   }
}

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

template<typename IntegerType>
IntegerType CalculateBinomialCoefficient(IntegerType n, IntegerType k) {
   
   IntegerType r = 1;
   
   for (IntegerType d = 1; d <= k; d++) {
       r *= n--;
       r /= d;
   }
   
   return r;
   
}

template<typename IntegerType>
int64_t CalculateNumCombination(const std::vector<IntegerType> &list) {
   
   std::unordered_map<IntegerType, int64_t> u_map;
   
   for (const auto &it: list) {
      u_map[it]++;
   }
   
   int64_t result = 1;
   int64_t size_list = list.size();
   
   for (const auto &it: u_map) {
      result *= CalculateBinomialCoefficient<int64_t>(size_list, it.second);
      size_list -= it.second;
   }
   
   return result;
   
}

int64_t CalculatePower(int base, int exponent) {
   int64_t a = 1;
   for (int i = 0; i < exponent; ++i) {
      a *= base;
   }
   return a;
}

template<typename IntegerType>
void NthPermutation(std::vector<IntegerType> *vec, int64_t target_num) {
   
   std::unordered_map<IntegerType, int64_t> u_map;
   
   for (const auto &it: *vec) {
      u_map[it]++;
   }
   
   int64_t size_vec = vec->size();
   
   for (int64_t i = 0; i < size_vec; ++i) {
      int64_t temp1 = 0;
      for (auto &&it: u_map) {
         if (it.second > 0) {
            it.second -= -1;
            int64_t temp2 = 1;
            int64_t size = size_vec - (i + 1);
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


} // namespace lattice
} // namespace compnal

#endif /* model_utility_hpp */
