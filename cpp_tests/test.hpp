//
//  test.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2022/01/08.
//

#ifndef TEST_TEST_HPP_
#define TEST_TEST_HPP_

#include "../src/sparse_matrix/all.hpp"
#include <unordered_map>

namespace compnal {
namespace test {

template<class MapType>
bool ExpectEQ(const MapType &lhs, const MapType &rhs) {
   if (lhs.size() != rhs.size()) {
      return false;
   }
   for (const auto &it: lhs) {
      if (rhs.count(it.first) == 1) {
         if (rhs.at(it.first) != it.second) {
            return false;
         }
      }
      else {
         return false;
      }
   }
   for (const auto &it: rhs) {
      if (lhs.count(it.first) == 1) {
         if (lhs.at(it.first) != it.second) {
            return false;
         }
      }
      else {
         return false;
      }
   }
   return true;
}

template<typename IndexType, typename ValueType>
bool ExpectNear(const std::unordered_map<IndexType, ValueType> &lhs, const std::unordered_map<IndexType, ValueType> &rhs, const double threshold) {
   if (lhs.size() != rhs.size()) {
      return false;
   }
   for (const auto &it: lhs) {
      if (rhs.count(it.first) == 1) {
         if (std::abs(rhs.at(it.first) - it.second) > threshold) {
            return false;
         }
      }
      else {
         return false;
      }
   }
   for (const auto &it: rhs) {
      if (lhs.count(it.first) == 1) {
         if (std::abs(lhs.at(it.first) - it.second) > threshold) {
            return false;
         }
      }
      else {
         return false;
      }
   }
   return true;
}

template<typename RealType>
bool ExpectNear(const compnal::sparse_matrix::CRS<RealType> &lhs, const compnal::sparse_matrix::CRS<RealType> &rhs, const RealType threshold) {
   if (lhs.row_dim != rhs.row_dim) {
      return false;
   }
   if (lhs.col_dim != rhs.col_dim) {
      return false;
   }
   if (lhs.row.size() != rhs.row.size()) {
      return false;
   }
   if (lhs.tag != rhs.tag) {
      return false;
   }
   for (std::size_t i = 0; i < lhs.row.size(); ++i) {
      if (lhs.row[i] != rhs.row[i]) {
         return false;
      }
   }
   if (lhs.col.size() != rhs.col.size()) {
      return false;
   }
   for (std::size_t i = 0; i < lhs.col.size(); ++i) {
      if (lhs.col[i] != rhs.col[i]) {
         return false;
      }
   }
   if (lhs.val.size() != rhs.val.size()) {
      return false;
   }
   for (std::size_t i = 0; i < lhs.val.size(); ++i) {
      if (std::abs(lhs.val[i] - rhs.val[i]) > threshold) {
         return false;
      }
   }
   return true;
}

}
} //Test



#endif /* TEST_TEST_HPP_ */
