//
//  hash.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/02.
//

#ifndef COMPNAL_UTILITY_HASH_HPP_
#define COMPNAL_UTILITY_HASH_HPP_

namespace compnal {
namespace utility {

struct PairHash {
   template<class T1, class T2>
   std::size_t operator() (const std::pair<T1, T2>& p) const {
      std::size_t lhs = std::hash<T1>()(p.first), rhs = std::hash<T2>()(p.second);
      return lhs^(rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2));
   }
};

struct VectorHash {
   template<class T>
   std::size_t operator() (const std::vector<T> &V) const {
      std::size_t hash = V.size();
      for (auto &i : V) {
         hash ^= std::hash<T>()(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      return hash;
   }
};

struct VectorIntHash {
   template<class T1, class T2>
   std::size_t operator() (const std::pair<T1, std::vector<T2>>& p) const {
      std::size_t hash = p.second.size();
      for (auto &i: p.second) {
         hash ^= std::hash<T2>()(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      std::size_t lhs = std::hash<T1>()(p.first);
      hash ^= lhs^(0x9e3779b9 + (lhs << 6));
      return hash;
   }
};


}
}

#endif /* COMPNAL_UTILITY_HASH_HPP_ */
