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

struct pair_hash {
   template <class T1, class T2>
   std::size_t operator() (const std::pair<T1, T2>& p) const {
      std::size_t lhs = std::hash<T1>()(p.first), rhs = std::hash<T2>()(p.second);
      return lhs^(rhs+0x9e3779b9+(lhs<<6)+(lhs>>2));
   }
};

}
}

#endif /* COMPNAL_UTILITY_HASH_HPP_ */
