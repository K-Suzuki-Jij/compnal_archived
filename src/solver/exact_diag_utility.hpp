//
//  exact_diag_utility.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef exact_diag_utility_hpp
#define exact_diag_utility_hpp

#include <vector>
#include <unordered_map>

namespace compnal {
namespace solver {

template<typename RealType>
struct ExactDiagMatrixElements {
   std::vector<RealType> val;
   std::vector<int64_t>  basis_affected;
   std::vector<int>      basis_onsite;
   std::vector<int64_t>  site_constant;
   std::unordered_map<int64_t, int64_t> inv_basis_affected;
   RealType zero_precision = 0.0;
};




}
}


#endif /* exact_diag_utility_hpp */
