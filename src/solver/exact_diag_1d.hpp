//
//  exact_diag_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#ifndef exact_diag_1d_hpp
#define exact_diag_1d_hpp

#include "sparse_matrix.hpp"
#include "model.hpp"
#include "exact_diag_utility.hpp"
#include <unordered_map>

namespace compnal {
namespace solver {

template<typename ModelClass1D>
class ExactDiag1D {
   
   using RealType = typename ModelClass1D::ValueType;
   
   using CRS = sparse_matrix::CRS<RealType>;

   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   
   const ModelClass1D model;
   
   explicit ExactDiag1D(const ModelClass1D &model_input): model(model_input) {}
   
   const BraketVector &GetGSVector() const {
      return gs_vector;
   }
   
   RealType CalculateExpectationValue(const CRS &M, int site) {
      return gs_vector*CalculateMatrixVectorProduct<RealType>(M, site, gs_vector, bases[model.GetTotal2Sz()], bases[model.GetTotal2Sz()]);
   }
   
   
private:
   std::unordered_map<int, std::vector<int64_t>> bases;
   BraketVector gs_vector;
   
};

} // namespace solver
} // namespace compnal



#endif /* exact_diag_1d_hpp */
