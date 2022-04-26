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
//  Created by Kohei Suzuki on 2022/04/27.
//

#ifndef COMPNAL_BLAS_TDMA_HPP_
#define COMPNAL_BLAS_TDMA_HPP_

namespace compnal {
namespace blas {

template<typename RealType>
void Tdma(std::vector<RealType> *vector_out,
          const std::vector<RealType> &diag,
          const std::vector<RealType> &off_diag,
          const std::vector<RealType> &vector_in,
          const RealType diag_add = 0.0) {
   
   if (off_diag.size() + 1 != diag.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "diag size=" << diag.size() << ", off_diag size=" << off_diag.size() << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   if (vector_in.size() != diag.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in "<< __FILE__ << std::endl;
      ss << "Matrix vector product (Ax=b) cannot be defined." << std::endl;
      throw std::runtime_error(ss.str());
   }
   
   const std::int64_t dim = static_cast<std::int64_t>(diag.size());
   
   if (dim == 1) {
      vector_out->resize(1);
      (*vector_out)[0] = diag[0]/vector_in[0];
      return;
   }
   
   std::vector<RealType> ppp(dim - 1);
   std::vector<RealType> qqq(dim);
   
   ppp[0] = -off_diag[0]/(diag[0] + diag_add);
   qqq[0] = vector_in[0]/(diag[0] + diag_add);
   
   for (std::int64_t i = 1; i < dim - 1;) {
      ppp[i] = -off_diag[i]/(diag[i]+ diag_add + off_diag[i - 1]*ppp[i - 1]);
      qqq[i] = (vector_in[i] - off_diag[i - 1]*qqq[i - 1])/(diag[i]+ diag_add + off_diag[i - 1]*ppp[i - 1]);
   }
   qqq[dim - 1] = (vector_in[dim - 1] - off_diag[dim - 2]*qqq[dim - 2])/(diag[dim - 1]+ diag_add + off_diag[dim - 2]*ppp[dim - 2]);
   
   vector_out->resize(dim);
   
   (*vector_out)[dim - 1] = qqq[dim - 1];
   
   for (std::int64_t i = dim - 2; i >= 0; --i) {
      (*vector_out)[i] = ppp[i]*(*vector_out)[i + 1] + qqq[i];
   }
   
}


} // namespace blas
} // namespace compnel

#endif /* COMPNAL_BLAS_TDMA_HPP_ */
