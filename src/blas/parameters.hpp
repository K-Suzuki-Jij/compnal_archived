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
//  Created by Kohei Suzuki on 2022/03/05.
//

#ifndef COMPNAL_BLAS_PARAMETERS_HPP_
#define COMPNAL_BLAS_PARAMETERS_HPP_

namespace compnal {
namespace blas {

const int TIME_UNIT_CONSTANT = 1000*1000;

template<typename RealType>
struct ParametersLanczos {
   int min_step = 0;
   int max_step = 1000;
   RealType acc = std::pow(10, -14);
   bool flag_use_initial_vec = false;
   bool flag_store_vec       = false;
   bool flag_display_info    = true;
   bool flag_symmetric_crs   = false;
};

template<typename RealType>
struct ParametersCG {
   int max_step = 1000;
   RealType acc = std::pow(10, -7);
   bool flag_use_initial_vec = false;
   bool flag_display_info    = true;
   bool flag_symmetric_crs   = false;
};

template<typename RealType>
struct ParametersII {
   ParametersII() {
      cg.flag_use_initial_vec = true;
   }
   
   int      max_step = 3;
   RealType acc      = std::pow(10, -7);
   RealType diag_add = std::pow(10, -11);
   bool     flag_display_info = true;
   ParametersCG<RealType> cg;
};

template<typename RealType>
struct Parameters {
   ParametersLanczos<RealType> lanczos;
   ParametersII<RealType> ii;
};




} // namespace blas
} // namespace compnel


#endif /* COMPNAL_BLAS_PARAMETERS_HPP_ */
