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
//  Created by Kohei Suzuki on 2021/11/06.
//

#ifndef COMPNAL_UTILITY_MODEL_HPP_
#define COMPNAL_UTILITY_MODEL_HPP_

#include <cmath>

namespace compnal {
namespace utility {

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


void CheckHalfInteger(double s) {
   s = 2*s;
   if (std::floor(s) != s) {
      throw std::runtime_error("The input number is not half-integer");
   }
}

int DoubleHalfInteger(double s) {
   CheckHalfInteger(s);
   return static_cast<int>(2*s);
};

}
}


#endif /* COMPNAL_UTILITY_MODEL_HPP_ */
