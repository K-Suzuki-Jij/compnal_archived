//
//  model.hpp
//  compnal
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

int DoubleTheNumber(double s) {
   CheckHalfInteger(s);
   return static_cast<int>(2*s);
};

}
}


#endif /* COMPNAL_UTILITY_MODEL_HPP_ */
