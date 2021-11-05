//
//  model.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/11/06.
//

#ifndef COMPNAL_UTILITY_MODEL_HPP_
#define COMPNAL_UTILITY_MODEL_HPP_

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

}
}


#endif /* COMPNAL_UTILITY_MODEL_HPP_ */
