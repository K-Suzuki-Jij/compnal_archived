//
//  model_utility.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/20.
//

#ifndef model_utility_hpp
#define model_utility_hpp

namespace compnal {
namespace model {

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


} // namespace lattice
} // namespace compnal

#endif /* model_utility_hpp */
