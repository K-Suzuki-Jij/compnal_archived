//
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
//  random_number.hpp
//  compnal
//
//  Created by kohei on 2022/10/02.
//  
//

#ifndef COMPNAL_UTILITY_RANDOM_NUMBER_HPP_
#define COMPNAL_UTILITY_RANDOM_NUMBER_HPP_

#include <random>

namespace compnal {
namespace utility {

class Xorshift {
public:
   using result_type = std::uint_fast32_t;
   
   static constexpr std::uint32_t min() { return 0u; }
   static constexpr std::uint32_t max() { return UINT_MAX; }
   std::uint32_t operator()() {
      std::uint32_t t = x_ ^ (x_ << 11);
      x_ = y_;
      y_ = z_;
      z_ = w_;
      return w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));
   }
   
   Xorshift() {
      std::random_device rd;
      w_ = rd();
   }
   
   Xorshift(std::uint32_t s) { w_ = s; }
   
private:
   std::uint32_t x_ = 123456789u;
   std::uint32_t y_ = 362436069u;
   std::uint32_t z_ = 521288629u;
   std::uint32_t w_;
};


} // namespace utility
} // namespace compnal

#endif /* COMPNAL_UTILITY_RANDOM_NUMBER_HPP_ */
