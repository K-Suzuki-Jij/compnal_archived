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
//  type.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_UTILITY_TYPE_HPP_
#define COMPNAL_UTILITY_TYPE_HPP_

namespace compnal {
namespace utility {

using IntStrType = std::variant<std::int32_t, std::string>;

using IndexType = std::variant<std::int32_t, std::string, std::vector<IntStrType>>;

using SpinType = std::int8_t;


} // namespace utility
} // namespace compnal

#endif /* COMPNAL_UTILITY_TYPE_HPP_ */
