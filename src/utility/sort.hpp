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
//  Created by Kohei Suzuki on 2021/06/17.
//

#ifndef COMPNAL_UTILITY_SORT_HPP_
#define COMPNAL_UTILITY_SORT_HPP_

#include <vector>

namespace compnal {
namespace utility {

template<typename T1, typename T2>
void QuickSort(std::vector<T1> *array_comp, std::vector<T2> *array2, const std::int64_t left, const std::int64_t right) {
   
   if (right - left <= 1) {
      return;
   }
   
   const std::int64_t index_pivot = (left + right)/2;
   const T1          val_pivot   = (*array_comp)[index_pivot];
   
   std::swap((*array_comp)[index_pivot], (*array_comp)[right - 1]);
   std::swap((*array2)    [index_pivot], (*array2)    [right - 1]);
   
   std::int64_t index = left;
   for (std::int64_t i = left; i < right - 1; ++i) {
      if ((*array_comp)[i] < val_pivot) {
         std::swap((*array_comp)[index], (*array_comp)[i]);
         std::swap((*array2)    [index], (*array2)    [i]);
         ++index;
      }
   }

   std::swap((*array_comp)[index], (*array_comp)[right - 1]);
   std::swap((*array2)    [index], (*array2)    [right - 1]);
   
   QuickSort(array_comp, array2, left     , index);
   QuickSort(array_comp, array2, index + 1, right);
   
}


} // namespace utility
} // namespace compnel


#endif /* COMPNAL_UTILITY_SORT_HPP_ */
