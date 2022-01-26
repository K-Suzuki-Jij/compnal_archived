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
//  Created by Kohei Suzuki on 2022/01/26.
//

#ifndef COMPNAL_TEST_UTILITY_TYPE_ALIAS_HPP_
#define COMPNAL_TEST_UTILITY_TYPE_ALIAS_HPP_

#include "../../src/model/all.hpp"

namespace compnal {
namespace test {

using compnal::test::ExpectEQ;
using compnal::test::ExpectNear;
using compnal::sparse_matrix::CRS;
using compnal::sparse_matrix::CRSTag;
using compnal::model::BaseU1Electron_1D;
using compnal::model::BaseU1Spin_1D;
using compnal::model::BaseU1SpinElectron_1D;
using compnal::model::BaseU1SpinMultiElectrons_1D;
using compnal::sparse_matrix::CRS;

}
}


#endif /* COMPNAL_TEST_UTILITY_TYPE_ALIAS_HPP_ */
