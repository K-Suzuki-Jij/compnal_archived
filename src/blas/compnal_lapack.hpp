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
//  Created by Kohei Suzuki on 2021/05/22.
//

#ifndef COMPNAL_BLAS_COMPNAL_LAPACK_HPP_
#define COMPNAL_BLAS_COMPNAL_LAPACK_HPP_

#include <cmath>
#include <limits>
#include <random>
#include <sstream>

#include "compressed_row_storage.hpp"
#include "tdma.hpp"

namespace compnal {
namespace blas {

template <typename RealType>
RealType GetAccuracy() {
   return 100 * std::numeric_limits<RealType>::epsilon();
}

template <typename RealType>
std::int64_t CountSignFlip(const RealType value, const std::vector<RealType> &diag,
                           const std::vector<RealType> &off_diag, const RealType threshold) {
   if (off_diag.size() + 1 != diag.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "diag size=" << diag.size() << ", off_diag size=" << off_diag.size() << std::endl;
      throw std::runtime_error(ss.str());
   }

   const std::int64_t dim = static_cast<std::int64_t>(diag.size());
   std::int64_t count = 0;
   RealType p0 = 1.0;
   RealType p1 = value - diag[0];
   RealType p2 = 0.0;

   for (std::int64_t i = 0; i < dim - 1; ++i) {
      p2 = (value - diag[i + 1]) * p1 - off_diag[i] * off_diag[i] * p0;

      const bool c1 = (p0 < 0 && p1 > 0) || (p0 > 0 && p1 < 0);
      const bool c2 = (std::abs(p1) < threshold) && ((p0 < 0 && p2 > 0) || (p0 > 0 && p2 < 0));

      if (c1 || c2) {
         count++;
      }

      p1 = p2 / p1;
      p0 = 1.0;
   }

   if ((p0 < 0 && p1 > 0) || (p0 > 0 && p1 < 0)) {
      count++;
   }

   return count;
}

template <typename RealType, typename IntegerType>
void Dstev(RealType *eigenvalue, std::vector<RealType> *eigenvector, const IntegerType n,
           const std::vector<RealType> &diag, const std::vector<RealType> &off_diag) {
   if (off_diag.size() + 1 != diag.size()) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "diag size=" << diag.size() << ", off_diag size=" << off_diag.size() << std::endl;
      throw std::runtime_error(ss.str());
   }

   if (n < 0 || n >= static_cast<IntegerType>(diag.size())) {
      std::stringstream ss;
      ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
      ss << "n=" << n << "is invalid" << std::endl;
      throw std::runtime_error(ss.str());
   }

   const int inverse_iteration_repeat = 3;
   const std::int64_t dim = static_cast<std::int64_t>(diag.size());
   const std::int64_t nn = static_cast<std::int64_t>(dim - 1 - n);

   if (dim == 1) {
      *eigenvalue = diag[0];
      eigenvector->resize(1);
      (*eigenvector)[0] = 1;
      return;
   }
   if (dim == 2) {
      const auto a = diag[0];
      const auto b = off_diag[0];
      const auto c = diag[1];
      const auto root = std::sqrt((a - c) * (a - c) + 4 * b * b);
      const auto e1 = (a + c - std::sqrt((a - c) * (a - c) + 4 * b * b)) * RealType{0.5};
      const auto e2 = (a + c + std::sqrt((a - c) * (a - c) + 4 * b * b)) * RealType{0.5};
      if (n == 0) {
         *eigenvalue = e1;
         eigenvector->resize(2);
         (*eigenvector)[0] = (a - c - root) / std::sqrt(4 * b * b - (a - c - root));
         (*eigenvector)[1] = 2 * b / std::sqrt(4 * b * b - (a - c - root));
      } else if (n == 1) {
         *eigenvalue = e2;
         eigenvector->resize(2);
         (*eigenvector)[0] = (a - c + root) / std::sqrt(4 * b * b + a - c - root);
         (*eigenvector)[1] = 2 * b / std::sqrt(4 * b * b + a - c - root);
      } else {
         std::stringstream ss;
         ss << "Error at " << __LINE__ << " in " << __func__ << " in " << __FILE__ << std::endl;
         ss << "n=" << n << " is invalid" << std::endl;
         throw std::runtime_error(ss.str());
      }
      return;
   }

   RealType x0 = std::numeric_limits<RealType>::max();
   RealType x1 = std::numeric_limits<RealType>::lowest();

   for (std::int64_t i = 0; i < dim - 2; ++i) {
      x0 = std::min(x0, diag[i + 1] - (std::abs(off_diag[i]) + std::abs(off_diag[i + 1])));
      x1 = std::max(x1, diag[i + 1] + (std::abs(off_diag[i]) + std::abs(off_diag[i + 1])));
   }

   const RealType th = GetAccuracy<RealType>();
   const RealType coeff = RealType{0.5};

   while (x1 - x0 > th) {
      const RealType pp = x1 - x0;
      const RealType piv = (x0 + x1) * coeff;
      if (CountSignFlip(piv, diag, off_diag, th) <= nn) {
         x1 = piv;
      } else {
         x0 = piv;
      }
      if (std::abs(pp - (x1 - x0)) < th) {
         break;
      }
   }

   *eigenvalue = (x1 + x0) * coeff;

   // Calculate eigenvectors
   eigenvector->resize(dim);

   std::uniform_real_distribution<RealType> uniform_rand(-1, 1);
   std::mt19937 random_number_engine;
   random_number_engine.seed(std::random_device()());
   for (std::int64_t i = 0; i < dim; ++i) {
      (*eigenvector)[i] = uniform_rand(random_number_engine);
   }

   auto normalize = [](std::vector<RealType> *vector, const std::int64_t dim) {
      RealType norm = 0.0;
      for (std::int64_t i = 0; i < dim; ++i) {
         norm += (*vector)[i] * (*vector)[i];
      }
      norm = std::sqrt(norm);
      for (std::int64_t i = 0; i < dim; ++i) {
         (*vector)[i] = (*vector)[i] / norm;
      }
   };

   for (int i = 0; i < inverse_iteration_repeat; ++i) {
      Tdma(eigenvector, diag, off_diag, *eigenvector, -(*eigenvalue) + th);
      normalize(eigenvector, dim);
   }
}

}  // namespace blas
}  // namespace compnal

#endif /* CCOMPNAL_BLAS_COMPNAL_LAPACK_HPP_ */
