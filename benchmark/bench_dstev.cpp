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
//  Created by Kohei Suzuki on 2022/04/27.
//

#include "../src/blas/all.hpp"
#include <random>

int main(void) {
   
   int max_dim = 1000;
   
   for (int dim = 10; dim < max_dim; dim += 100) {
      std::vector<double> diag(dim);
      std::vector<double> off_diag(dim - 1);
      
      std::uniform_real_distribution<double> uniform_rand(-1, 1);
      const unsigned int seed = std::random_device()();
      std::mt19937 random_number_engine;
      random_number_engine.seed(seed);
      
      for (int i = 0; i < dim; ++i) {
         diag[i] = uniform_rand(random_number_engine);
      }
      for (int i = 0; i < dim - 1; ++i) {
         off_diag[i] = uniform_rand(random_number_engine);
      }
      
      double val_lapack  = 0;
      double val_compnal = 0;
      std::vector<double> vec_lapack(dim);
      std::vector<double> vec_compnal(dim);
      
      auto start = std::chrono::system_clock::now();
      compnal::blas::LapackDstev(&val_lapack, &vec_lapack, diag, off_diag);
      auto time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      double time_sec   = static_cast<double>(time_count)/(1000*1000);
      std::cout << "Time Lapack :" << time_sec << "[sec]: " << val_lapack << std::endl;
      
      
      start = std::chrono::system_clock::now();
      compnal::blas::Dstev(&val_compnal, &vec_compnal, 0, diag, off_diag);
      time_count = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count();
      time_sec   = static_cast<double>(time_count)/(1000*1000);
      std::cout << "Time Compnal:" << time_sec << "[sec]: " << val_compnal << std::endl;
   }

}
