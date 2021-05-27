//
//  ed_heisengerg_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/27.
//

#ifndef ed_heisengerg_1d_hpp
#define ed_heisengerg_1d_hpp

#include "sparse_matrix.hpp"
#include "model.hpp"
#include <unordered_map>
#include <cstdint>

namespace compnal {
namespace solver {

template<typename RealType>
class EdHeisenberg1D {
  
   using CRS = sparse_matrix::CRS<RealType>;
   
   using BraketVector = sparse_matrix::BraketVector<RealType>;
   
public:
   const model::Heisenberg1D<RealType> model;
   
   explicit EdHeisenberg1D(const model::Heisenberg1D<RealType> model_in): model(model_in) {
      SetDim();
   }
   
   int64_t GetDim() const {
      int system_size   = model.GetSystemSize();
      int total_2sz     = model.GetTotal2Sz();
      int max_total_2sz = system_size*model.GetMagnitude2Spin();
      if (dim_[system_size - 1][(total_2sz + max_total_2sz)/2] < 0) {
         throw std::runtime_error("Maybe too large Hilbert space. Overflow detected.");
      }
      return dim_[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   int64_t GetDim(double target_sz) const {
      model.CheckInteger(target_sz);
      int system_size   = model.GetSystemSize();
      int total_2sz     = static_cast<int>(target_sz*2);
      int max_total_2sz = system_size*model.GetMagnitude2Spin();
      if (dim_[system_size - 1][(total_2sz + max_total_2sz)/2] < 0) {
         throw std::runtime_error("Maybe too large Hilbert space. Overflow detected.");
      }
      if (std::abs(total_2sz) > max_total_2sz) {
         return 0;
      }
      return dim_[system_size - 1][(total_2sz + max_total_2sz)/2];
   }
   
   
private:
   std::vector<int64_t> basis_;
   BraketVector gs_vector_;
   std::vector<std::vector<int64_t>> dim_;
   
   void SetDim() {
      int system_size = model.GetSystemSize();
      int magnitude_2spin = model.GetMagnitude2Spin();
      int max_total_2sz = system_size*magnitude_2spin;
      dim_.assign(system_size, std::vector<int64_t>(max_total_2sz + 1));
      
      for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
         dim_[0][(s + magnitude_2spin)/2] = 1;
      }
      
      for (int site = 1; site < system_size; site++) {
         for (int s = -magnitude_2spin; s <= magnitude_2spin; s += 2) {
            for (int s_prev = -magnitude_2spin*site; s_prev <= magnitude_2spin*site; s_prev += 2) {
               dim_[site][(s + s_prev + magnitude_2spin*(site + 1))/2] += dim_[site - 1][(s_prev + magnitude_2spin*site)/2];
            }
         }
      }
   }
   
   
};

} // namespace solver
} // namespace compnal


#endif /* ed_heisengerg_1d_hpp */
