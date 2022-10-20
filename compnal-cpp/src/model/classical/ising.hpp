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
//  ising.hpp
//  compnal
//
//  Created by kohei on 2022/08/13.
//  
//

#ifndef COMPNAL_MODEL_ISING_HPP_
#define COMPNAL_MODEL_ISING_HPP_

#include "../../lattice/all.hpp"
#include "../../utility/type.hpp"
#include "../../interaction/quadratic_any_interaction.hpp"
#include <vector>

namespace compnal {
namespace model {

//! @brief Class for representing the classical Ising models.
//! @tparam LatticeType The lattice type.
//! @tparam RealType The value type, which must be floating point type.
template<class LatticeType, typename RealType>
class Ising {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   //! @brief The value type.
   using ValueType = RealType;
   
   //! @brief The index type.
   using IndexType = std::int32_t;
   
   //! @brief The operator type, which here is the type of Ising spin -1 or +1.
   using OPType = utility::SpinType;
   
   //! @brief The linear interaction type.
   using LinearType = RealType;
   
   //! @brief The quadratic interaction type.
   using QuadraticType = RealType;
   
   //! @brief Constructor for Ising class.
   //! @param linear The linear interaction.
   //! @param quadratic The quadratic interaction.
   Ising(const LatticeType &lattice,
         const LinearType linear,
         const QuadraticType quadratic):
   lattice_(lattice), linear_(linear), quadratic_(quadratic) {}
   
   //! @brief Get linear interaction.
   //! @return The linear interaction.
   LinearType GetLinear() const {
      return linear_;
   }
   
   //! @brief Get quadratic interaction.
   //! @return The quadratic interaction.
   QuadraticType GetQuadratic() const {
      return quadratic_;
   }
   
   //! @brief Get the system size.
   //! @return The system size.
   std::int32_t GetSystemSize() const {
      return lattice_.GetSystemSize();
   }
   
   //! @brief Get the boundary condition.
   //! @return The boundary condition.
   lattice::BoundaryCondition GetBoundaryCondition() const {
      return lattice_.GetBoundaryCondition();
   }
   
   //! @brief Generate index list.
   //! @return The index list.
   auto GenerateIndexList() const {
      return lattice_.GenerateIndexList();
   }
   
   //! @brief Get the degree of the interactions.
   //! @return The degree.
   std::int32_t GetDegree() const {
      if (std::abs(quadratic_) > std::numeric_limits<RealType>::epsilon()) {
         return 2;
      }
      else if (std::abs(linear_) > std::numeric_limits<RealType>::epsilon()) {
         return 1;
      }
      else {
         return 0;
      }
   }
   
   //! @brief Get the lattice.
   //! @return The lattice.
   const LatticeType &GetLattice() const {
      return lattice_;
   }
   
   //! @brief Calculate energy corresponding to the spin configuration.
   //! @param spin_configuration The spin configuration.
   //! @return The energy.
   RealType CalculateEnergy(const std::vector<OPType> &spin_configuration) const {
      return CalculateEnergy(lattice_, spin_configuration);
   }
   
   //! @brief Calculate n-th moment for the list of spin configuration.
   //! @param spin_configurations The list of spin configuration.
   //! @param degree The degree of the moment.
   //! @param num_threads The number of threads.
   //! @return The moment.
   RealType CalculateMoment(const std::vector<std::vector<OPType>> &spin_configurations,
                            const std::int32_t degree,
                            const std::int32_t num_threads = utility::DEFAULT_NUM_THREADS) const {
      if (degree <= 0) {
         throw std::runtime_error("degree must be lager than 0.");
      }
      
      RealType val = 0;
#pragma omp parallel for schedule(guided) reduction(+: val) num_threads(num_threads)
      for (std::int32_t i = 0; i < static_cast<std::int32_t>(spin_configurations.size()); ++i) {
         RealType avg = 0;
         for (std::size_t j = 0; j < spin_configurations[i].size(); ++j) {
            avg += spin_configurations[i][j];
         }
         
         RealType prod = 1;
         for (std::int32_t j = 0; j < degree; ++j) {
            prod = prod*avg;
         }
         val += prod;
      }
      return val/spin_configurations.size();
   }
   
   //! @brief Calculate the on-site expectation value.
   //! @param spin_configurations The list of spin configuration.
   //! @param index The index of the spin.
   //! @return The on-site expectation value.
   RealType CalculateOnsiteAverage(const std::vector<std::vector<OPType>> &spin_configurations,
                                   const IndexType index) const {
      if (index < 0) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::size_t i = 0; i < spin_configurations.size(); ++i) {
         val += spin_configurations[i][index];
      }
      return val/spin_configurations.size();
   }
   
   //! @brief Calculate the spin-spin correlation.
   //! @param spin_configurations The list of spin configuration.
   //! @param index1 The index of the spin.
   //! @param index2 The index of the spin.
   //! @return The spin-spin correlation.
   RealType CalculateCorrelation(const std::vector<std::vector<OPType>> &spin_configurations,
                                 const IndexType index1,
                                 const IndexType index2) const {
      const std::int32_t size = static_cast<std::int32_t>(spin_configurations.size());
      if (index1 < 0 || index2 < 0 || index1 >= size || index2 >= size) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::int32_t i = 0; i < size; ++i) {
         val += spin_configurations[i][index1]*spin_configurations[i][index2];
      }
      return val/size;
   }
   
private:
   //! @brief The linear interaction.
   const LatticeType lattice_;
   
   //! @brief The linear interaction.
   const LinearType linear_ = 0;
   
   //! @brief The quadratic interaction.
   const QuadraticType quadratic_ = 0;
   
   //! @brief Calculate energy corresponding to the spin configuration on the one-dimensional chain.
   //! @param lattice The one-dimensional chain.
   //! @param spin_configuration The spin configuration.
   //! @return The energy.
   RealType CalculateEnergy(const lattice::Chain &lattice,
                            const std::vector<OPType> &spin_configuration) const {
      RealType energy = 0;
      const std::int32_t system_size = lattice.GetSystemSize();
      
      if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
         for (std::int32_t index = 0; index < system_size - 1; ++index) {
            energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
         }
         energy += quadratic_*spin_configuration[system_size - 1]*spin_configuration[0] + linear_*spin_configuration[system_size - 1];
      }
      else if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
         for (std::int32_t index = 0; index < system_size - 1; ++index) {
            energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
         }
         energy += linear_*spin_configuration[system_size - 1];
      }
      else {
         throw std::runtime_error("Unsupported BinaryCondition");
      }
      
      return energy;
   }
   
   //! @brief Calculate energy corresponding to the spin configuration on the two-dimensional square lattice.
   //! @param lattice The two-dimensional square lattice.
   //! @param spin_configuration The spin configuration.
   //! @return The energy.
   RealType CalculateEnergy(const lattice::Square &lattice,
                            const std::vector<OPType> &spin_configuration) const {
      RealType energy = 0;
      const std::int32_t x_size = lattice.GetXSize();
      const std::int32_t y_size = lattice.GetYSize();
      
      if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
               const std::int32_t index = coo_y*x_size + coo_x;
               energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
            }
            energy += quadratic_*spin_configuration[coo_y*x_size + x_size - 1]*spin_configuration[coo_y*x_size + 0] + linear_*spin_configuration[coo_y*x_size + x_size - 1];
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
               const std::int32_t index = coo_y*x_size + coo_x;
               const std::int32_t index_p1 = (coo_y + 1)*x_size + coo_x;
               energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
            }
            energy += quadratic_*spin_configuration[(y_size - 1)*x_size + coo_x]*spin_configuration[0*x_size + coo_x] + linear_*spin_configuration[(y_size - 1)*x_size + coo_x];
         }
         
      }
      else if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
               const std::int32_t index = coo_y*x_size + coo_x;
               energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
            }
            energy += linear_*spin_configuration[coo_y*x_size + x_size - 1];
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
               const std::int32_t index = coo_y*x_size + coo_x;
               const std::int32_t index_p1 = (coo_y + 1)*x_size + coo_x;
               energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
            }
            energy += linear_*spin_configuration[(y_size - 1)*x_size + coo_x];
         }
      }
      else {
         throw std::runtime_error("Unsupported BinaryCondition");
      }
      
      return energy;
   }
   
   //! @brief Calculate energy corresponding to the spin configuration on the three-dimensional cubic lattice.
   //! @param lattice The three-dimensional cubic lattice.
   //! @param spin_configuration The spin configuration.
   //! @return The energy.
   RealType CalculateEnergy(const lattice::Cubic &lattice,
                            const std::vector<OPType> &spin_configuration) const {
      
      RealType energy = 0;
      const std::int32_t x_size = lattice.GetXSize();
      const std::int32_t y_size = lattice.GetYSize();
      const std::int32_t z_size = lattice.GetZSize();
      
      if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::PBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
               for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
               }
               energy +=
               quadratic_*spin_configuration[coo_z*x_size*y_size + coo_y*x_size + x_size - 1]*spin_configuration[coo_z*x_size*y_size + coo_y*x_size + 0] +
               linear_*spin_configuration[coo_z*x_size*y_size + coo_y*x_size + x_size - 1];
            }
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
               for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  const std::int32_t index_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
               }
               energy +=
               quadratic_*spin_configuration[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x]*spin_configuration[coo_z*x_size*y_size + 0*x_size + coo_x] +
               linear_*spin_configuration[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x];
            }
         }
         
         // z-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
               for (std::int32_t coo_z = 0; coo_z < z_size - 1; ++coo_z) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  const std::int32_t index_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
               }
               energy += quadratic_*spin_configuration[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x]*spin_configuration[0*x_size*y_size + coo_y*x_size + coo_x] +
               linear_*spin_configuration[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x];
            }
         }
      }
      else if (lattice.GetBoundaryCondition() == lattice::BoundaryCondition::OBC) {
         // x-direction
         for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
            for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
               for (std::int32_t coo_x = 0; coo_x < x_size - 1; ++coo_x) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index + 1] + linear_*spin_configuration[index];
               }
               energy += linear_*spin_configuration[coo_z*x_size*y_size + coo_y*x_size + x_size - 1];
            }
         }
         
         // y-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_z = 0; coo_z < z_size; ++coo_z) {
               for (std::int32_t coo_y = 0; coo_y < y_size - 1; ++coo_y) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  const std::int32_t index_p1 = coo_z*x_size*y_size + (coo_y + 1)*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
               }
               energy += linear_*spin_configuration[coo_z*x_size*y_size + (y_size - 1)*x_size + coo_x];
            }
         }
         
         // z-direction
         for (std::int32_t coo_x = 0; coo_x < x_size; ++coo_x) {
            for (std::int32_t coo_y = 0; coo_y < y_size; ++coo_y) {
               for (std::int32_t coo_z = 0; coo_z < z_size - 1; ++coo_z) {
                  const std::int32_t index = coo_z*x_size*y_size + coo_y*x_size + coo_x;
                  const std::int32_t index_p1 = (coo_z + 1)*x_size*y_size + coo_y*x_size + coo_x;
                  energy += quadratic_*spin_configuration[index]*spin_configuration[index_p1] + linear_*spin_configuration[index];
               }
               energy += linear_*spin_configuration[(z_size - 1)*x_size*y_size + coo_y*x_size + coo_x];
            }
         }
      }
      else {
         throw std::runtime_error("Unsupported BoundaryCondition");
      }
      
      return energy;
   }
   
   //! @brief Calculate energy corresponding to the spin configuration on the infinite range lattice.
   //! @param lattice The infinite range lattice.
   //! @param spin_configuration The spin configuration.
   //! @return The energy.
   RealType CalculateEnergy(const lattice::InfiniteRange &lattice,
                            const std::vector<OPType> &spin_configuration) const {
      RealType energy = 0;
      const std::int32_t system_size = lattice.GetSystemSize();
      for (std::int32_t i = 0; i < system_size; ++i) {
         energy += linear_*spin_configuration[i];
         for (std::int32_t j = i + 1; j < system_size; ++j) {
            energy += quadratic_*spin_configuration[i]*spin_configuration[j];
         }
      }
      return energy;
   }
   
};

template<typename RealType>
class Ising<lattice::AnyLattice, RealType> {
   static_assert(std::is_floating_point<RealType>::value, "Template parameter RealType must be floating point type");
   
public:
   using ValueType = RealType;
   using OPType = utility::SpinType;
   using IndexType = typename interaction::QuadraticAnyInteraction<RealType>::IndexType;
   using IndexHash = typename interaction::QuadraticAnyInteraction<RealType>::IndexHash;
   using LinearType = typename interaction::QuadraticAnyInteraction<RealType>::LinearType;
   using QuadraticType = typename interaction::QuadraticAnyInteraction<RealType>::QuadraticType;
   
   Ising(const lattice::AnyLattice &lattice,
         const LinearType &linear,
         const QuadraticType &quadratic):
   lattice_(lattice), interaction_(linear, quadratic) {}
   
   const std::vector<IndexType> &GenerateIndexList() const {
      return interaction_.GetIndexList();
   }
   
   RealType GetConstant() const {
      return interaction_.GetConstant();
   }
   
   const std::vector<RealType> &GetLinear() const {
      return interaction_.GetLinear();
   }
   
   const std::unordered_map<IndexType, std::int32_t, IndexHash> &GetIndexMap() const {
      return interaction_.GetIndexMap();
   }
   
   const std::vector<std::int64_t> &GetRowPtr() const {
      return interaction_.GetRowPtr();
   }
   
   const std::vector<std::int32_t> &GetColPtr() const {
      return interaction_.GetColPtr();
   }
   
   const std::vector<RealType> &GetValPtr() const {
      return interaction_.GetValPtr();
   }
   
   std::int32_t GetSystemSize() const {
      return interaction_.GetSystemSize();
   }
   
   std::int32_t GetDegree() const {
      return interaction_.GetDegree();
   }
   
   lattice::BoundaryCondition GetBoundaryCondition() const {
      return lattice_.GetBoundaryCondition();
   }
   
   RealType CalculateEnergy(const std::vector<OPType> &spin_configuration) const {
      RealType val = interaction_.GetConstant();
      return val;
   }
   
   RealType CalculateMoment(const std::vector<std::vector<OPType>> &spin_configurations,
                            const std::int32_t degree,
                            const std::int32_t num_threads = utility::DEFAULT_NUM_THREADS) const {
      if (degree <= 0) {
         throw std::runtime_error("degree must be lager than 0.");
      }
      
      RealType val = 0;
#pragma omp parallel for schedule(guided) reduction(+: val) num_threads(num_threads)
      for (std::int32_t i = 0; i < static_cast<std::int32_t>(spin_configurations.size()); ++i) {
         RealType avg = CalculateMagnetization(spin_configurations[i]);
         RealType prod = 1;
         for (std::int32_t j = 0; j < degree; ++j) {
            prod = prod*avg;
         }
         val += prod;
      }
      return val/spin_configurations.size();
   }
   
   RealType CalculateCorrelation(const std::vector<std::vector<OPType>> &spin_configurations,
                                 const IndexType ind1,
                                 const IndexType ind2) const {
      const std::unordered_map<IndexType, std::int32_t, IndexHash> &index_map = interaction_.GetIndexMap();
      if (index_map.count(ind1) == 0 || index_map.count(ind2) == 0) {
         throw std::runtime_error("The index is out of range.");
      }
      RealType val = 0;
      for (std::size_t i = 0; i < spin_configurations.size(); ++i) {
         val += spin_configurations[i][index_map.at(ind1)]*spin_configurations[i][index_map.at(ind2)];
      }
      return val/spin_configurations.size();
   }
   
private:
   const interaction::QuadraticAnyInteraction<RealType> interaction_;
   const lattice::AnyLattice lattice_;
   
   RealType CalculateMagnetization(const std::vector<OPType> &spin_configuration) const {
      RealType val = 0;
      for (std::size_t i = 0; i < spin_configuration.size(); ++i) {
         val += spin_configuration[i];
      }
      return val/spin_configuration.size();
   }
   
};

template<class LatticeType, typename RealType>
auto make_ising(const LatticeType &lattice,
                const typename Ising<LatticeType, RealType>::LinearType &linear,
                const typename Ising<LatticeType, RealType>::QuadraticType &quadratic) {
   return Ising<LatticeType, RealType>{lattice, linear, quadratic};
}

} // namespace model
} // namespace compnal

#endif /* COMPNAL_MODEL_ISING_HPP_ */
