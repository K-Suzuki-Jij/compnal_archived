//
//  base_u1_spin_multi_electrons_1d.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/12/20.
//

#ifndef COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_
#define COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_

namespace compnal {
namespace model {

template<typename RealType>
class BaseU1SpinMultiElectrons_1D {
   
   //! @brief Alias of compressed row strage (CRS) with RealType.
   using CRS = sparse_matrix::CRS<RealType>;
   
public:
   
   //! @brief The type of real values.
   using ValueType = RealType;
   
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   BaseU1SpinMultiElectrons_1D() {
      
   }
   
   //! @brief Generate the annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}_{\alpha, \uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{c}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorCUp(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         int num_electron = 0;
         for (int o = 0; o < orbital; ++o) {
            const int basis_electron_onsite = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
            num_electron += CalculateNumElectron(basis_electron_onsite);
         }
         int sign = 1;
         if (num_electron%2 == 1) {
            sign = -1;
         }
         
         for (int col = 0; col < dim_onsite; col++) {
            bool c_1 = true;
            for (int o = 0; o < num_orbital; ++o) {
               const int basis_1 = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
               const int basis_2 = CalculateBasisOnsiteElectron(col, magnitude_lspin, o, num_orbital);
               if (o != orbital && basis_1 != basis_2) {
                  c_1 = false;
               }
            }
            
            int basis_row_electron = CalculateBasisOnsiteElectron(row, magnitude_lspin, orbital, num_orbital);
            int basis_col_electron = CalculateBasisOnsiteElectron(col, magnitude_lspin, orbital, num_orbital);
            
            const bool c_2 = (basis_row_electron == 0 && basis_col_electron == 1);
            const bool c_3 = (basis_row_electron == 2 && basis_col_electron == 3);
            
            if (c_1 && row%dim_onsite_lspin == col%dim_onsite_lspin && (c_2 || c_3)) {
               matrix.val.push_back(sign);
               matrix.col.push_back(col);
            }
      
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}_{\alpha, \downarrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{c}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorCDown(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      
      //--------------------------------
      // # <->  [Cherge  ] -- (N,  2*sz)
      // 0 <->  [        ] -- (0,  0   )
      // 1 <->  [up      ] -- (1,  1   )
      // 2 <->  [down    ] -- (1, -1   )
      // 3 <->  [up&down ] -- (2,  0   )
      //--------------------------------
      
      CRS matrix(dim_onsite, dim_onsite);
      
      for (int row = 0; row < dim_onsite; ++row) {
         int num_electron = 0;
         for (int o = 0; o < orbital; ++o) {
            const int basis_electron_onsite = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
            num_electron += CalculateNumElectron(basis_electron_onsite);
         }
         int sign_1 = 1;
         if (num_electron%2 == 1) {
            sign_1 = -1;
         }
         
         for (int col = 0; col < dim_onsite; col++) {
            bool c_1 = true;
            for (int o = 0; o < num_orbital; ++o) {
               const int basis_1 = CalculateBasisOnsiteElectron(row, magnitude_lspin, o, num_orbital);
               const int basis_2 = CalculateBasisOnsiteElectron(col, magnitude_lspin, o, num_orbital);
               if (o != orbital && basis_1 != basis_2) {
                  c_1 = false;
               }
            }
            
            int basis_row_electron = CalculateBasisOnsiteElectron(row, magnitude_lspin, orbital, num_orbital);
            int basis_col_electron = CalculateBasisOnsiteElectron(col, magnitude_lspin, orbital, num_orbital);
            
            const bool c_2 = (basis_row_electron == 0 && basis_col_electron == 2);
            const bool c_3 = (basis_row_electron == 1 && basis_col_electron == 3);
            
            int sign_2 = 1;
            if (c_3) {
               sign_2 = -1;
            }
            
            if (c_1 && row%dim_onsite_lspin == col%dim_onsite_lspin && (c_2 || c_3)) {
               matrix.val.push_back(sign_1*sign_2);
               matrix.col.push_back(col);
            }
      
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{z}\f$.
   static CRS CreateOnsiteOperatorSzL(const double magnitude_lspin, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; ++row) {
         for (int col = 0; col < dim_onsite; ++col) {
            const int basis_row_lspin = row%dim_onsite_lspin;
            const int basis_col_lspin = col%dim_onsite_lspin;
            const int basis_row_total_electron = row/dim_onsite_lspin;
            const int basis_col_total_electron = col/dim_onsite_lspin;
            if (basis_row_total_electron == basis_col_total_electron && basis_row_lspin == basis_col_lspin) {
               matrix.val.push_back(magnitude_lspin - basis_row_lspin);
               matrix.col.push_back(col);
            }
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{+}\f$.
   static CRS CreateOnsiteOperatorSpL(const double magnitude_lspin, const int num_orbital) {
      const int magnitude_2lspin = utility::DoubleTheNumber(magnitude_lspin);
      const int dim_onsite_lspin = magnitude_2lspin + 1;
      const int dim_onsite = dim_onsite_lspin*static_cast<int>(std::pow(4, num_orbital));
      CRS matrix(dim_onsite, dim_onsite);
      for (int row = 0; row < dim_onsite; ++row) {
         for (int col = 0; col < dim_onsite; ++col) {
            const int a = row%dim_onsite_lspin + 1;
            const int b = col%dim_onsite_lspin + 1;
            const int basis_row_total_electron = row/dim_onsite_lspin;
            const int basis_col_total_electron = col/dim_onsite_lspin;
            if (basis_row_total_electron == basis_col_total_electron && a + 1 == b) {
               matrix.val.push_back(std::sqrt( (magnitude_lspin + 1.0)*(a + b - 1) - a*b ) );
               matrix.col.push_back(col);
            }
         }
         matrix.row.push_back(matrix.col.size());
      }
      return matrix;
   }
   
   //! @brief Generate the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorCUpDagger(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorCDownDagger(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin
   //! \f$ \hat{n}_{\alpha, \uparrow}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha, \uparrow}\f$.
   static CRS CreateOnsiteOperatorNCUp(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin
   //! \f$ \hat{n}_{\alpha, \downarrow}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha, \downarrow}\f$.
   static CRS CreateOnsiteOperatorNCDown(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the number operator for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{n}_{\alpha}=\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{n}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorNC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorNCUp(magnitude_lspin, orbital, num_orbital) + CreateOnsiteOperatorNCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the spin operator for the x-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{x}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //!  + \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{x}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSxC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital) + CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the spin operator for the y-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ i\hat{s}^{y}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //!  - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{s}^{y}_{\alpha}\f$.
   static CRS CreateOnsiteOperatoriSyC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital) - CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the spin operator for the z-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{z}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}
   //!  - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow})\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{z}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSzC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorNCUp(magnitude_lspin, orbital, num_orbital) - CreateOnsiteOperatorNCDown(magnitude_lspin, orbital, num_orbital));
   }
   
   //! @brief Generate the raising operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{+}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{+}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSpC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCUpDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCDown(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the lowering operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{-}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow}\f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{s}^{-}_{\alpha}\f$.
   static CRS CreateOnsiteOperatorSmC(const double magnitude_lspin, const int orbital, const int num_orbital) {
      return CreateOnsiteOperatorCDownDagger(magnitude_lspin, orbital, num_orbital)*CreateOnsiteOperatorCUp(magnitude_lspin, orbital, num_orbital);
   }
   
   //! @brief Generate the spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$
   //! @return The matrix of \f$ \hat{S}^{-}\f$.
   static CRS CreateOnsiteOperatorSmL(const double magnitude_lspin, const int num_orbital) {
      return sparse_matrix::CalculateTransposedMatrix(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ \hat{S}^{x}\f$.
   static CRS CreateOnsiteOperatorSxL(const double magnitude_lspin, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital) + CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate the spin-\f$ S\f$ operator of the local spin
   //!  for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @return The matrix of \f$ i\hat{S}^{y}\f$.
   static CRS CreateOnsiteOperatoriSyL(const double magnitude_lspin, const int num_orbital) {
      return 0.5*(CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital) - CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital));
   }
   
   //! @brief Generate \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}=
   //! \hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}+\hat{s}^{z}_{\alpha}\hat{S}^{z}\f$
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}\f$.
   static CRS CreateOnsiteOperatorSCSL(const double magnitude_lspin, const int orbital, const int num_orbital) {
      const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin, orbital, num_orbital);
      const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin, orbital, num_orbital);
      const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin, orbital, num_orbital);
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin, num_orbital);
      return szc*szl + 0.5*(spc*sml + smc*spl);
   }
   
   //! @brief Generate \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}=
   //! \sum_{\alpha}\left(\hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}
   //!  +\hat{s}^{z}_{\alpha}\hat{S}^{z}\right)\f$
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The matrix of \f$ \hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{S}}\f$.
   static CRS CreateOnsiteOperatorSCSL(const double magnitude_lspin, const int num_orbital) {
      const CRS spl = CreateOnsiteOperatorSpL(magnitude_lspin, num_orbital);
      const CRS sml = CreateOnsiteOperatorSmL(magnitude_lspin, num_orbital);
      const CRS szl = CreateOnsiteOperatorSzL(magnitude_lspin, num_orbital);
      CRS out(szl.row_dim, szl.col_dim);
      for (int o = 0; o < num_orbital; ++o) {
         const CRS spc = CreateOnsiteOperatorSpC(magnitude_lspin, o, num_orbital);
         const CRS smc = CreateOnsiteOperatorSmC(magnitude_lspin, o, num_orbital);
         const CRS szc = CreateOnsiteOperatorSzC(magnitude_lspin, o, num_orbital);
         out = out + szc*szl + 0.5*(spc*sml + smc*spl);
      }
      return out;
   }
   
   
protected:
   //! @brief The annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin, \f$ \hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_c_up_;
   
   //! @brief The annihilation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_c_down_;
   
   //! @brief The creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin \f$ \hat{c}^{\dagger}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_c_up_dagger_;
   
   //! @brief The the creation operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin \f$ \hat{c}^{\dagger}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_c_down_dagger_;
   
   //! @brief The number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the up spin,
   //! \f$ \hat{n}_{\alpha, \uparrow}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_nc_up_;
   
   //! @brief The number operator for the electrons
   //! with the orbital \f$ \alpha \f$ and the down spin,
   //! \f$ \hat{n}_{\alpha, \downarrow}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_nc_down_;
   
   //! @brief The number operator for the electrons with the orbital
   //! \f$ \alpha \f$, \f$ \hat{n}_{\alpha}=\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_nc_;
   
   //! @brief The number operator for the electrons,
   //! \f$ \hat{n}=\sum_{\alpha}\left(\hat{n}_{\alpha, \uparrow} + \hat{n}_{\alpha, \downarrow}\right)\f$.
   CRS onsite_operator_nc_tot_;
   
   //! @brief The spin operator for the x-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{x}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //! + \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   std::vector<CRS> onsite_operator_sxc_;
   
   //! @brief The spin operator for the y-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ i\hat{s}^{y}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}
   //! - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow})\f$.
   //! Here \f$ i=\sqrt{-1}\f$ is the the imaginary unit.
   std::vector<CRS> onsite_operator_isyc_;
   
   //! @brief The spin operator for the z-direction for the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{z}_{\alpha}=\frac{1}{2}(\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \uparrow}
   //! - \hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \downarrow})\f$.
   std::vector<CRS> onsite_operator_szc_;
   
   //! @brief The raising operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{+}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \uparrow}\hat{c}_{\alpha, \downarrow}\f$.
   std::vector<CRS> onsite_operator_spc_;
   
   //! @brief The lowering operator for spin of the electrons with the orbital \f$ \alpha \f$,
   //! \f$ \hat{s}^{-}_{\alpha}=\hat{c}^{\dagger}_{\alpha, \downarrow}\hat{c}_{\alpha, \uparrow}\f$.
   std::vector<CRS> onsite_operator_smc_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the x-direction \f$ \hat{S}^{x}\f$.
   CRS onsite_operator_sxl_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the y-direction \f$ i\hat{S}^{y}\f$ with \f$ i\f$ being the imaginary unit.
   CRS onsite_operator_isyl_;
   
   //! @brief The spin-\f$ S\f$ operator of the local spin for the z-direction \f$ \hat{S}^{z}\f$.
   CRS onsite_operator_szl_;
   
   //! @brief The spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{+}\f$.
   CRS onsite_operator_spl_;
   
   //! @brief The spin-\f$ S\f$ raising operator of the local spin \f$ \hat{S}^{-}\f$.
   CRS onsite_operator_sml_;
   
   //! @brief The correlation between the electron with the orbital \f$ \alpha \f$ spin and local spin
   //! \f$ \hat{\boldsymbol{s}}_{\alpha}\cdot\hat{\boldsymbol{S}}=
   //! \hat{s}^{x}_{\alpha}\hat{S}^{x}+\hat{s}^{y}_{\alpha}\hat{S}^{y}+\hat{s}^{z}_{\alpha}\hat{S}^{z}\f$
   std::vector<CRS> onsite_operator_scsl_;
   
   //! @brief The dimension of the local Hilbert space for the electrons, 4.
   const int dim_onsite_electron_ = 4;
   
   //! @brief The total electron at each orbital \f$ \alpha \f$, \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$.
   std::vector<int> total_electron_ = {0};
   
   //! @brief The number of electron orbitals \f$ n_{\rm o}\f$.
   int num_electron_orbital_ = static_cast<int>(total_electron_.size());
   
   //! @brief The dimension of the local Hilbert space for the total electrons, \f$ 4^{n_{\rm o}}\f$.
   int dim_onsite_total_electron_ = static_cast<int>(std::pow(dim_onsite_electron_, num_electron_orbital_));
   
   //! @brief The system size.
   int system_size_ = 0;
   
   //! @brief Twice the number of the total sz
   //! \f$ 2\langle\hat{S}^{z}_{\rm tot}\rangle =
   //! 2\sum^{N}_{i=1}\left(\hat{S}^{z}_{i} + \sum_{\alpha}\hat{s}^{z}_{i,\alpha}\right)\f$.
   int total_2sz_ = 0;
   
   //! @brief The magnitude of the local spin \f$ S\f$.
   int magnitude_2lspin_ = 1;
   
   //! @brief The dimension of the local Hilbert space, \f$ 2S + 1\f$.
   int dim_onsite_lspin_ = magnitude_2lspin_ + 1;
   
   //! @brief The dimension of the local Hilbert space, \f$ 4^{n_{\rm o}}\times (2S + 1) \f$.
   int dim_onsite_ = dim_onsite_total_electron_*dim_onsite_lspin_;
   
   //! @brief The calculated eigenvectors and eigenvalues.
   std::unordered_set<int> calculated_eigenvector_set_;
   
   //! @brief Bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons at each orbital \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! The first value of std::vector<int> stores twice the number of the total sz and remaining values correspond to the orbitals of the electrons.
   std::unordered_map<std::vector<int>, std::vector<std::int64_t>, utility::VectorHash> bases_;
   
   //! @brief Inverse bases of the target Hilbert space specified by
   //! the system size \f$ N\f$, the magnitude of the local spin \f$ S\f$,
   //! the number of the total electrons at each orbital \f$ \langle\hat{N}_{{\rm e}, \alpha}\rangle\f$,
   //! and the total sz \f$ \langle\hat{S}^{z}_{\rm tot}\rangle \f$.
   //! The first value of std::vector<int> stores twice the number of the total sz and remaining values correspond to the orbitals of the electrons.
   std::unordered_map<std::vector<int>, std::unordered_map<std::int64_t, std::int64_t>, utility::VectorHash> bases_inv_;
      
   //! @brief Calculate onsite basis for the electrons from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @return The onsite basis for the electrons.
   inline int CalculateBasisOnsiteElectron(const int basis_onsite, const int orbital) const {
      const int num_inner_electron = num_electron_orbital_ - orbital - 1;
      int basis_electron_onsite = basis_onsite/dim_onsite_lspin_;
      for (int i = 0; i < num_inner_electron; ++i) {
         basis_electron_onsite /= dim_onsite_electron_;
      }
      return basis_electron_onsite%dim_onsite_electron_;
   }
   
   //! @brief Calculate onsite basis for the electrons from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @param magnitude_lspin The magnitude of the local spin \f$ S \f$.
   //! @param orbital The electron orbital \f$ \alpha \f$.
   //! @param num_orbital The number of the orbitals of the electrons \f$ n_{\rm o}\f$.
   //! @return The onsite basis for the electrons.
   inline static int CalculateBasisOnsiteElectron(const int basis_onsite, const double magnitude_lspin, const int orbital, const int num_orbital) {
      const int dim_onsite_lspin = utility::DoubleTheNumber(magnitude_lspin) + 1;
      const int num_inner_electron = num_orbital - orbital - 1;
      int dim_onsite_electron = 1;
      for (int i = 0; i < num_orbital; ++i) {
         dim_onsite_electron *= 4;
      }
      int basis_electron_onsite = basis_onsite/dim_onsite_lspin;
      for (int i = 0; i < num_inner_electron; ++i) {
         basis_electron_onsite /= dim_onsite_electron;
      }
      return basis_electron_onsite%dim_onsite_electron;
   }
   
   //! @brief Calculate onsite basis for the loca spins from an onsite basis.
   //! @param basis_onsite The onsite basis.
   //! @return The onsite basis for the local spins.
   inline int CalculateBasisOnsiteLSpin(const int basis_onsite) const {
      return basis_onsite%dim_onsite_lspin_;
   }
   
   //! @brief Calculate the number of electrons from the onsite electron basis.
   //! @param basis_electron_onsite The onsite electron basis.
   //! @return The number of electrons
   inline static int CalculateNumElectron(const int basis_electron_onsite) {
      if (basis_electron_onsite == 0) {
         return 0;
      }
      else if (basis_electron_onsite == 1 || basis_electron_onsite == 2) {
         return 1;
      }
      else if (basis_electron_onsite == 3) {
         return 2;
      }
      else {
         std::stringstream ss;
         ss << "Unknown error detected in " << __FUNCTION__ << std::endl;
         throw std::runtime_error(ss.str());
      }
   }
   
};

}
}


#endif /* COMPNAL_MODEL_BASE_U1_SPIN_MULTI_ELECTRONS_1D_HPP_ */
