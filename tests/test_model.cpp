//
//  test_model.cpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/05/23.
//

#include "model.hpp"
#include <gtest/gtest.h>

using namespace compnal::model;
using namespace compnal::sparse_matrix;

template<typename RealType>
void CRSTest(const CRS<RealType> &m1, const CRS<RealType> &m2) {
   EXPECT_EQ(m1.GetRowDim(), m2.GetRowDim());
   EXPECT_EQ(m1.GetColDim(), m2.GetColDim());
   
   EXPECT_EQ(m1.GetRow().size(), m2.GetRow().size());
   EXPECT_EQ(m1.GetCol().size(), m2.GetCol().size());
   EXPECT_EQ(m1.GetVal().size(), m2.GetVal().size());
   
   for (std::size_t i = 0; i < m1.GetRow().size(); ++i) {
      EXPECT_EQ(m1.GetRow().at(i), m2.GetRow().at(i));
   }
   
   for (std::size_t i = 0; i < m1.GetCol().size(); ++i) {
      EXPECT_EQ(m1.GetCol().at(i), m2.GetCol().at(i));
   }
   
   for (std::size_t i = 0; i < m1.GetVal().size(); ++i) {
      EXPECT_DOUBLE_EQ(m1.GetVal().at(i), m2.GetVal().at(i));
   }
}

void StateTestSpinOneHalf(const Heisenberg1D<double> &model) {
   
   EXPECT_EQ(model.GetBoundaryCondition()   , BoundaryCondition::OBC);
   EXPECT_EQ(model.GetSystemSize()          , 3                     );
   EXPECT_EQ(model.GetDimOnsite()           , 2                     );
   EXPECT_EQ(model.GetMagnitudeSpin()       , 0.5                   );
   EXPECT_EQ(model.GetTotalSz()             , 0                     );
   EXPECT_EQ(model.GetNumConservedQuantity(), 1                     );
   
   EXPECT_EQ(model.GetJz().size() , 1);
   EXPECT_DOUBLE_EQ(model.GetJz().at(0), 1.0);
   
   EXPECT_EQ(model.GetJxy().size(), 1);
   EXPECT_DOUBLE_EQ(model.GetJxy().at(0), 1.0);
   
   EXPECT_DOUBLE_EQ(model.GetHz(), 0.0);
   EXPECT_DOUBLE_EQ(model.GetDz(), 0.0);
   EXPECT_EQ(model.GetOperatorSx() , CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.5}, {+0.5, +0.0}}));
   EXPECT_EQ(model.GetOperatoriSy(), CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.5}, {-0.5, +0.0}}));
   EXPECT_EQ(model.GetOperatorSz() , CRS<double>(std::vector<std::vector<double>>{{+0.5, +0.0}, {+0.0, -0.5}}));
   EXPECT_EQ(model.GetOperatorSp() , CRS<double>(std::vector<std::vector<double>>{{+0.0, +1.0}, {+0.0, +0.0}}));
   EXPECT_EQ(model.GetOperatorSm() , CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.0}, {+1.0, +0.0}}));
   EXPECT_EQ(model.GetOperatorHam(), CRS<double>(2,2));
   
}

void StateTestSpinOne(const Heisenberg1D<double> &model) {
   
   EXPECT_EQ(model.GetBoundaryCondition()   , BoundaryCondition::OBC);
   EXPECT_EQ(model.GetSystemSize()          , 3                     );
   EXPECT_EQ(model.GetDimOnsite()           , 3                     );
   EXPECT_EQ(model.GetMagnitudeSpin()       , 1.0                   );
   EXPECT_EQ(model.GetTotalSz()             , 0                     );
   EXPECT_EQ(model.GetNumConservedQuantity(), 1                     );
   
   EXPECT_EQ(model.GetJz().size() , 1);
   EXPECT_DOUBLE_EQ(model.GetJz().at(0), 1.0);
   
   EXPECT_EQ(model.GetJxy().size(), 1);
   EXPECT_DOUBLE_EQ(model.GetJxy().at(0), 1.0);
   
   EXPECT_DOUBLE_EQ(model.GetHz(), 0.0);
   EXPECT_DOUBLE_EQ(model.GetDz(), 0.0);
   
   double root2 = std::sqrt(2);
   
   std::vector<std::vector<double>> Sx = {
      {+0.0      , +1.0/root2, +0.0      },
      {+1.0/root2, +0.0      , +1.0/root2},
      {+0.0      , +1.0/root2, +0.0      }
   };
   
   std::vector<std::vector<double>> iSy = {
      {+0.0      , +1.0/root2, +0.0      },
      {-1.0/root2, +0.0      , +1.0/root2},
      {+0.0      , -1.0/root2, +0.0      }
   };
   
   std::vector<std::vector<double>> Sz = {
      {+1.0, 0.0, +0.0},
      {+0.0, 0.0, +0.0},
      {+0.0, 0.0, -1.0}
   };
   
   std::vector<std::vector<double>> Sp = {
      {+0.0, +root2, +0.0  },
      {+0.0, +0.0  , +root2},
      {+0.0, +0.0  , +0.0  }
   };
   
   std::vector<std::vector<double>> Sm = {
      {+0.0  , +0.0  , +0.0},
      {+root2, +0.0  , +0.0},
      {+0.0  , +root2, +0.0}
   };
   
   CRSTest(model.GetOperatorSx() , CRS<double>(Sx));
   CRSTest(model.GetOperatoriSy(), CRS<double>(iSy));
   CRSTest(model.GetOperatorSz() , CRS<double>(Sz));
   CRSTest(model.GetOperatorSp() , CRS<double>(Sp));
   CRSTest(model.GetOperatorSm() , CRS<double>(Sm));
   CRSTest(model.GetOperatorHam(), CRS<double>(3,3));
      
}

TEST(ModelHeisenberg, Constructor1) {
   Heisenberg1D<double> model(3);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor2_1) {
   Heisenberg1D<double> model(3, 0.5);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor2_2) {
   Heisenberg1D<double> model(3, 1.0);
   StateTestSpinOne(model);
}

TEST(ModelHeisenberg, Constructor3) {
   Heisenberg1D<double> model(3, BoundaryCondition::OBC);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor4) {
   Heisenberg1D<double> model(3, 0.5, BoundaryCondition::OBC);
   StateTestSpinOneHalf(model);
   
   Heisenberg1D<double>::CreateOperatorSz(0.5);
}

