//
//  test_model.hpp
//  compnal
//
//  Created by Kohei Suzuki on 2021/06/14.
//

#ifndef test_model_hpp
#define test_model_hpp

#include "model.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

template<typename RealType>
void CRSTest(const compnal::sparse_matrix::CRS<RealType> &m1, const compnal::sparse_matrix::CRS<RealType> &m2) {
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

void StateTestSpinOneHalf(const compnal::model::Heisenberg1D<double> &model) {
   
   EXPECT_EQ(model.GetBoundaryCondition()   , compnal::model::BoundaryCondition::OBC);
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
   EXPECT_EQ(model.GetOperatorSx() , compnal::sparse_matrix::CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.5}, {+0.5, +0.0}}));
   EXPECT_EQ(model.GetOperatoriSy(), compnal::sparse_matrix::CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.5}, {-0.5, +0.0}}));
   EXPECT_EQ(model.GetOperatorSz() , compnal::sparse_matrix::CRS<double>(std::vector<std::vector<double>>{{+0.5, +0.0}, {+0.0, -0.5}}));
   EXPECT_EQ(model.GetOperatorSp() , compnal::sparse_matrix::CRS<double>(std::vector<std::vector<double>>{{+0.0, +1.0}, {+0.0, +0.0}}));
   EXPECT_EQ(model.GetOperatorSm() , compnal::sparse_matrix::CRS<double>(std::vector<std::vector<double>>{{+0.0, +0.0}, {+1.0, +0.0}}));
   EXPECT_EQ(model.GetOperatorHam(), compnal::sparse_matrix::CRS<double>(2,2));
   
}

void StateTestSpinOne(const compnal::model::Heisenberg1D<double> &model) {
   
   EXPECT_EQ(model.GetBoundaryCondition()   , compnal::model::BoundaryCondition::OBC);
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
   
   CRSTest(model.GetOperatorSx() , compnal::sparse_matrix::CRS<double>(Sx));
   CRSTest(model.GetOperatoriSy(), compnal::sparse_matrix::CRS<double>(iSy));
   CRSTest(model.GetOperatorSz() , compnal::sparse_matrix::CRS<double>(Sz));
   CRSTest(model.GetOperatorSp() , compnal::sparse_matrix::CRS<double>(Sp));
   CRSTest(model.GetOperatorSm() , compnal::sparse_matrix::CRS<double>(Sm));
   CRSTest(model.GetOperatorHam(), compnal::sparse_matrix::CRS<double>(3,3));
      
}

TEST(ModelHeisenberg, Constructor1) {
   compnal::model::Heisenberg1D<double> model(3);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor2_1) {
   compnal::model::Heisenberg1D<double> model(3, 0.5);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor2_2) {
   compnal::model::Heisenberg1D<double> model(3, 1.0);
   StateTestSpinOne(model);
}

TEST(ModelHeisenberg, Constructor3) {
   compnal::model::Heisenberg1D<double> model(3, compnal::model::BoundaryCondition::OBC);
   StateTestSpinOneHalf(model);
}

TEST(ModelHeisenberg, Constructor4) {
   compnal::model::Heisenberg1D<double> model(3, 0.5, compnal::model::BoundaryCondition::OBC);
   StateTestSpinOneHalf(model);
   
   compnal::model::Heisenberg1D<double>::CreateOperatorSz(0.5);
}


#endif /* test_model_hpp */
