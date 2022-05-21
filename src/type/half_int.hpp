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
//  Created by Kohei Suzuki on 2022/01/30.
//

#ifndef COMPNAL_TYPE_HALF_INT_HPP_
#define COMPNAL_TYPE_HALF_INT_HPP_

#include <cmath>
#include <stdexcept>

namespace compnal {
namespace type {

//! @brief Class to represent half-integer.
class HalfInt {
  public:
   //------------------------------------------------------------------
   //---------------------------Constructors---------------------------
   //------------------------------------------------------------------
   //! @brief Constructor of HalfInt class.
   HalfInt(){};

   //! @brief Constructor of HalfInt class.
   //! @tparam ValueType Value type.
   //! @param value The value to be assigned. It must be half-integer.
   template <typename ValueType>
   HalfInt(const ValueType value) {
      if (std::floor(2 * value) != 2 * value) {
         std::stringstream ss;
         ss << "The input number " << value << " is not half-integer." << std::endl;
         throw std::runtime_error(ss.str());
      }
      integer_ = static_cast<int>(2 * value);
   };

   //------------------------------------------------------------------
   //----------------------Public Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Get integer (two times the actual value).
   int GetInteger() const { return integer_; }

   //------------------------------------------------------------------
   //-----------------------Operator Overloading-----------------------
   //------------------------------------------------------------------
   //! @brief Operator overloading: unary plus operator.
   HalfInt operator+() const { return *this; }

   //! @brief Operator overloading: unary negation operator.
   HalfInt operator-() const { return HalfInt(-0.5 * integer_); }

   //! @brief Operator overloading: compound assignment plus operator.
   //! @tparam ValueType Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template <typename ValueType>
   HalfInt &operator+=(const ValueType rhs) {
      return *this = *this + rhs;
   }

   //! @brief Operator overloading: compound assignment subtraction operator.
   //! @tparam ValueType Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template <typename ValueType>
   HalfInt &operator-=(const ValueType rhs) {
      return *this = *this - rhs;
   }

   //! @brief Operator overloading: compound assignment multiplication operator.
   //! @tparam ValueType Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template <typename ValueType>
   HalfInt &operator*=(const ValueType rhs) {
      return *this = *this * rhs;
   }

   //! @brief Operator overloading: compound assignment division operator.
   //! @tparam ValueType Value type of the right-hand side.
   //! @param rhs The value of the right-hand side.
   template <typename ValueType>
   HalfInt &operator/=(const ValueType rhs) {
      return *this = *this / rhs;
   }

   //! @brief Operator overloading: casting operator to ValueType.
   //! @tparam ValueType Value type.
   template <typename ValueType>
   operator ValueType() const noexcept {
      return static_cast<ValueType>(0.5 * integer_);
   }

  private:
   //------------------------------------------------------------------
   //---------------------Private Member variables---------------------
   //------------------------------------------------------------------
   //! @brief Store two times the actual half-integer as int.
   int integer_ = 0;
};

//------------------------------------------------------------------
//-------------Operator overloading: Addition Operator--------------
//------------------------------------------------------------------
//! @brief Operator overloading: addition operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator+(const ValueType lhs, const HalfInt &rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{lhs + 0.5 * rhs.GetInteger()};
   } else {
      return lhs + ValueType{0.5} * rhs.GetInteger();
   }
}

//! @brief Operator overloading: addition operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator+(const HalfInt &lhs, const ValueType rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{0.5 * lhs.GetInteger() + rhs};
   } else {
      return ValueType{0.5} * lhs.GetInteger() + rhs;
   }
}

//! @brief Operator overloading: addition operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
HalfInt operator+(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt{0.5 * (lhs.GetInteger() + rhs.GetInteger())};
}

//------------------------------------------------------------------
//-----------Operator overloading: Subtraction Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: subtraction operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator-(const ValueType lhs, const HalfInt &rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{lhs - 0.5 * rhs.GetInteger()};
   } else {
      return lhs - ValueType{0.5} * rhs.GetInteger();
   }
}

//! @brief Operator overloading: subtraction operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator-(const HalfInt &lhs, const ValueType rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{0.5 * lhs.GetInteger() - rhs};
   } else {
      return ValueType{0.5} * lhs.GetInteger() - rhs;
   }
}

//! @brief Operator overloading: subtraction operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
HalfInt operator-(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt{0.5 * (lhs.GetInteger() - rhs.GetInteger())};
}

//------------------------------------------------------------------
//----------Operator overloading: Multiplication Operator-----------
//------------------------------------------------------------------
//! @brief Operator overloading: multiplication operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator*(const ValueType lhs, const HalfInt &rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{lhs * 0.5 * rhs.GetInteger()};
   } else {
      return lhs * ValueType{0.5} * rhs.GetInteger();
   }
}

//! @brief Operator overloading: multiplication operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator*(const HalfInt &lhs, const ValueType rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return HalfInt{0.5 * lhs.GetInteger() * rhs};
   } else {
      return ValueType{0.5} * lhs.GetInteger() * rhs;
   }
}

//! @brief Operator overloading: multiplication operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
double operator*(const HalfInt &lhs, const HalfInt &rhs) { return 0.5 * lhs.GetInteger() * 0.5 * rhs.GetInteger(); }

//------------------------------------------------------------------
//-------------Operator overloading: Division Operator--------------
//------------------------------------------------------------------
//! @brief Operator overloading: division operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator/(const ValueType lhs, const HalfInt &rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return double{lhs / (0.5 * rhs.GetInteger())};
   } else {
      return lhs / (ValueType{0.5} * rhs.GetInteger());
   }
}

//! @brief Operator overloading: division operator.
//! @tparam ValueType Value type of the left-hand side.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
auto operator/(const HalfInt &lhs, const ValueType rhs) {
   if constexpr (std::is_integral<ValueType>::value) {
      return double{0.5 * lhs.GetInteger() / rhs};
   } else {
      return ValueType{0.5} * lhs.GetInteger() / rhs;
   }
}

//! @brief Operator overloading: division operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
double operator/(const HalfInt &lhs, const HalfInt &rhs) {
   return static_cast<double>(lhs.GetInteger()) / rhs.GetInteger();
}

//------------------------------------------------------------------
//-------------Operator overloading: Equality Operator--------------
//------------------------------------------------------------------
//! @brief Operator overloading: equality operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator==(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() == rhs.GetInteger(); }

//! @brief Operator overloading: equality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator==(const ValueType lhs, const HalfInt &rhs) {
   return lhs == 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: equality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator==(const HalfInt &lhs, const ValueType rhs) {
   return 0.5 * lhs.GetInteger() == rhs;
}

//------------------------------------------------------------------
//------------Operator overloading: Inequality Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: inequality operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator!=(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() != rhs.GetInteger(); }

//! @brief Operator overloading: inequality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator!=(const ValueType lhs, const HalfInt &rhs) {
   return lhs != 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: inequality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator!=(const HalfInt &lhs, const ValueType rhs) {
   return 0.5 * lhs.GetInteger() != rhs;
}

//------------------------------------------------------------------
//------------Operator overloading: Comparison Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: comparison operator (less than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator<(const ValueType lhs, const HalfInt &rhs) {
   return lhs < 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (less than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator<(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5 * lhs.GetInteger() < rhs;
}

//! @brief Operator overloading: comparison operator (less than).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator<(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() < rhs.GetInteger(); }

//! @brief Operator overloading: comparison operator (greater than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator>(const ValueType lhs, const HalfInt &rhs) {
   return lhs > 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator>(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5 * lhs.GetInteger() > rhs;
}

//! @brief Operator overloading: comparison operator (greater than).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator>(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() > rhs.GetInteger(); }

//! @brief Operator overloading: comparison operator (less than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator<=(const ValueType lhs, const HalfInt &rhs) {
   return lhs <= 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (less than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator<=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5 * lhs.GetInteger() <= rhs;
}

//! @brief Operator overloading: comparison operator (less than or equal).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator<=(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() <= rhs.GetInteger(); }

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator>=(const ValueType lhs, const HalfInt &rhs) {
   return lhs >= 0.5 * rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template <typename ValueType>
bool operator>=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5 * lhs.GetInteger() >= rhs;
}

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator>=(const HalfInt &lhs, const HalfInt &rhs) { return lhs.GetInteger() >= rhs.GetInteger(); }

//------------------------------------------------------------------
//----------------Operator overloading: I/O Stream------------------
//------------------------------------------------------------------
//! @brief Operator overloading: output operator.
//! @param os Ostream object.
//! @param half_int Half-integer.
std::ostream &operator<<(std::ostream &os, const HalfInt &half_int) {
   os << half_int.GetInteger() * 0.5;
   return os;
}

}  // namespace type
}  // namespace compnal

#endif /* COMPNAL_TYPE_HALF_INT_HPP_ */
