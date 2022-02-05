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
   HalfInt() {};
   
   //! @brief Constructor of HalfInt class.
   //! @tparam ValueType Value type.
   //! @param value The value to be assigned. It must be half-integer.
   template<typename ValueType>
   HalfInt(ValueType value) {
      value = 2*value;
      if (std::floor(value) != value) {
         throw std::runtime_error("The input number is not half-integer");
      }
      integer_ = static_cast<int>(value);
   };
   
   //------------------------------------------------------------------
   //----------------------Public Member Functions---------------------
   //------------------------------------------------------------------
   //! @brief Get integer (two times the actual value).
   int GetInteger() const {
      return integer_;
   }
   
   //------------------------------------------------------------------
   //-----------------------Operator Overloading-----------------------
   //------------------------------------------------------------------
   //! @brief Operator overloading: unary plus operator.
   HalfInt operator+() const {
      return *this;
   }
   
   //! @brief Operator overloading: unary negation operator.
   HalfInt operator-() const {
      return HalfInt(-0.5*this->integer_);
   }
      
   //! @brief Operator overloading: casting operator to ValueType.
   //! @tparam ValueType Value type.
   template<typename ValueType>
   operator ValueType() const noexcept {
      return static_cast<ValueType>(0.5*integer_);
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
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator+(const ValueType lhs, const HalfInt &rhs) {
   return static_cast<ReturnType>(lhs + 0.5*rhs.GetInteger());
}

template<> HalfInt operator+<char, HalfInt>(const char, const HalfInt&);
template<> HalfInt operator+<int, HalfInt>(const int, const HalfInt&);
template<> HalfInt operator+<long, HalfInt>(const long, const HalfInt&);
template<> HalfInt operator+<long long, HalfInt>(const long long, const HalfInt&);
template<> HalfInt operator+<unsigned char, HalfInt>(const unsigned char, const HalfInt&);
template<> HalfInt operator+<unsigned int, HalfInt>(const unsigned int, const HalfInt&);
template<> HalfInt operator+<unsigned long, HalfInt>(const unsigned long, const HalfInt&);
template<> HalfInt operator+<unsigned long long, HalfInt>(const unsigned long long, const HalfInt&);

//! @brief Operator overloading: addition operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator+(const HalfInt &lhs, const ValueType rhs) {
   return static_cast<ReturnType>(0.5*lhs.GetInteger() + rhs);
}

template<> HalfInt operator+<char, HalfInt>(const HalfInt&, const char);
template<> HalfInt operator+<int, HalfInt>(const HalfInt&, const int);
template<> HalfInt operator+<long, HalfInt>(const HalfInt&, const long);
template<> HalfInt operator+<long long, HalfInt>(const HalfInt&, const long long);
template<> HalfInt operator+<unsigned char, HalfInt>(const HalfInt&, const unsigned char);
template<> HalfInt operator+<unsigned int, HalfInt>(const HalfInt&, const unsigned int);
template<> HalfInt operator+<unsigned long, HalfInt>(const HalfInt&, const unsigned long);
template<> HalfInt operator+<unsigned long long, HalfInt>(const HalfInt&, const unsigned long long);

//! @brief Operator overloading: addition operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
HalfInt operator+(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt(0.5*(lhs.GetInteger() + rhs.GetInteger()));
}

//------------------------------------------------------------------
//-----------Operator overloading: Subtraction Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: subtraction operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator-(const ValueType lhs, const HalfInt &rhs) {
   return static_cast<ReturnType>(lhs - 0.5*rhs.GetInteger());
}

template<> HalfInt operator-<char, HalfInt>(const char, const HalfInt&);
template<> HalfInt operator-<int, HalfInt>(const int, const HalfInt&);
template<> HalfInt operator-<long, HalfInt>(const long, const HalfInt&);
template<> HalfInt operator-<long long, HalfInt>(const long long, const HalfInt&);
template<> HalfInt operator-<unsigned char, HalfInt>(const unsigned char, const HalfInt&);
template<> HalfInt operator-<unsigned int, HalfInt>(const unsigned int, const HalfInt&);
template<> HalfInt operator-<unsigned long, HalfInt>(const unsigned long, const HalfInt&);
template<> HalfInt operator-<unsigned long long, HalfInt>(const unsigned long long, const HalfInt&);

//! @brief Operator overloading: subtraction operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator-(const HalfInt &lhs, const ValueType rhs) {
   return static_cast<ReturnType>(0.5*lhs.GetInteger() - rhs);
}

template<> HalfInt operator-<char, HalfInt>(const HalfInt&, const char);
template<> HalfInt operator-<int, HalfInt>(const HalfInt&, const int);
template<> HalfInt operator-<long, HalfInt>(const HalfInt&, const long);
template<> HalfInt operator-<long long, HalfInt>(const HalfInt&, const long long);
template<> HalfInt operator-<unsigned char, HalfInt>(const HalfInt&, const unsigned char);
template<> HalfInt operator-<unsigned int, HalfInt>(const HalfInt&, const unsigned int);
template<> HalfInt operator-<unsigned long, HalfInt>(const HalfInt&, const unsigned long);
template<> HalfInt operator-<unsigned long long, HalfInt>(const HalfInt&, const unsigned long long);

//! @brief Operator overloading: subtraction operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
HalfInt operator-(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt(0.5*(lhs.GetInteger() - rhs.GetInteger()));
}

//------------------------------------------------------------------
//----------Operator overloading: Multiplication Operator-----------
//------------------------------------------------------------------
//! @brief Operator overloading: multiplication operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator*(const ValueType lhs, const HalfInt &rhs) {
   return static_cast<ReturnType>(lhs*0.5*rhs.GetInteger());
}

template<> HalfInt operator*<char, HalfInt>(const char, const HalfInt&);
template<> HalfInt operator*<int, HalfInt>(const int, const HalfInt&);
template<> HalfInt operator*<long, HalfInt>(const long, const HalfInt&);
template<> HalfInt operator*<long long, HalfInt>(const long long, const HalfInt&);
template<> HalfInt operator*<unsigned char, HalfInt>(const unsigned char, const HalfInt&);
template<> HalfInt operator*<unsigned int, HalfInt>(const unsigned int, const HalfInt&);
template<> HalfInt operator*<unsigned long, HalfInt>(const unsigned long, const HalfInt&);
template<> HalfInt operator*<unsigned long long, HalfInt>(const unsigned long long, const HalfInt&);

//! @brief Operator overloading: multiplication operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator*(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger()*rhs;
}

template<> HalfInt operator*<char, HalfInt>(const HalfInt&, const char);
template<> HalfInt operator*<int, HalfInt>(const HalfInt&, const int);
template<> HalfInt operator*<long, HalfInt>(const HalfInt&, const long);
template<> HalfInt operator*<long long, HalfInt>(const HalfInt&, const long long);
template<> HalfInt operator*<unsigned char, HalfInt>(const HalfInt&, const unsigned char);
template<> HalfInt operator*<unsigned int, HalfInt>(const HalfInt&, const unsigned int);
template<> HalfInt operator*<unsigned long, HalfInt>(const HalfInt&, const unsigned long);
template<> HalfInt operator*<unsigned long long, HalfInt>(const HalfInt&, const unsigned long long);

//! @brief Operator overloading: multiplication operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
double operator*(const HalfInt &lhs, const HalfInt &rhs) {
   return 0.5*lhs.GetInteger()*0.5*rhs.GetInteger();
}

//------------------------------------------------------------------
//-------------Operator overloading: Division Operator--------------
//------------------------------------------------------------------
//! @brief Operator overloading: division operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator/(const ValueType lhs, const HalfInt &rhs) {
   return static_cast<ReturnType>(lhs/(0.5*rhs.GetInteger()));
}

template<> double operator/<char, double>(char, const HalfInt&);
template<> double operator/<int, double>(int, const HalfInt&);
template<> double operator/<long, double>(long, const HalfInt&);
template<> double operator/<long long, double>(long long, const HalfInt&);
template<> double operator/<unsigned char, double>(unsigned char, const HalfInt&);
template<> double operator/<unsigned int, double>(unsigned int, const HalfInt&);
template<> double operator/<unsigned long, double>(unsigned long, const HalfInt&);
template<> double operator/<unsigned long long, double>(unsigned long long, const HalfInt&);

//! @brief Operator overloading: division operator.
//! @tparam ValueType Value type.
//! @tparam ReturnType Return type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType, typename ReturnType = ValueType>
ReturnType operator/(const HalfInt &lhs, const ValueType rhs) {
   return static_cast<ReturnType>(0.5*lhs.GetInteger()/rhs);
}

template<> double operator/<char, double>(const HalfInt&, char);
template<> double operator/<int, double>(const HalfInt&, int);
template<> double operator/<long, double>(const HalfInt&, long);
template<> double operator/<long long, double>(const HalfInt&, long long);
template<> double operator/<unsigned char, double>(const HalfInt&, unsigned char);
template<> double operator/<unsigned int, double>(const HalfInt&, unsigned int);
template<> double operator/<unsigned long, double>(const HalfInt&, unsigned long);
template<> double operator/<unsigned long long, double>(const HalfInt&, unsigned long long);

//! @brief Operator overloading: division operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
double operator/(const HalfInt &lhs, const HalfInt &rhs) {
   return static_cast<double>(lhs.GetInteger()/rhs.GetInteger());
}

//------------------------------------------------------------------
//-------------Operator overloading: Equality Operator--------------
//------------------------------------------------------------------
//! @brief Operator overloading: equality operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator==(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() == rhs.GetInteger();
}

//! @brief Operator overloading: equality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator==(const ValueType lhs, const HalfInt &rhs) {
   return lhs == 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: equality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator==(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger() == rhs;
}

//------------------------------------------------------------------
//------------Operator overloading: Inequality Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: inequality operator.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator!=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() != rhs.GetInteger();
}

//! @brief Operator overloading: inequality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator!=(const ValueType lhs, const HalfInt &rhs) {
   return lhs != 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: inequality operator.
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator!=(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger() != rhs;
}

//------------------------------------------------------------------
//------------Operator overloading: Comparison Operator-------------
//------------------------------------------------------------------
//! @brief Operator overloading: comparison operator (less than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator<(const ValueType lhs, const HalfInt &rhs) {
   return lhs < 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (less than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator<(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() < rhs;
}

//! @brief Operator overloading: comparison operator (less than).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator<(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() < rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator>(const ValueType lhs, const HalfInt &rhs) {
   return  lhs > 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator>(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() > rhs;
}

//! @brief Operator overloading: comparison operator (greater than).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator>(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() > rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (less than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator<=(const ValueType lhs, const HalfInt &rhs) {
   return lhs <= 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (less than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator<=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() <= rhs;
}

//! @brief Operator overloading: comparison operator (less than or equal).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator<=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() <= rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator>=(const ValueType lhs, const HalfInt &rhs) {
   return lhs >= 0.5*rhs.GetInteger();
}

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @tparam ValueType Value type.
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
template<typename ValueType>
bool operator>=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() >= rhs;
}

//! @brief Operator overloading: comparison operator (greater than or equal).
//! @param lhs The value of the left-hand side.
//! @param rhs The value of the right-hand side.
bool operator>=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() >= rhs.GetInteger();
}

//------------------------------------------------------------------
//----------------Operator overloading: I/O Stream------------------
//------------------------------------------------------------------
//! @brief Operator overloading: output operator.
//! @param os Ostream object.
//! @param half_int Half-integer.
std::ostream& operator<<(std::ostream &os, const HalfInt &half_int) {
   os << half_int.GetInteger()*0.5;
   return os;
}

} // namespace type
} // namespace compnal

#endif /* COMPNAL_TYPE_HALF_INT_HPP_ */
