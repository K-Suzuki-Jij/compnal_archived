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

#ifndef COMPNAL_UTILITY_TYPE_HPP_
#define COMPNAL_UTILITY_TYPE_HPP_

#include <iostream>

namespace compnal {

class HalfInt {
   
public:
   HalfInt() {};
   
   template<typename ValueType>
   HalfInt(ValueType value) {
      value = 2*value;
      if (std::floor(value) != value) {
         throw std::runtime_error("The input number is not half-integer");
      }
      integer_ = static_cast<int>(value);
   };
   
   int GetInteger() const {
      return integer_;
   }
   
   HalfInt operator-() const {
      return HalfInt(-0.5*this->integer_);
   }
   
   HalfInt operator+() const {
      return *this;
   }
   
   operator int   () const noexcept { return static_cast<int>(0.5*integer_); }
   operator double() const noexcept { return 0.5*integer_; }

private:
   int integer_ = 0;
      
};

// Operator overloading: addition
double operator+(const double lhs, const HalfInt &rhs) {
   return lhs + 0.5*rhs.GetInteger();
}

double operator+(const HalfInt &lhs, const double rhs) {
   return 0.5*lhs.GetInteger() + rhs;
}

HalfInt operator+(const int lhs, const HalfInt &rhs) {
   return HalfInt(lhs + 0.5*rhs.GetInteger());
}

HalfInt operator+(const HalfInt &lhs, const int rhs) {
   return HalfInt(0.5*lhs.GetInteger() + rhs);
}

HalfInt operator+(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt(0.5*(lhs.GetInteger() + rhs.GetInteger()));
}


// Operator overloading: subtraction
double operator-(const double lhs, const HalfInt &rhs) {
   return lhs - 0.5*rhs.GetInteger();
}

double operator-(const HalfInt &lhs, const double rhs) {
   return 0.5*lhs.GetInteger() - rhs;
}

HalfInt operator-(const int lhs, const HalfInt &rhs) {
   return HalfInt(lhs - 0.5*rhs.GetInteger());
}

HalfInt operator-(const HalfInt &lhs, const int rhs) {
   return HalfInt(0.5*lhs.GetInteger() - rhs);
}

HalfInt operator-(const HalfInt &lhs, const HalfInt &rhs) {
   return HalfInt(0.5*(lhs.GetInteger() - rhs.GetInteger()));
}

// Operator overloading: multiplication
double operator*(const double lhs, const HalfInt &rhs) {
   return lhs*0.5*rhs.GetInteger();
}

double operator*(const HalfInt &lhs, const double rhs) {
   return 0.5*lhs.GetInteger()*rhs;
}

double operator*(const HalfInt &lhs, const HalfInt &rhs) {
   return 0.5*lhs.GetInteger()*0.5*rhs.GetInteger();
}

HalfInt operator*(const int lhs, const HalfInt &rhs) {
   return HalfInt(lhs*0.5*rhs.GetInteger());
}

HalfInt operator*(const HalfInt &lhs, const int rhs) {
   return HalfInt(0.5*lhs.GetInteger()*rhs);
}

// Operator overloading: division
template<typename ValueType>
double operator/(const ValueType lhs, const HalfInt &rhs) {
   return lhs/(0.5*rhs.GetInteger());
}

template<typename ValueType>
double operator/(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger()/rhs;
}

double operator/(const HalfInt &lhs, const HalfInt &rhs) {
   return static_cast<double>(lhs.GetInteger()/rhs.GetInteger());
}

// Operator overloading: Equality Compare
bool operator==(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() == rhs.GetInteger();
}

template<typename ValueType>
bool operator==(const ValueType lhs, const HalfInt &rhs) {
   return lhs == 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator==(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger() == rhs;
}

bool operator!=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() != rhs.GetInteger();
}

template<typename ValueType>
bool operator!=(const ValueType lhs, const HalfInt &rhs) {
   return lhs != 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator!=(const HalfInt &lhs, const ValueType rhs) {
   return 0.5*lhs.GetInteger() != rhs;
}

// Operator overloading: Compare
template<typename ValueType>
bool operator<(const ValueType lhs, const HalfInt &rhs) {
   return lhs < 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator<(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() < rhs;
}

bool operator<(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() < rhs.GetInteger();
}

template<typename ValueType>
bool operator>(const ValueType lhs, const HalfInt &rhs) {
   return  lhs > 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator>(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() > rhs;
}

bool operator>(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() > rhs.GetInteger();
}

template<typename ValueType>
bool operator<=(const ValueType lhs, const HalfInt &rhs) {
   return lhs <= 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator<=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() <= rhs;
}

bool operator<=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() <= rhs.GetInteger();
}

template<typename ValueType>
bool operator>=(const ValueType lhs, const HalfInt &rhs) {
   return lhs >= 0.5*rhs.GetInteger();
}

template<typename ValueType>
bool operator>=(const HalfInt &lhs, const ValueType &rhs) {
   return 0.5*lhs.GetInteger() >= rhs;
}

bool operator>=(const HalfInt &lhs, const HalfInt &rhs) {
   return lhs.GetInteger() >= rhs.GetInteger();
}

//User-defined literals
HalfInt operator"" _hi(const char* value) {
   return HalfInt(std::stold(value, nullptr));
}

std::ostream& operator<<(std::ostream &os, const HalfInt &half_int) {
   os << half_int.GetInteger()*0.5;
   return os;
}

} // namespace compnal

#endif /* COMPNAL_UTILITY_TYPE_HPP_ */
