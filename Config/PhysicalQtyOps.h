// -*- C++ -*-
//
// PhysicalQtyOps.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2019 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_Ops_H
#define Physical_Qty_Ops_H
#include "PhysicalQty.h"
#include <cmath>

/** @file PhysicalQtyOps.h 
 * Overloads for mathematical operations on physical quantities.
 */

namespace ThePEG {
/// @name Overloads for mathematical operations on physical quantities.
//@{
// qty = qty * qty
template<typename T, typename U>
inline constexpr typename BinaryOpTraits<T,U>::MulT
operator*(T q1, U q2) {
  typedef typename BinaryOpTraits<T,U>::MulT RetT;
  return RetT{RetT::baseunit(), q1.rawValue()*q2.rawValue()};
}

// qty = qty / qty
template<typename T, typename U>
inline constexpr typename BinaryOpTraits<T,U>::DivT
operator/(T q1, U q2) {
  typedef typename BinaryOpTraits<T,U>::DivT RetT;
  return RetT{RetT::baseunit(), q1.rawValue()/q2.rawValue()};
}

// qty = qty + qty
template<typename T, typename U>
inline enable_if_same_qty<T,T,U>
operator+(T q1, U q2) {
  q1 += q2;
  return q1;
}

// qty = qty - qty
template<typename T, typename U>
inline enable_if_same_qty<T,T,U>
operator-(T q1, U q2) {
  q1 -= q2;
  return q1;
}

// qty == qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator==(T q1, U q2) {
  return q1.rawValue()==q2.rawValue();
}

// qty != qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator!=(T q1, U q2) {
  return q1.rawValue()!=q2.rawValue();
}

// qty < qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator<(T q1, U q2) {
  return q1.rawValue()<q2.rawValue();
}

// qty <= qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator<=(T q1, U q2) {
  return q1.rawValue()<=q2.rawValue();
}

// qty > qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator>(T q1, U q2) {
  return q1.rawValue()>q2.rawValue();
}

// qty >= qty
template<typename T, typename U>
inline constexpr enable_if_same_qty<bool,T,U>
operator>=(T q1, U q2) {
  return q1.rawValue()>=q2.rawValue();
}

// comparisons with ZERO
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator==(T q1, ZeroUnit) {
  return q1.rawValue() == 0.0;
}
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator!=(T q1, ZeroUnit) {
  return q1.rawValue() != 0.0;
}
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator<(T q1, ZeroUnit) {
  return q1.rawValue() < 0.0;
}
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator>(T q1, ZeroUnit) {
  return q1.rawValue() > 0.0;
}
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator<=(T q1, ZeroUnit) {
  return q1.rawValue() <= 0.0;
}
template<typename T>
inline constexpr enable_if_same_qty<bool, T>
operator>=(T q1, ZeroUnit) {
  return q1.rawValue() >= 0.0;
}

// qty = qty * double
template<typename T>
inline constexpr enable_if_same_qty<T, T>
operator*(T q,double x) {
  return T{q,x};
}

// qty = double * qty
template<typename T>
inline constexpr enable_if_same_qty<T, T>
operator*(double x,T q) {
  return T{q,x};
}

// qty = qty / double
template<typename T>
inline constexpr enable_if_same_qty<T, T>
operator/(T q,double x) {
  return T{q, 1./x};
}

// qty = double / qty
template<typename T>
inline constexpr enable_if_same_qty<typename T::Inverse, T>
operator/(double x, T q) {
  typedef typename T::Inverse RetT;
  return RetT{RetT::baseunit(), x/q.rawValue()};
}

// qty = -qty
template<typename T>
inline constexpr enable_if_same_qty<T, T>
operator-(T q) {
  typedef T RetT;
  return RetT{q, -1.0};
}

// qty = sqrt(qty) // std::sqrt is not constexpr
template<typename T>
inline enable_if_same_qty<typename T::Sqrt, T>
sqrt(T q) {
  typedef typename T::Sqrt RetT;
  return RetT{RetT::baseunit(), std::sqrt(q.rawValue())};
}

// double = atan2(y,x)
template<typename T, typename U>
inline constexpr enable_if_same_qty<double,T,U>
atan2(T y, U x) {
  return std::atan2(y.rawValue(), x.rawValue());
}

// qty = abs(qty)
template<typename T>
inline constexpr enable_if_same_qty<T, T>
abs(T q) {
  return T{T::baseunit(), std::abs(q.rawValue())};
}

// qty = pow<P,R>(qty)
template<long int Num, long int Den, typename T>
inline constexpr enable_if_same_qty<typename T::template Power<Num,Den>, T>
pow(T q) {
  typedef typename T::template Power<Num,Den> RetT;
  return RetT{RetT::baseunit(), std::pow(q.rawValue(),double(Num)/double(Den))};
}
  
// max for T,U types
template<typename T, typename U>
inline T max(const T & t, const U & u) {
  const T & utmp = u;
  return std::max(t, utmp);
}

// ZeroUnit in front should take U type
template<typename U>
inline U max(const ZeroUnit & t, const U & u) {
  const U & ttmp = t;
  return std::max(ttmp, u);
}

// min for T,U types
template<typename T, typename U>
inline T min(const T & t, const U & u) {
  const T & utmp = u;
  return std::min(t, utmp);
}

// ZeroUnit in front should take U type
template<typename U>
inline U min(const ZeroUnit & t, const U & u) {
  const U & ttmp = t;
  return std::min(ttmp, u);
}


//@}
}

#endif
