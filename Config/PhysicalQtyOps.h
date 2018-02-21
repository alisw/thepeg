// -*- C++ -*-
//
// PhysicalQtyOps.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2017 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_Ops_H
#define Physical_Qty_Ops_H
#include <cmath>

/** @file PhysicalQtyOps.h 
 * Overloads for mathematical operations on physical quantities.
 */

namespace ThePEG {
/// @name Overloads for mathematical operations on physical quantities.
//@{
// qty = qty * qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> 
operator*(Qty<L1,E1,Q1,DL1,DE1,DQ1> q1, Qty<L2,E2,Q2,DL2,DE2,DQ2> q2) {
  typedef
    Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> RetT;
  return RetT{RetT::baseunit(), q1.rawValue()*q2.rawValue()};
}


// qty = qty / qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> 
operator/(Qty<L1,E1,Q1,DL1,DE1,DQ1> q1, Qty<L2,E2,Q2,DL2,DE2,DQ2> q2) {
  typedef
    Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> RetT;
  return RetT{RetT::baseunit(), q1.rawValue()/q2.rawValue()};
}

// qty = qty + qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline Qty<L,E,Q,DL,DE,DQ> 
operator+(Qty<L,E,Q,DL,DE,DQ> q1, 
	  Qty<QtyInt<L,DL,DL2>::I,
	      QtyInt<E,DE,DE2>::I,
	      QtyInt<Q,DQ,DQ2>::I, 
	  DL2,DE2,DQ2> q2) {
  Qty<L,E,Q,DL,DE,DQ> q = q1;
  q += q2;
  return q;
}

// qty = qty - qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline Qty<L,E,Q,DL,DE,DQ> 
operator-(Qty<L,E,Q,DL,DE,DQ> q1, 
	  Qty<QtyInt<L,DL,DL2>::I,
	      QtyInt<E,DE,DE2>::I,
	      QtyInt<Q,DQ,DQ2>::I, 
	  DL2,DE2,DQ2> q2) {
  Qty<L,E,Q,DL,DE,DQ> q = q1;
  q -= q2;
  return q;
}

// qty == qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator==(Qty<L,E,Q,DL,DE,DQ> q1,
	   Qty<QtyInt<L,DL,DL2>::I,
	       QtyInt<E,DE,DE2>::I,
	       QtyInt<Q,DQ,DQ2>::I, 
	       DL2,DE2,DQ2> q2) {
  return q1.rawValue()==q2.rawValue();
}

// qty != qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator!=(Qty<L,E,Q,DL,DE,DQ> q1,
	   Qty<QtyInt<L,DL,DL2>::I,
	       QtyInt<E,DE,DE2>::I,
	       QtyInt<Q,DQ,DQ2>::I, 
	       DL2,DE2,DQ2> q2) {
  return q1.rawValue()!=q2.rawValue();
}

// qty < qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator<(Qty<L,E,Q,DL,DE,DQ> q1,
          Qty<QtyInt<L,DL,DL2>::I,
              QtyInt<E,DE,DE2>::I,
              QtyInt<Q,DQ,DQ2>::I, 
              DL2,DE2,DQ2> q2) {
  return q1.rawValue()<q2.rawValue();
}

// qty <= qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator<=(Qty<L,E,Q,DL,DE,DQ> q1,
	   Qty<QtyInt<L,DL,DL2>::I,
	       QtyInt<E,DE,DE2>::I,
	       QtyInt<Q,DQ,DQ2>::I, 
	       DL2,DE2,DQ2> q2) {
  return q1.rawValue()<=q2.rawValue();
}

// qty > qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator>(Qty<L,E,Q,DL,DE,DQ> q1,
	  Qty<QtyInt<L,DL,DL2>::I,
	      QtyInt<E,DE,DE2>::I,
              QtyInt<Q,DQ,DQ2>::I, 
              DL2,DE2,DQ2> q2) {
  return q1.rawValue()>q2.rawValue();
}

// qty >= qty
template<int L, int E, int Q, int DL, int DE, int DQ, int DL2, int DE2, int DQ2>
inline constexpr bool
operator>=(Qty<L,E,Q,DL,DE,DQ> q1,
	   Qty<QtyInt<L,DL,DL2>::I,
	       QtyInt<E,DE,DE2>::I,
               QtyInt<Q,DQ,DQ2>::I, 
               DL2,DE2,DQ2> q2) {
  return q1.rawValue()>=q2.rawValue();
}

// comparisons with ZERO
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator==(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() == 0.0;
}
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator!=(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() != 0.0;
}
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator<(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() < 0.0;
}
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator>(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() > 0.0;
}
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator<=(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() <= 0.0;
}
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr bool
operator>=(Qty<L,E,Q,DL,DE,DQ> q1, ZeroUnit) {
  return q1.rawValue() >= 0.0;
}

// qty = qty * double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<L,E,Q,DL,DE,DQ>
operator*(Qty<L,E,Q,DL,DE,DQ> q,double x) {
  return Qty<L,E,Q,DL,DE,DQ>{q,x};
}

// qty = double * qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<L,E,Q,DL,DE,DQ>
operator*(double x,Qty<L,E,Q,DL,DE,DQ> q) {
  return Qty<L,E,Q,DL,DE,DQ>{q,x};
}

// qty = qty / double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<L,E,Q,DL,DE,DQ>
operator/(Qty<L,E,Q,DL,DE,DQ> q,double x) {
  return Qty<L,E,Q,DL,DE,DQ>{q, 1./x};
}

// qty = double / qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<-L,-E,-Q,DL,DE,DQ>
operator/(double x, Qty<L,E,Q,DL,DE,DQ> q) {
  typedef Qty<-L,-E,-Q,DL,DE,DQ> RetT;
  return RetT{RetT::baseunit(), x/q.rawValue()};
}

// qty = -qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<L,E,Q,DL,DE,DQ>
operator-(Qty<L,E,Q,DL,DE,DQ> q) {
  typedef Qty<L,E,Q,DL,DE,DQ> RetT;
  return RetT{RetT::baseunit(), -q.rawValue()};
}

// qty = sqr(qty)
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<2*L,2*E,2*Q,DL,DE,DQ>
sqr(Qty<L,E,Q,DL,DE,DQ> q) {
  return q*q;
}

// qty = sqrt(qty) // std::sqrt is not constexpr
template<int L, int E, int Q, int DL, int DE, int DQ>
inline Qty<L,E,Q,DL*2,DE*2,DQ*2>
sqrt(Qty<L,E,Q,DL,DE,DQ> q) {
  typedef Qty<L,E,Q,DL*2,DE*2,DQ*2> RetT;
  return RetT(std::sqrt(q.rawValue())*RetT::baseunit());
}

// double = atan2(y,x)
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr double
atan2(Qty<L,E,Q,DL,DE,DQ> y, Qty<L,E,Q,DL,DE,DQ> x) {
  return std::atan2(y.rawValue(), x.rawValue());
}

// qty = abs(qty)
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<L,E,Q,DL,DE,DQ>
abs(Qty<L,E,Q,DL,DE,DQ> q) {
  typedef Qty<L,E,Q,DL,DE,DQ> RetT;
  return RetT(std::abs(q.rawValue())*RetT::baseunit());
}

// qty = pow<P,R>(qty)
template<int P, int R, int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr Qty<P*L,P*E,P*Q,R*DL,R*DE,R*DQ> 
pow(Qty<L,E,Q,DL,DE,DQ> q) {
  typedef Qty<P*L,P*E,P*Q,R*DL,R*DE,R*DQ> RetT;
  return RetT(std::pow(q.rawValue(),double(P)/double(R))*RetT::baseunit());
}
  
// max for T,U types
template<typename T, typename U>
inline
T max(const T & t, const U & u) {
  const T & utmp = u;
  return std::max(t, utmp);
}

// ZeroUnit in front should take U type
template<typename U>
inline
U max(const ZeroUnit & t, const U & u) {
  const U & ttmp = t;
  return std::max(ttmp, u);
}

// min for T,U types
template<typename T, typename U>
inline
T min(const T & t, const U & u) {
  const T & utmp = u;
  return std::min(t, utmp);
}

// ZeroUnit in front should take U type
template<typename U>
inline
U min(const ZeroUnit & t, const U & u) {
  const U & ttmp = t;
  return std::min(ttmp, u);
}


//@}
}

#endif
