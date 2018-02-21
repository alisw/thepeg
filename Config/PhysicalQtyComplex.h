// -*- C++ -*-
//
// PhysicalQtyComplex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2017 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_Complex_H
#define Physical_Qty_Complex_H
#include <complex>

/** @file PhysicalQtyComplex.h 
 * Overloads for operations on complex physical quantities.
 */

namespace std {
  /**
   *  Template specialization for std::complex<Qty<0,0,0> > 
   *  with conversions to complex<double>
   */
  template<int DL, int DE, int DQ>
  class complex<ThePEG::Qty<0,0,0,DL,DE,DQ> >
  {
  public:
    /// Default constructor
    constexpr complex(double r=0.0, double i=0.0) 
      : rawValue_(r,i) {}

    /// Constructor from complex<double>
    constexpr complex(complex<double> C)
      : rawValue_(C) {}

    /**
     * The internal representation of the dimensionful quantity.
     * Using this will break dimension-consistency.
     */ 
    constexpr complex<double> rawValue() const { return rawValue_; }

    /// Real part
    constexpr double real() const { return rawValue_.real(); }

    /// Imaginary part
    constexpr double imag() const { return rawValue_.imag(); }
   
    /// Cast to complex<double>
    constexpr operator complex<double>() const {
      return rawValue_;
    }
    
    /// Addition-assignment
    complex<ThePEG::Qty<0,0,0,DL,DE,DQ> > & 
    operator+=(const complex<ThePEG::Qty<0,0,0,DL,DE,DQ> > x) { 
      rawValue_ += x.rawValue(); 
      return *this; 
    }
    
    /// Subtraction-assignment
    complex<ThePEG::Qty<0,0,0,DL,DE,DQ> > & 
    operator-=(const complex<ThePEG::Qty<0,0,0,DL,DE,DQ> > x) { 
      rawValue_ -= x.rawValue(); 
      return *this; 
    }

  private:
    /// Internal value of the dimensioned quantity
    complex<double> rawValue_;
  };
}
// =========================================

namespace ThePEG {

/// @name Overloads for mathematical operations
//@{
// complex qty = complex qty * complex qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
			DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator*(std::complex<Qty<L1,E1,Q1,DL1,DE1,DQ1> > q1, 
	  std::complex<Qty<L2,E2,Q2,DL2,DE2,DQ2> > q2) {
  typedef std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
    DL1*DL2,DE1*DE2,DQ1*DQ2> > RetT;
  return RetT(q1.real()*q2.real() - q1.imag()*q2.imag(),
	      q1.real()*q2.imag() + q1.imag()*q2.real());
}

// complex qty = complex qty * complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> >
operator*(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  typedef std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> > RetT;
  return RetT(q1.real()*q2.real() - q1.imag()*q2.imag(),
	      q1.real()*q2.imag() + q1.imag()*q2.real());
}

// complex qty = complex double - complex qty
template<int DL, int DE, int DQ>
inline constexpr std::complex<double>
operator-(std::complex<double> q1, 
	  std::complex<Qty<0,0,0,DL,DE,DQ> > q2) {
  typedef std::complex<double> RetT;
  return RetT(q1.real()-q2.real(),q1.imag()-q2.imag());
}

// complex qty = complex double * complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator*(std::complex<double> q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  typedef std::complex<Qty<L,E,Q,DL,DE,DQ> > RetT;
  return RetT(q1.real()*q2.real() - q1.imag()*q2.imag(),
	      q1.real()*q2.imag() + q1.imag()*q2.real());
}

// complex qty = complex double / complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<Qty<-L,-E,-Q,DL,DE,DQ> >
operator/(std::complex<double> q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  typedef std::complex<Qty<-L,-E,-Q,DL,DE,DQ> > RetT;
  std::complex<Qty<L,E,Q,DL,DE,DQ> > tmp = q1*conj(q2);
  Qty<2*L,2*E,2*Q,DL,DE,DQ> norm = (q2*conj(q2)).real();
  return RetT(tmp.real()/norm,tmp.imag()/norm);
}

// complex qty = complex double / qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<-L,-E,-Q,DL,DE,DQ> >
operator/(std::complex<double> q1, 
	  Qty<L,E,Q,DL,DE,DQ>  q2) {
  typedef std::complex<Qty<-L,-E,-Q,DL,DE,DQ> > RetT;
  return RetT(q1.real()/q2,q1.imag()/q2);
}

// complex qty = complex qty / complex double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator/(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  std::complex<double> q2) {
  std::complex<Qty<L,E,Q,DL,DE,DQ> > tmp = q1*conj(q2);
  double norm = (q2*conj(q2)).real();
  return std::complex<Qty<L,E,Q,DL,DE,DQ> >(tmp.real()/norm,tmp.imag()/norm);
}

// complex qty = qty / complex double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator/(Qty<L,E,Q,DL,DE,DQ> q1, 
	  std::complex<double> q2) {
  std::complex<Qty<L,E,Q,DL,DE,DQ> > tmp = q1*conj(q2);
  double norm = (q2*conj(q2)).real();
  return std::complex<Qty<L,E,Q,DL,DE,DQ> >(tmp.real()/norm,tmp.imag()/norm);
}

// complex double = complex qty / complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<double>
operator/(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> > tmp = q1*conj(q2);
  Qty<2*L,2*E,2*Q,DL,DE,DQ> norm = (q2*conj(q2)).real();
  return std::complex<double>(tmp.real()/norm,tmp.imag()/norm);
}

// complex double = qty / complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<double>
operator/(Qty<L,E,Q,DL,DE,DQ> q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> > tmp = q1*conj(q2);
  Qty<2*L,2*E,2*Q,DL,DE,DQ> norm = (q2*conj(q2)).real();
  return std::complex<double>(tmp.real()/norm,tmp.imag()/norm);
}

// complex double = complex qty / qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<double>
operator/(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  Qty<L,E,Q,DL,DE,DQ> q2) {
  return std::complex<double>(q1.real()/q2,q1.imag()/q2);
}

// complex qty = complex qty / complex qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator/(std::complex<Qty<L1,E1,Q1,DL1,DE1,DQ1> > q1, 
	  std::complex<Qty<L2,E2,Q2,DL2,DE2,DQ2> > q2) {
  typedef std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,
                           DL1*DL2,DE1*DE2,DQ1*DQ2> > RetT;
  std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
    DL1*DL2,DE1*DE2,DQ1*DQ2> > tmp = q1*conj(q2);
  Qty<2*L2,2*E2,2*Q2,DL2,DE2,DQ2> norm = (q2*conj(q2)).real();
  return RetT(tmp.real()/norm,tmp.imag()/norm);
}

// complex qty = qty / complex qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator/(Qty<L1,E1,Q1,DL1,DE1,DQ1> q1, 
	  std::complex<Qty<L2,E2,Q2,DL2,DE2,DQ2> > q2) {
  typedef std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,
                           DL1*DL2,DE1*DE2,DQ1*DQ2> > RetT;
  std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
    DL1*DL2,DE1*DE2,DQ1*DQ2> > tmp = q1*conj(q2);
  Qty<2*L2,2*E2,2*Q2,DL2,DE2,DQ2> norm = (q2*conj(q2)).real();
  return RetT(tmp.real()/norm,tmp.imag()/norm);
}

// complex qty = complex qty / qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator/(std::complex<Qty<L1,E1,Q1,DL1,DE1,DQ1> > q1, 
	  Qty<L2,E2,Q2,DL2,DE2,DQ2> q2) {
  typedef std::complex<Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,
                           DL1*DL2,DE1*DE2,DQ1*DQ2> > RetT;
  return RetT(q1.real()/q2,q1.imag()/q2);
}


// complex qty = complex qty * complex double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator*(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  std::complex<double> q2) {
  return q2 * q1;
}


// complex qty = qty * complex qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
			DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator*(Qty<L1,E1,Q1,DL1,DE1,DQ1> q1, 
	  std::complex<Qty<L2,E2,Q2,DL2,DE2,DQ2> > q2) {
  typedef std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
    DL1*DL2,DE1*DE2,DQ1*DQ2> > RetT;
  return RetT(q1*q2.real(), q1*q2.imag());
}

// complex qty = qty * complex qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> >
operator*(Qty<L,E,Q,DL,DE,DQ> q1, 
	  std::complex<Qty<L,E,Q,DL,DE,DQ> > q2) {
  typedef std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> > RetT;
  return RetT(q1*q2.real(), q1*q2.imag());
}

// complex qty = qty * complex double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator*(Qty<L,E,Q,DL,DE,DQ> q1, 
	  std::complex<double> q2) {
  typedef std::complex<Qty<L,E,Q,DL,DE,DQ> > RetT;
  return RetT(q1*q2.real(), q1*q2.imag());
}

// complex qty = complex double * qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<L,E,Q,DL,DE,DQ> >
operator*(std::complex<double> q1,
	  Qty<L,E,Q,DL,DE,DQ> q2) {
  return q2 * q1;
}


// complex qty = complex qty * qty
template<int L1, int L2, int E1, int E2, int Q1, int Q2,
  int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
inline constexpr std::complex<Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
			DL1*DL2,DE1*DE2,DQ1*DQ2> >
operator*(std::complex<Qty<L1,E1,Q1,DL1,DE1,DQ1> > q1, 
	  Qty<L2,E2,Q2,DL2,DE2,DQ2> q2) {
  return q2 * q1;
}

// complex qty = complex qty * qty
template<int L, int E, int Q, int DL, int DE, int DQ>
inline constexpr std::complex<Qty<2*L,2*E,2*Q,DL,DE,DQ> >
operator*(std::complex<Qty<L,E,Q,DL,DE,DQ> > q1, 
	  Qty<L,E,Q,DL,DE,DQ> q2) {
  return q2 * q1;
}

// // complex qty *= complex double
// template<int L, int E, int Q, int DL, int DE, int DQ>
// inline std::complex<Qty<L,E,Q,DL,DE,DQ> > &
// operator*=(std::complex<Qty<L,E,Q,DL,DE,DQ> > & q1,
// 	   std::complex<double> q2) {
//   q1 = q1 * q2;
//   return q1;
// }

// complex qty *= double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<Qty<L,E,Q,DL,DE,DQ> > &
operator*=(std::complex<Qty<L,E,Q,DL,DE,DQ> > & q1,
	   double q2) {
  q1 = q1 * q2;
  return q1;
}

// // complex qty /= complex double
// template<int L, int E, int Q, int DL, int DE, int DQ>
// inline std::complex<Qty<L,E,Q,DL,DE,DQ> > &
// operator/=(std::complex<Qty<L,E,Q,DL,DE,DQ> > & q1,
// 	   std::complex<double> q2) {
//   q1 = q1 / q2;
//   return q1;
// }

// complex qty /= double
template<int L, int E, int Q, int DL, int DE, int DQ>
inline std::complex<Qty<L,E,Q,DL,DE,DQ> > &
operator/=(std::complex<Qty<L,E,Q,DL,DE,DQ> > & q1,
	   double q2) {
  q1 = q1 / q2;
  return q1;
}
//@}
}

#endif
