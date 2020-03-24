// -*- C++ -*-
//
// PhysicalQtyComplex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2019 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_Complex_H
#define Physical_Qty_Complex_H
#include "PhysicalQty.h"
#include "PhysicalQtyOps.h"
#include <complex>

/** @file PhysicalQtyComplex.h 
 * Overloads for operations on complex physical quantities.
 */

namespace std {
  /**
   *  Template specialization for std::complex<Qty<0,0,0> > 
   *  with conversions to complex<double>
   */
  template<>
  class complex<ThePEG::QtyDouble>
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
    complex<ThePEG::QtyDouble> & 
    operator+=(const complex<ThePEG::QtyDouble> x) { 
      rawValue_ += x.rawValue(); 
      return *this; 
    }
    
    /// Subtraction-assignment
    complex<ThePEG::QtyDouble> & 
    operator-=(const complex<ThePEG::QtyDouble> x) { 
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
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline constexpr auto
operator*(std::complex<Qty<L1,E1,Q1>> q1, 
          std::complex<Qty<L2,E2,Q2>> q2) 
-> std::complex<decltype(q1.real()*q2.real())>
{
  return {q1.real()*q2.real() - q1.imag()*q2.imag(),
          q1.real()*q2.imag() + q1.imag()*q2.real()};
}

// complex qty = complex qty * complex qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<typename Qty<L,E,Q>::Squared>
operator*(std::complex<Qty<L,E,Q>> q1, 
          std::complex<Qty<L,E,Q>> q2) 
{
  return {q1.real()*q2.real() - q1.imag()*q2.imag(),
          q1.real()*q2.imag() + q1.imag()*q2.real()};
}

// complex qty = complex double - complex qty
inline constexpr std::complex<double>
operator-(std::complex<double> q1, std::complex<QtyDouble> q2) {
  return {q1.real()-q2.real(), q1.imag()-q2.imag()};
}

// complex qty = complex double + complex qty
inline constexpr std::complex<double>
operator+(std::complex<double> q1, std::complex<QtyDouble> q2) {
  return {q1.real()+q2.real(), q1.imag()+q2.imag()};
}

// complex qty = complex double * complex qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>>
operator*(std::complex<double> q1, std::complex<Qty<L,E,Q>> q2) {
  return {q1.real()*q2.real() - q1.imag()*q2.imag(),
          q1.real()*q2.imag() + q1.imag()*q2.real()};
}

// complex qty = complex double / complex qty
template<typename L, typename E, typename Q>
inline std::complex<typename Qty<L,E,Q>::Inverse>
operator/(std::complex<double> q1, std::complex<Qty<L,E,Q>> q2) {
  auto tmp  =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex qty = complex double / qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<typename Qty<L,E,Q>::Inverse>
operator/(std::complex<double> q1, Qty<L,E,Q> q2) {
  return {q1.real()/q2, q1.imag()/q2};
}

// complex qty = complex qty / complex double
template<typename L, typename E, typename Q>
inline std::complex<Qty<L,E,Q>>
operator/(std::complex<Qty<L,E,Q>> q1, std::complex<double> q2) {
  auto tmp  =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex qty = qty / complex double
template<typename L, typename E, typename Q>
inline std::complex<Qty<L,E,Q>>
operator/(Qty<L,E,Q> q1, std::complex<double> q2) {
  auto tmp  =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex double = complex qty / complex qty
template<typename L, typename E, typename Q>
inline std::complex<double>
operator/(std::complex<Qty<L,E,Q>> q1, 
          std::complex<Qty<L,E,Q>> q2) {
  auto tmp  =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex double = qty / complex qty
template<typename L, typename E, typename Q>
inline std::complex<double>
operator/(Qty<L,E,Q> q1, std::complex<Qty<L,E,Q>> q2) {
  auto tmp = q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex double = complex qty / qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<double>
operator/(std::complex<Qty<L,E,Q>> q1, Qty<L,E,Q> q2) {
  return {q1.real()/q2, q1.imag()/q2};
}

// complex qty = complex qty / complex qty
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline auto
operator/(std::complex<Qty<L1,E1,Q1>> q1, 
          std::complex<Qty<L2,E2,Q2>> q2) 
-> std::complex<decltype(q1.real()/q2.real())>
{
  auto  tmp =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex qty = qty / complex qty
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline auto 
operator/(Qty<L1,E1,Q1> q1, 
          std::complex<Qty<L2,E2,Q2>> q2) 
-> std::complex<decltype(q1/q2.real())>
{
  auto  tmp =  q1*conj(q2);
  auto norm = (q2*conj(q2)).real();
  return {tmp.real()/norm, tmp.imag()/norm};
}

// complex qty = complex qty / qty
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline constexpr auto
operator/(std::complex<Qty<L1,E1,Q1>> q1, Qty<L2,E2,Q2> q2) 
-> std::complex<decltype(q1.real()/q2)>
{
  return {q1.real()/q2, q1.imag()/q2};
}


// complex qty = complex qty * complex double
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>>
operator*(std::complex<Qty<L,E,Q>> q1, std::complex<double> q2) {
  return q2 * q1;
}


// complex qty = qty * complex qty
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline constexpr auto
operator*(Qty<L1,E1,Q1> q1, std::complex<Qty<L2,E2,Q2>> q2) 
-> std::complex<decltype(q1*q2.real())>
{
  return {q1*q2.real(), q1*q2.imag()};
}

// complex qty = qty * complex qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<typename Qty<L,E,Q>::Squared>
operator*(Qty<L,E,Q> q1, std::complex<Qty<L,E,Q>> q2) {
  return {q1*q2.real(), q1*q2.imag()};
}

// complex qty = qty * complex double
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>>
operator*(Qty<L,E,Q> q1, std::complex<double> q2) {
  return {q1*q2.real(), q1*q2.imag()};
}

// complex qty = complex double * qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>>
operator*(std::complex<double> q1, Qty<L,E,Q> q2) {
  return q2 * q1;
}


// complex qty = complex qty * qty
template<typename L1, typename E1, typename Q1,
         typename L2, typename E2, typename Q2>
inline constexpr auto
operator*(std::complex<Qty<L1,E1,Q1>> q1, Qty<L2,E2,Q2> q2) 
-> decltype(q2*q1)
{
  return q2 * q1;
}

// complex qty = complex qty * qty
template<typename L, typename E, typename Q>
inline constexpr std::complex<typename Qty<L,E,Q>::Squared>
operator*(std::complex<Qty<L,E,Q>> q1, Qty<L,E,Q> q2) {
  return q2 * q1;
}

// complex qty *= double
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>> &
operator*=(std::complex<Qty<L,E,Q>> & q1, double q2) {
  return (q1 = q1 * q2);
}

// complex qty /= double
template<typename L, typename E, typename Q>
inline constexpr std::complex<Qty<L,E,Q>> &
operator/=(std::complex<Qty<L,E,Q>> & q1, double q2) {
  return (q1 = q1 / q2);
}
//@}
}

#endif
