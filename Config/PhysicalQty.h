// -*- C++ -*-
//
// PhysicalQty.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2011 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_H
#define Physical_Qty_H
#include "TemplateTools.h"
#include <sstream>

/** @file 
 *
 * The PhysicalQty class allows compile-time checking of dimensional
 * correctness. Mathematical operations that are inconsistent are
 * flagged as type errors.
 *
 * Do not use the classes directly in ThePEG, use the wrappers defined
 * in Units.h or Phys_Qty.h instead.
 */

namespace ThePEG {

/// Helper class to construct zero unitful quantities.
struct ZeroUnit {
  /** Automatic conversion to double. */
  operator double() const { return 0.0; }
};

/// ZERO can be used as zero for any unitful quantity.
const ZeroUnit ZERO = ZeroUnit();

/// Helper classes to extend or shorten fractions
//@{
/**
 * Template to help with fractional powers of dimensions 
 */
template <int M, int II>
struct QtyHelper 
{
  /// The numerator, indicating failure.
  static const int I = -999999;
};

/**
 * Template to help with fractional powers of dimensions
 */
template <int II>
struct QtyHelper<0,II> 
{
  /// The new numerator.
  static const int I = II;
};

/**
 * Template to help with fractional powers of dimensions 
 */
template <int II, int DI, int DI2>
struct QtyInt 
{
  /// The new numerator.
  static const int I = QtyHelper<(DI2*II)%DI,(DI2*II)/DI>::I;
};
//@}

/**
 * This template class allows the compiler to check calculations with
 * physical quantities for dimensional correctness. A quantity can be
 * composed of arbitrary fractional powers of length L, energy E and
 * charge Q. Commonly used quantities should be typedef'ed (see Units.h).
 *
 * Some member functions can break dimensional consistency if the user
 * is not careful; these are marked explicitly.
 *
 * Do not use this class directly in ThePEG, use the pre-defined quantities
 * from Units.h or the wrapper in Phys_Qty.h instead.
 */
template<int L, int E, int Q, int DL = 1, int DE = 1, int DQ = 1>
class Qty
{
private:
  /// Constructor from raw values. Breaks consistency.
  Qty(double val) : rawValue_(val) {}

public:

  /// The name of the class for persistent IO
  static std::string className() {
    std::ostringstream os;
    os << "Qty<" 
       <<  L << ','
       <<  E << ','
       <<  Q << ','
       << DL << ','
       << DE << ','
       << DQ << '>';
      
    return os.str();
  }


  /// The squared type.
  typedef Qty<2*L,2*E,2*Q,DL,DE,DQ> Squared;

  /// Basic unit of this quantity.
  static Qty<L,E,Q,DL,DE,DQ> baseunit() 
  {
    return Qty<L,E,Q,DL,DE,DQ>(1.0);
  }

  /// Default constructor to 0.
  Qty() : rawValue_(0.0) {}

  /// Default constructor to 0.
  Qty(ZeroUnit) : rawValue_(0.0) {}

  /// Constructor from a compatible quantity
  template <int DL2, int DE2, int DQ2>
  Qty(const Qty<QtyInt<L,DL,DL2>::I,
                QtyInt<E,DE,DE2>::I,
                QtyInt<Q,DQ,DQ2>::I,
                DL2,DE2,DQ2> & q)
    : rawValue_(q.rawValue()) {}

  /// Access to the raw value. Breaks consistency.
  double rawValue() const { return rawValue_; }

  /// Assignment multiplication by dimensionless number.
  Qty<L,E,Q,DL,DE,DQ> & operator*=(double x) { rawValue_ *= x; return *this; }

  /// Assignment division by dimensionless number.
  Qty<L,E,Q,DL,DE,DQ> & operator/=(double x) { rawValue_ /= x; return *this; }

  /// Assignment addition with compatible quantity.
  template <int DL2, int DE2, int DQ2>
  Qty<L,E,Q,DL,DE,DQ> & 
  operator+=(const Qty<QtyInt<L,DL,DL2>::I,
	               QtyInt<E,DE,DE2>::I,
	               QtyInt<Q,DQ,DQ2>::I,
	               DL2,DE2,DQ2> x) 
  { 
    rawValue_ += x.rawValue(); 
    return *this; 
  }

  /// Assignment subtraction with compatible quantity.
  template <int DL2, int DE2, int DQ2>
  Qty<L,E,Q,DL,DE,DQ> & 
  operator-=(const Qty<QtyInt<L,DL,DL2>::I,
	               QtyInt<E,DE,DE2>::I,
	               QtyInt<Q,DQ,DQ2>::I,
  	               DL2,DE2,DQ2> x) 
  { 
    rawValue_ -= x.rawValue(); 
    return *this; 
  }

private:
  /// The raw value in units of Qty::baseunit().
  double rawValue_;
};

/// Specialization of Qty for <0,0,0> with conversions to double.
template<int DL, int DE, int DQ>
class Qty<0,0,0,DL,DE,DQ>
{
public:
  /// The squared type.
  typedef Qty<0,0,0,DL,DE,DQ> Squared;

  /// Basic unit of this quantity.
  static double baseunit() {
    return 1.0;
  }

  /// Default constructor to 0.
  Qty(ZeroUnit) : rawValue_(0.0) {}

  /// Default constructor from a double.
  Qty(double x = 0.0) : rawValue_(x) {}

  /// Constructor from a compatible quantity
  template <int DL2, int DE2, int DQ2>
  Qty(const Qty<0,0,0,DL2,DE2,DQ2> & q) : rawValue_(q.rawValue()) {}

  /// Access to the raw value.
  double rawValue() const { return rawValue_; }

  /// Cast to double.
  operator double() const { return rawValue_; }

  /// Assignment multiplication by dimensionless number.
  Qty<0,0,0,DL,DE,DQ> & operator*=(double x) { rawValue_ *= x; return *this; }

  /// Assignment division by dimensionless number.
  Qty<0,0,0,DL,DE,DQ> & operator/=(double x) { rawValue_ /= x; return *this; }

  /// Assignment addition with compatible quantity.
  template <int DL2, int DE2, int DQ2>
  Qty<0,0,0,DL,DE,DQ> & operator+=(const Qty<0,0,0,DL2,DE2,DQ2> x) { 
    rawValue_ += x.rawValue(); 
    return *this; 
  }

  /// Assignment subtraction with compatible quantity.
  template <int DL2, int DE2, int DQ2>
  Qty<0,0,0,DL,DE,DQ> & operator-=(const Qty<0,0,0,DL2,DE2,DQ2> x) { 
    rawValue_ -= x.rawValue(); 
    return *this; 
  }

  /// Assignment addition with double.
  Qty<0,0,0,DL,DE,DQ> & operator+=(double x) { 
    rawValue_ += x; 
    return *this; 
  }

  /// Assignment subtraction with double.
  Qty<0,0,0,DL,DE,DQ> & operator-=(double x) { 
    rawValue_ -= x; 
    return *this; 
  }

private:
  /// The raw value.
  double rawValue_;
};

/// @name Result types for binary operations.
//@{
/**
 * BinaryOpTraits should be specialized with typdefs called MulT and
 * DivT which gives the type resulting when multiplying and dividing
 * the template argument types respectively.
 */
template <typename T, typename U> 
struct BinaryOpTraits;

/** @cond TRAITSPECIALIZATIONS */

template<int L1, int L2, int E1, int E2, int Q1, int Q2,
	 int DL1, int DL2, int DE1, int DE2, int DQ1, int DQ2>
struct BinaryOpTraits<Qty<L1,E1,Q1,DL1,DE1,DQ1>,
		      Qty<L2,E2,Q2,DL2,DE2,DQ2> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Qty<L1*DL2+L2*DL1,E1*DE2+E2*DE1,Q1*DQ2+Q2*DQ1,
              DL1*DL2,DE1*DE2,DQ1*DQ2> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Qty<L1*DL2-L2*DL1,E1*DE2-E2*DE1,Q1*DQ2-Q2*DQ1,
              DL1*DL2,DE1*DE2,DQ1*DQ2> DivT;
};


template<int L1, int E1, int Q1, int DL1, int DE1, int DQ1>
struct BinaryOpTraits<Qty<L1,E1,Q1,DL1,DE1,DQ1>,
		      Qty<L1,E1,Q1,DL1,DE1,DQ1> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Qty<2*L1,2*E1,2*Q1,
              DL1,DE1,DQ1> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef double DivT;
};

/**
 *  Multiplication template
 */
template<int L1, int E1, int Q1, int DL1, int DE1, int DQ1>
struct BinaryOpTraits<double,
		      Qty<L1,E1,Q1,DL1,DE1,DQ1> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Qty<L1,E1,Q1,
              DL1,DE1,DQ1> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Qty<-L1,-E1,-Q1,
              DL1,DE1,DQ1> DivT;
};

/**
 *  Multiplication template
 */
template<int L1, int E1, int Q1, int DL1, int DE1, int DQ1>
struct BinaryOpTraits<Qty<L1,E1,Q1,DL1,DE1,DQ1>,
		      double> {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Qty<L1,E1,Q1,
              DL1,DE1,DQ1> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Qty<L1,E1,Q1,
              DL1,DE1,DQ1> DivT;
};
//@}

/// @name Type traits for alternative code generation.
//@{
/** Type traits for alternative code generation*/
template <int L, int E, int Q, int DL, int DE, int DQ> 
struct TypeTraits<Qty<L,E,Q,DL,DE,DQ> >
{
  /** Enum for dimensions*/
  enum { hasDimension = true };
  /// Type switch set to dimensioned type.
  typedef DimensionT DimType;
  static const Qty<L,E,Q,DL,DE,DQ> baseunit;
};

/** Type traits for alternative code generation*/
template <int L, int E, int Q, int DL, int DE, int DQ> 
const Qty<L,E,Q,DL,DE,DQ> 
TypeTraits<Qty<L,E,Q,DL,DE,DQ> >::baseunit = Qty<L,E,Q,DL,DE,DQ>::baseunit();


/** Type traits for alternative code generation*/
template <int DL, int DE, int DQ> 
struct TypeTraits<Qty<0,0,0,DL,DE,DQ> >
{
  /** Enum for dimensions*/
  enum { hasDimension = false };
  /// Type switch set to standard type.
  typedef StandardT DimType;
  static const double baseunit;
};

/** Type traits for alternative code generation*/
template <int DL, int DE, int DQ> 
const double 
TypeTraits<Qty<0,0,0,DL,DE,DQ> >::baseunit = 1.0;
//@}

/** @endcond */

}

#endif
