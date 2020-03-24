// -*- C++ -*-
//
// PhysicalQty.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2019 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Physical_Qty_H
#define Physical_Qty_H
#include "TemplateTools.h"
#include <sstream>
#include <ratio>
#include <type_traits>

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
  constexpr operator double() const { return 0.0; }
};

/// ZERO can be used as zero for any unitful quantity.
constexpr ZeroUnit ZERO = ZeroUnit();

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

// only specialization is with std::ratio below
template <typename L, typename E, typename T>
class Qty;

template <typename T, typename U>
struct qty_equal {
  static constexpr bool value = false;
};

template<typename L1, typename L2, typename E1, typename E2, typename Q1, typename Q2>
struct qty_equal<Qty<L1,E1,Q1>, Qty<L2,E2,Q2>> {
  static constexpr bool value
    =  std::ratio_equal<L1,L2>::value 
    && std::ratio_equal<E1,E2>::value
    && std::ratio_equal<Q1,Q2>::value;
};

template <typename T>
struct is_qty {
  static constexpr bool value = qty_equal<T,T>::value;
};

template <typename ResultT, typename T, typename U = T>
using enable_if_same_qty = typename std::enable_if<qty_equal<T,U>::value, ResultT>::type;

template<long int L, long int E, long int Q, long int DL, long int DE, long int DQ>
class Qty<std::ratio<L,DL>, std::ratio<E,DE>, std::ratio<Q,DQ>>
{
private:
  /// Constructor from raw values. Breaks consistency.
  constexpr Qty(double val) : rawValue_(val) {}

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

  /// General power type
  template <long int Num, long int Den>
  using Power = Qty<typename std::ratio<Num*L,Den*DL>::type, 
                    typename std::ratio<Num*E,Den*DE>::type, 
                    typename std::ratio<Num*Q,Den*DQ>::type>;
  /// Our type
  using Type    = Power<1,1>;
  /// The squared type.
  using Squared = Power<2,1>;
  /// The inverse type.
  using Inverse = Power<-1,1>;
  /// The sqrt type.
  using Sqrt    = Power<1,2>;

  /// Basic unit of this quantity.
  static constexpr Type baseunit() 
  {
    return Type(1.0);
  }

  /// Default constructor to 0.
  constexpr Qty() : rawValue_(0.0) {}

  /// Default constructor to 0.
  constexpr Qty(ZeroUnit) : rawValue_(0.0) {}

  /// Constructor from a compatible quantity
  template <typename U>
  constexpr Qty(const U & q, double factor = 1.0, 
                enable_if_same_qty<void,Type,U> * = nullptr)
    : rawValue_(q.rawValue() * factor) {}

  /// Access to the raw value. Breaks consistency.
  constexpr double rawValue() const { return rawValue_; }

  /// Assignment multiplication by dimensionless number.
  Type & operator*=(double x) { rawValue_ *= x; return *this; }

  /// Assignment division by dimensionless number.
  Type & operator/=(double x) { rawValue_ /= x; return *this; }

  /// Assignment addition with compatible quantity.
  Type & operator+=(const Type & x) 
  { 
    rawValue_ += x.rawValue(); 
    return *this; 
  }

  /// Assignment subtraction with compatible quantity.
  Type & operator-=(const Type & x) 
  { 
    rawValue_ -= x.rawValue(); 
    return *this; 
  }

private:
  /// The raw value in units of Qty::baseunit().
  double rawValue_;
};

/// Specialization of Qty for <0,0,0> with conversions to double.
template<>
class Qty<std::ratio<0>,std::ratio<0>,std::ratio<0>>
{
public:
  /// Our type
  using Type    = Qty<std::ratio<0>,std::ratio<0>,std::ratio<0>>;
  /// General power type
  template <long int Num, long int Den>
  using Power   = Type;
  /// The squared type.
  using Squared = Type;
  /// The inverse type.
  using Inverse = Type;
  /// The sqrt type.
  using Sqrt    = Type;


  /// Basic unit of this quantity.
  static constexpr Type baseunit() {
    return 1.0;
  }

  /// Default constructor to 0.
  constexpr Qty(ZeroUnit) : rawValue_(0.0) {}

  /// Default constructor from a double.
  constexpr Qty(double x = 0.0, double factor=1.0) 
  	: rawValue_(x * factor) {}

  /// Constructor from a compatible quantity
  template <typename U>
  constexpr Qty(const U & q, double factor=1.0, 
                enable_if_same_qty<void,Type,U> * = nullptr)
   : rawValue_(q.rawValue() * factor) {}

  /// Access to the raw value.
  constexpr double rawValue() const { return rawValue_; }

  /// Cast to double.
  constexpr operator double() const { return rawValue_; }

  /// Assignment multiplication by dimensionless number.
  Type & operator*=(double x) { rawValue_ *= x; return *this; }

  /// Assignment division by dimensionless number.
  Type & operator/=(double x) { rawValue_ /= x; return *this; }

  /// Assignment addition with compatible quantity.
  Type & operator+=(const Type & x) { 
    rawValue_ += x.rawValue(); 
    return *this; 
  }

  /// Assignment subtraction with compatible quantity.
  Type & operator-=(const Type & x) { 
    rawValue_ -= x.rawValue(); 
    return *this; 
  }

  /// Assignment addition with double.
  Type & operator+=(double x) { 
    rawValue_ += x; 
    return *this; 
  }

  /// Assignment subtraction with double.
  Type & operator-=(double x) { 
    rawValue_ -= x; 
    return *this; 
  }

private:
  /// The raw value.
  double rawValue_;
};

using QtyDouble = Qty<std::ratio<0>,std::ratio<0>,std::ratio<0>>;

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

template<typename L1, typename L2, 
         typename E1, typename E2, 
         typename Q1, typename Q2>
struct BinaryOpTraits<Qty<L1,E1,Q1>, Qty<L2,E2,Q2>> {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Qty<std::ratio_add<L1,L2>,
              std::ratio_add<E1,E2>,
              std::ratio_add<Q1,Q2>> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Qty<std::ratio_subtract<L1,L2>,
              std::ratio_subtract<E1,E2>,
              std::ratio_subtract<Q1,Q2>> DivT;
};

/**
 *  Multiplication template
 */
template<typename L, typename E, typename Q>
struct BinaryOpTraits<double, Qty<L,E,Q>> {
  /** The type resulting from multiplication of the template type */
  typedef Qty<L,E,Q> MulT;
  /** The type resulting from division of the template type */
  typedef typename BinaryOpTraits<QtyDouble, Qty<L,E,Q>>::DivT DivT;
};

/**
 *  Multiplication template
 */
template<typename L, typename E, typename Q>
struct BinaryOpTraits<Qty<L,E,Q>, double> {
  /** The type resulting from multiplication of the template type */
  typedef Qty<L,E,Q> MulT;
  /** The type resulting from division of the template type */
  typedef Qty<L,E,Q> DivT;
};
//@}

/// @name Type traits for alternative code generation.
//@{
/** Type traits for alternative code generation*/
template <typename L, typename E, typename Q>
struct TypeTraits<Qty<L,E,Q>>
{
  /** Enum for dimensions*/
  enum { hasDimension = true };
  /// Type switch set to dimensioned type.
  typedef DimensionT DimType;
  /// Base unit
  static constexpr Qty<L,E,Q> baseunit() 
	{ return Qty<L,E,Q>::baseunit(); }
};

/** Type traits for alternative code generation*/
template <> 
struct TypeTraits<QtyDouble>
{
  /** Enum for dimensions*/
  enum { hasDimension = false };
  /// Type switch set to standard type.
  typedef StandardT DimType;
  /// Base unit
  static constexpr QtyDouble baseunit() { return 1.0; }
};

//@}

/** @endcond */

}

#endif
