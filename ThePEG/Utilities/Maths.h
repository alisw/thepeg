// -*- C++ -*-
//
// Maths.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Math_H
#define ThePEG_Math_H

#include <cmath>

namespace ThePEG {

/** The Math namespace includes the declaration of some useful
 *  mathematical functions. */
namespace Math {

/**
 * MathType is an empty non-polymorphic base class for all
 * mathematical function types.
 */
struct MathType {};

/** The gamma function */
double gamma(double);

/** The log of the gamma function */
double lngamma(double);

/** Return \f${\rm atanh}(x)\f$ */
double atanh(double);

/** Return \f$1-e^x\f$, with highest possible precision for
 *  \f$x\rightarrow 0\f$. */
double exp1m(double x);

/** Return \f$1\log(1-x)\f$, with highest possible precision for
 *  \f$x\rightarrow 0\f$. */
double log1m(double);

/** Return x rased to the integer power p, using recursion. */
double powi(double x, int p);

/** Return the integral of \f$x^p dx\f$ between xl and xu. */
inline double pIntegrate(double p, double xl, double xu) {
  return p == -1.0? log(xu/xl): (pow(xu, p + 1.0) - pow(xl, p + 1.0))/(p + 1.0);
}

/** Return the integral of \f$x^p dx\f$ between xl and xu. */
inline double pIntegrate(int p, double xl, double xu) {
  return p == -1? log(xu/xl): (powi(xu, p + 1) - powi(xl, p + 1))/double(p + 1);
}

/** Return the integral of \f$x^{e-1} dx\f$ between xl and xl+dx with
 *  highest possible precision for \f$dx\rightarrow 0\f$ and/or
 *  \f$e\rightarrow 0\f$. */
inline double pXIntegrate(double e, double xl, double dx) {
  return e == 0.0? log1m(-dx/xl): -pow(xl, e)*exp1m(e*log1m(-dx/xl))/e;
}

/** Generate an x between xl and xu distributed as \f$x^p\f$. */
inline double pGenerate(double p, double xl, double xu, double rnd) {
  return p == -1.0? xl*pow(xu/xl, rnd):
    pow((1.0 - rnd)*pow(xl, p + 1.0) + rnd*pow(xu, p + 1.0), 1.0/(1.0 + p));
}

/** Generate an x between xl and xu distributed as \f$x^p\f$. */
inline double pGenerate(int p, double xl, double xu, double rnd) {
  return p == -1? xl*pow(xu/xl, rnd):
    pow((1.0 - rnd)*powi(xl, p + 1) + rnd*powi(xu, p + 1), 1.0/double(1 + p));
}

/** Generate an x between xl and xl + dx distributed as \f$x^{e-1}\f$
 *  with highest possible precision for\f$dx\rightarrow 0\f$ and/or *
 *  \f$e\rightarrow 0\f$.
 * @param e the parameter defining the power in \f$x^{e-1}\f$.
 * @param xl the lower bound of the generation interval.
 * @param dx the interval.
 * @param rnd a flat random number in the interval ]0,1[. */
inline double pXGenerate(double e, double xl, double dx, double rnd) {
  return e == 0.0? -xl*exp1m(rnd*log1m(-dx/xl)):
    -exp1m(log1m(rnd*exp1m(e*log1m(-dx/xl)))/e)*xl;
}

/** Returns (x - y)/(|x| + |y|). */
template <typename FloatType>
inline double relativeError(FloatType x, FloatType y) {
  return ( x == y ? 0.0 : double((x - y)/(abs(x) + abs(y))) );
}

/** Return x if |x|<|y|, else return y. */
template <typename T>
inline T absmin(const T & x, const T & y) {
  return abs(x) < abs(y)? x: y;
}

/** Return x if |x|>|y|, else return y. */
template <typename T>
inline T absmax(const T & x, const T & y) {
  return abs(x) > abs(y)? x: y;
}

/** Transfer the sign of the second argument to the first.
 * @return \f$|x|\f$ if \f$y>0\f$ otherwise return \f$-|x|\f$.
 */
template <typename T, typename U>
inline T sign(T x, U y) {
  return y > U()? abs(x): -abs(x);
}

/** Templated class for calculating integer powers. */
//@{
/**
 *  Struct for powers
 */
template <int N, bool Inv>
struct Power: public MathType {};

/**
 *  Struct for powers
 */
template <int N>
struct Power<N,false> {
  /** Member for the power*/
  static double pow(double x) { return x*Power<N-1,false>::pow(x); }
};

/**
 *  Struct for powers
 */
template <int N>
struct Power<N,true> {
  /** Member for the power*/
  static double pow(double x) { return Power<N+1,true>::pow(x)/x; }
};

/**
 *  Struct for powers
 */
template <>
struct Power<0,true> {
  /** Member for the power*/
  static double pow(double) { return 1.0; }
};

/**
 *  Struct for powers
 */
template <>
struct Power<0,false> {
  /** Member for the power*/
  static double pow(double) { return 1.0; }
};
//@}

/** Templated function to calculate integer powers known at
 *  compile-time. */
template <int N>
inline double Pow(double x) { return Power<N, (N < 0)>::pow(x); }

/** This namespace introduces some useful function classes with known
 *  primitive and inverse primitive functions. Useful to sample
 *  corresponding distributions.*/
namespace Functions {

/** Class corresponding to functions of the form \f$x^N\f$ with integer N. */
template <int N>
struct PowX: public MathType {

  /** The primitive function. */
  static double primitive(double x) { 
    return Pow<N+1>(x)/double(N+1); 
  }

  /** Integrate function in a given interval. */
  static double integrate(double x0, double x1) {
    return primitive(x1) - primitive(x0);
  }

  /** Sample a distribution in a given interval given a flat random
   *  number R in the interval ]0,1[. */
  static double generate(double x0, double x1, double R) {
    return pow(primitive(x0) + R*integrate(x0, x1), 1.0/double(N+1));
  }

};

/** @cond TRAITSPECIALIZATIONS */

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<1>::generate(double x0, double x1, double R) {
  return std::sqrt(x0*x0 + R*(x1*x1 - x0*x0));
}

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<0>::generate(double x0, double x1, double R) {
  return x0 + R*(x1 - x0);
}

/**
 *  Template for generating according to a specific power
 */
template<>
inline double PowX<-1>::primitive(double x) {
  return log(x);
}

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<-1>::integrate(double x0, double x1) {
  return log(x1/x0);
}

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<-1>::generate(double x0, double x1, double R) {
  return x0*pow(x1/x0, R);
}

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<-2>::generate(double x0, double x1, double R) {
  return x0*x1/(x1 - R*(x1 - x0));
}

/**
 *  Template for generating according to a specific power
 */
template <>
inline double PowX<-3>::generate(double x0, double x1, double R) {
  return x0*x1/std::sqrt(x1*x1 - R*(x1*x1 - x0*x0));
}

/** @endcond */



/** Class corresponding to functions of the form \f$(1-x)^N\f$
 *  with integer N. */
template <int N>
struct Pow1mX: public MathType {

  /** The primitive function. */
  static double primitive(double x) {
    return -PowX<N>::primitive(1.0 - x);
  }

  /** Integrate function in a given interval. */
  static double integrate(double x0, double x1) {
    return PowX<N>::integrate(1.0 - x1, 1.0 - x0);
  }

  /** Sample a distribution in a given interval given a flat random
   *  number R in the interval ]0,1[. */
  static double generate(double x0, double x1, double R) {
    return 1.0 - PowX<N>::generate(1.0 - x1, 1.0 - x0, R);
  }

};

/** Class corresponding to functions of the form \f$1/(x(1-x))\f$ */
struct InvX1mX: public MathType {

  /** The primitive function. */
  static double primitive(double x) {
    return log(x/(1.0 - x));
  }

  /** Integrate function in a given interval. */
  static double integrate(double x0, double x1) {
    return log(x1*(1.0 - x0)/(x0*(1.0 - x1)));
  }

  /** Sample a distribution in a given interval given a flat random
   *  number R in the interval ]0,1[. */
  static double generate(double x0, double x1, double R) {
    double r = pow(x1*(1.0 - x0)/(x0*(1.0 - x1)), R)*x0/(1.0 - x0);
    return r/(1.0 + r);
  }

};

/** Class corresponding to functions of the form \f$e^x\f$ */
struct ExpX: public MathType {

  /** The primitive function. */
  static double primitive(double x) { 
    return exp(x);
  }

  /** Integrate function in a given interval. */
  static double integrate(double x0, double x1) {
    return exp(x1) - exp(x0);
  }

  /** Sample a distribution in a given interval given a flat random
   *  number R in the interval ]0,1[. */
  static double generate(double x0, double x1, double R) {
    return log(exp(x0) + R*(exp(x1) - exp(x0)));
  }

};  

/** Class corresponding to functions of the form \f$x^{N/D}\f$
 *  with integer N and D. */
template <int N, int D>
struct FracPowX: public MathType {

  /** The primitive function. */
  static double primitive(double x) {
    double r = double(N)/double(D) + 1.0;
    return pow(x, r)/r;
  }

  /** Integrate function in a given interval. */
  static double integrate(double x0, double x1) {
    return primitive(x1) - primitive(x0);
  }

  /** Sample a distribution in a given interval given a flat random
   *  number R in the interval ]0,1[. */
  static double generate(double x0, double x1, double R) {
    double r = double(N)/double(D) + 1.0;
    return pow(primitive(x0) + R*integrate(x0, x1), 1.0/r);
  }

};

}

}

}

#endif /* ThePEG_Math_H */
