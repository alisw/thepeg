// -*- C++ -*-
//
// UnitIO.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_UnitIO_H
#define ThePEG_UnitIO_H
// This is the declaration of the IUnit and OUnit classes and
// associated templated functions.

#include <complex>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>

namespace ThePEG {

/**
 * The OUnit< class is used to
 * facilitate output of unitful numbers to a
 * persistent stream. An Energy can hence be written like
 * this:<BR> <code>os
 * << ounit(x, GeV);</code><BR> Also containers of unitful
 * numbers may be written like this, as well as LorentzVector and
 * ThreeVector.
 *
 * @see PersistentOStream
 * @see PersistentIStream
 * 
 */
template <typename T, typename UT>
struct OUnit {

  /** Constructor given an object to be written assuming the given
   *  unit. */
  OUnit(const T & t, const UT & u): theX(t), theUnit(u) {}

  /** Copy constructor */
  OUnit(const OUnit<T,UT> & iu): theX(iu.theX), theUnit(iu.theUnit) {}

  /** Reference to the object to be written. */
  const T & theX;

  /** The unit assumed when writing the object. */
  const UT & theUnit;
};

/**
 * The IUnit class is used to facilitate input of unitful numbers from
 * and to a persistent stream. An Energy can hence be read like
 * this:<BR> <code>is >> iunit(x, GeV);</code><BR> Also containers of
 * unitful numbers may be read like this, as well as LorentzVector and
 * ThreeVector.
 *
 * @see PersistentOStream
 * @see PersistentIStream
 * 
 */
template <typename T, typename UT>
struct IUnit {

  /** Constructor given an object to be read assuming the given
   *  unit. */
  IUnit(T & t, const UT & u): theX(t), theUnit(u) {}

  /** Copy constructor */
  IUnit(const IUnit<T,UT> & iu): theX(iu.theX), theUnit(iu.theUnit) {}

  /** Reference to the object to be read. */
  T & theX;

  /** The unit assumed when reading the object. */
  const UT & theUnit;

};

/** Helper function creating a OUnit object given an object and a
 *  unit. */
template <typename T, typename UT>
inline OUnit<T,UT> ounit(const T & t, const UT & ut) {
  return OUnit<T,UT>(t, ut);
}

/** Helper function creating a IUnit object given an object and a
 *  unit. */
template <typename T, typename UT>
inline IUnit<T,UT> iunit(T & t, const UT & ut) {
  return IUnit<T,UT>(t, ut);
}

/** Helper function writing out an object with a given unit to an
 *  output stream. */
template <typename OStream, typename T, typename UT>
void ounitstream(OStream & os, const T & t, UT & u) {
  os << t/u;
}

/** Helper function reading an object with a given unit from an
 *  input stream. */
template <typename IStream, typename T, typename UT>
void iunitstream(IStream & is, T & t, UT & u) {
  double d;
  is >> d;
  t = d*u;;
}

/** Helper function reading a complex object with a given unit from an
 *  input stream. */
template <typename IStream, typename T, typename UT>
void iunitstream(IStream & is, std::complex<T> & t, UT & u) {
  std::complex<double> d;
  is >> d;
  t = d*u;;
}

/** Output an OUnit object to a stream. */
template <typename OStream, typename T, typename UT>
OStream & operator<<(OStream & os, const OUnit<T,UT> & u) {
  ounitstream(os, u.theX, u.theUnit);
  return os;
}

/** Input an IUnit object from a stream. */
template <typename IStream, typename T, typename UT>
IStream & operator>>(IStream & is, const IUnit<T,UT> & u) {
  iunitstream(is, u.theX, u.theUnit);
  return is;
}

/**
 * OUnitErr is used to write out unitful numbers with an error
 * estimate on a standard ostream. using the helper function ouniterr
 * an energy <code>e</code> with an error estimate <code>de</code> can
 * be written out as eg. <code>cout << ouniterr(e, de,
 * GeV);</code>. The result will be presented in scientific format
 * (with the exponent divisible by three) with the relevant number of
 * significant digits with a single digit in parenthesis indicating
 * the error in the least significant digit,
 * eg. <code>1.23(2)e+03</code>.
 */
template <typename T, typename UT>
struct OUnitErr {

  /** Constructor given an object to be written assuming the given
   *  unit. */
  OUnitErr(const T & t, const T & dt, const UT & u): x(t/u), dx(dt/u) {}

  /** The number to be written. */
  double x;

  /** The estimated error of the number to be written. */
  double dx;
  
};

/** Helper function creating a OUnitErr object. */
template <typename T, typename UT>
inline OUnitErr<T,UT> ouniterr(const T & t, const T & dt, const UT & ut) {
  return OUnitErr<T,UT>(t, dt, ut);
}

/** Helper function creating a OUnitErr object. */
inline OUnitErr<double,double> ouniterr(double t, double dt) {
  return OUnitErr<double,double>(t, dt, 1.0);
}

/** Output an OUnitErr object to a stream. */
template <typename OStream, typename T, typename UT>
OStream & operator<<(OStream & os, const OUnitErr<T,UT> & u) {
  if ( ! isfinite(u.x) ) return os << u.x;
  if ( ! isfinite(u.dx) ) {
    ostringstream out;
    out << u.x << '(' << u.dx << ')';
    return os << out.str();
  }
  double dx = min(u.dx, abs(u.x));
  if ( dx <= 0.0 ) return os << u.x;
  double x = abs(u.x);
  ostringstream osse;
  osse << std::scientific << setprecision(0) << dx;
  string sse = osse.str();
  string::size_type ee = sse.find('e');
  long m = static_cast<long>(round(abs(x)/std::pow(10.0,std::atoi(sse.substr(ee + 1).c_str()))));
  int powx = m <= 0? os.precision(): int(log10(double(m)));
  if ( m <= 0 || powx > os.precision() ) sse[0]='0';  
  ostringstream oss;
  oss << std::scientific << setprecision(powx) << x;
  string ss = oss.str();
  string::size_type e = ss.find('e');
  ostringstream out;
  int pp = std::atoi(ss.substr(e + 1).c_str());
  if ( pp%3 == 0 )
    out << ss.substr(0, e) << "(" << sse[0] << ")" << ss.substr(e);
  else if ( (pp - 1)%3 == 0 ) {
    ostringstream oss;
    oss << std::scientific << setprecision(powx) << x/10.0;
    string ss = oss.str();
    string::size_type e = ss.find('e');
    if ( powx == 0 )
      out << ss.substr(0, e) << "0(" << sse[0] << "0)" << ss.substr(e);
    else if ( powx == 1 )
      out << ss.substr(0, ss.find('.'))
	  << ss.substr(ss.find('.') + 1, e - ss.find('.') - 1)
	  << "(" << sse[0] << ")" << ss.substr(e);
    else {
      swap(ss[ss.find('.')], ss[ss.find('.') + 1]);
      out << ss.substr(0, e) << "(" << sse[0] << ")" << ss.substr(e);
    }
  }
  else {
    ostringstream oss;
    oss << std::scientific << setprecision(powx) << x*10.0;
    string ss = oss.str();
    string::size_type e = ss.find('e');
    if ( powx == 0 )
      out << "0." << ss.substr(0, e) << "(" << sse[0] << ")" << ss.substr(e);
    else {
      swap(ss[ss.find('.')], ss[ss.find('.') - 1]);
      out << ss.substr(0, ss.find('.')) << "0" << ss.substr(ss.find('.'), e)
	  << "(" << sse[0] << ")" << ss.substr(e);
    }
  }
  string res = out.str();
  if ( u.x < 0.0 )
    res = "-" + res;
  return os << res;
}

/**
 * The IUnitErr class is used to facilitate input of unitful numbers
 * with error estimates written out using the OUnitErr class.
 * 
 */
template <typename T, typename UT>
struct IUnitErr {

  /** Constructor given an object to be read assuming the given
   *  unit. */
  IUnitErr(T & t, T & dt, const UT & u): x(t), dx(dt), ut(u) {}

  /** Reference to the object to be read. */
  T & x;

  /** The estimated error of the number to be read. */
  T & dx;
  
  /** The unit assumed when reading the object. */
  UT ut;

};

/** Helper function creating a IUnitErr object. */
template <typename T, typename UT>
inline IUnitErr<T,UT> iuniterr(T & t, T & dt, const UT & ut) {
  return IUnitErr<T,UT>(t, dt, ut);
}

/** Helper function creating a OUnitErr object. */
inline IUnitErr<double,double> iuniterr(double & t, double & dt) {
  return IUnitErr<double,double>(t, dt, 1.0);
}

/** Input an IUnit object from a stream. */
template <typename IStream, typename T, typename UT>
IStream & operator>>(IStream & is, const IUnitErr<T,UT> & u) {
  string s;
  double x = 0.0;
  double dx = 0.0;
  double ex = 1.0;
  is >> s;
  string::size_type open = s.find('(');
  string::size_type close = s.find(')');
  string se = "0";
  string sp = "1";
  double pe = 1.0;
  if ( open != string::npos && close != string::npos ) {
    se = s.substr(open + 1);
    sp += s.substr(close + 1);
    string::size_type dot = s.find('.');
    if ( dot != string::npos && dot < open ) pe = std::pow(10.0, 1.0 - (open - dot));
  }

  istringstream(s) >> x;
  istringstream(se) >> dx;
  istringstream(sp) >> ex;

  u.x = x*ex*u.ut;
  u.dx = dx*ex*pe*u.ut;

  return is;
}

}

#endif /* ThePEG_UnitIO_H */
