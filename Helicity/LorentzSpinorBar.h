// -*- C++ -*-
//
// LorentzSpinorBar.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzSpinorBar_H
#define ThePEG_LorentzSpinorBar_H
// This is the declaration of the LorentzSpinorBar class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "HelicityDefinitions.h"
#include "LorentzSpinor.fh"
#include "LorentzSpinorBar.fh"

namespace ThePEG {
namespace Helicity {

/**
 *  The LorentzSpinorBar class implements the storage of a barred
 *  LorentzSpinor. The design is based on that of the LorentzSpinor
 *  class where the details of the implemented are discussed in more
 *  detail.
 *
 * @see LorentzSpinor
 *
 * @author Peter Richardson
 */

template<typename Value>
class LorentzSpinorBar {
public:

  /** @name Standard constructors. */
  //@{
  /**
   * Default zero constructor, optionally specifying \a t, the type
   */
  LorentzSpinorBar(SpinorType t = SpinorType::unknown) : _type(t), _spin() {}

  /**
   * Constructor with complex numbers specifying the components,
   * optionally specifying \a t, the type
   */
  LorentzSpinorBar(complex<Value> a, complex<Value> b,
		   complex<Value> c, complex<Value> d,
		   SpinorType t = SpinorType::unknown)
    : _type(t), _spin{{a,b,c,d}} {}

  template <typename U>
  LorentzSpinorBar(const LorentzSpinorBar<U> & other)
    : _type(other._type), _spin(other._spin) {}

  //@}

  /** @name Access the components. */
  //@{
  /**
   * Subscript operator to return spinor components
   */
  complex<Value> operator[](int i) const {
    assert( i>= 0 && i <= 3 );
    return _spin[i];
  }

  /**
   * Subscript operator to return spinor components
   */
  complex<Value> operator()(int i) const {
    assert( i>= 0 && i <= 3 );
    return _spin[i];
  }

  /**
   * Set components by index.
   */
  complex<Value> & operator()(int i) {
    assert( i>= 0 && i <= 3 );
    return _spin[i];
  }

  /**
   * Set components by index.
   */
  complex<Value> & operator[](int i) {
    assert( i>= 0 && i <= 3 );
    return _spin[i];
  }

  /**
   * Get first component.
   */
  complex<Value> s1() const {return _spin[0];}

  /**
   * Get second component.
   */
  complex<Value> s2() const {return _spin[1];}

  /**
   * Get third component.
   */
  complex<Value> s3() const {return _spin[2];}

  /**
   * Get fourth component.
   */
  complex<Value> s4() const {return _spin[3];}

  /**
   * Set first component.
   */
  void setS1(complex<Value> in) {_spin[0]=in;}

  /**
   * Set second component.
   */
  void setS2(complex<Value> in) {_spin[1]=in;}

  /**
   * Set third component.
   */
  void setS3(complex<Value> in) {_spin[2]=in;}

  /**
   * Set fourth component.
   */
  void setS4(complex<Value> in) {_spin[3]=in;}
  //@}

  /// @name Mathematical assignment operators.
  //@{
  template <typename ValueB>
  LorentzSpinorBar<Value> & operator+=(const LorentzSpinorBar<ValueB> & a) {
    for(unsigned int ix=0;ix<4;++ix) _spin[ix] += a._spin[ix];
    return *this;
  }
  
  template <typename ValueB>
  LorentzSpinorBar<Value> & operator-=(const LorentzSpinorBar<ValueB> & a) {
    for(unsigned int ix=0;ix<4;++ix) _spin[ix] -= a._spin[ix];
    return *this;
  }
  
  LorentzSpinorBar<Value> & operator*=(double a) {
    for(unsigned int ix=0;ix<4;++ix) _spin[ix] *=a;
    return *this;
  }

  LorentzSpinorBar<Value> & operator/=(double a) {
    for(unsigned int ix=0;ix<4;++ix) _spin[ix] /=a;
    return *this;
  }
  //@}

  /** @name Transformations. */
  //@{
  /**
   * Return the barred spinor
   */
  LorentzSpinor<Value> bar() const;

  /**
   * Return the conjugated spinor \f$u_c=C\bar{u}^T\f$. This operation
   * transforms u-spinors to v-spinors and vice-versa and is required when
   * dealing with majorana particles.
   */
  LorentzSpinorBar conjugate() const;

  /**
   * Standard Lorentz boost specifying the components of the beta vector.
   */
  LorentzSpinorBar & boost(double,double,double);

  /**
   * Standard Lorentz boost specifying the beta vector.
   */
  LorentzSpinorBar & boost(const Boost &);

  /**
   * General Lorentz transformation
   */
  LorentzSpinorBar & transform(const SpinHalfLorentzRotation &) ;
  
  /**
   * General Lorentz transformation
   */
  LorentzSpinorBar & transform(const LorentzRotation & r) {
    transform(r.half());
    return *this;
  }
  //@}

  /** @name Functions related to type. */
  //@{
  /**
   * Return the type of the spinor.
   */
  SpinorType Type() const {return _type;}
  //@}

  /**
   *  @name Functions to apply the projection operator
   */
  //@{
  /**
   *   Apply \f$p\!\!\!\!\!\not\,\,\,+m\f$
   */
  template<typename ValueB> 
  auto projectionOperator(const LorentzVector<ValueB> & p, 
                          const ValueB & m) const 
  -> LorentzSpinorBar<decltype(m*Value())>
  {
    typedef decltype(m*Value()) ResultT;
    LorentzSpinorBar<ResultT> spin;
    static const Complex ii(0.,1.);
    complex<ValueB> p0pp3=p.t()+p.z();
    complex<ValueB> p0mp3=p.t()-p.z();
    complex<ValueB> p1pp2=p.x()+ii*p.y();
    complex<ValueB> p1mp2=p.x()-ii*p.y();
    spin.setS1(m*s1()+p0pp3*s3()+p1pp2*s4());
    spin.setS2(m*s2()+p0mp3*s4()+p1mp2*s3());
    spin.setS3(m*s3()+p0mp3*s1()-p1pp2*s2());
    spin.setS4(m*s4()+p0pp3*s2()-p1mp2*s1());
    return spin;
  }

  /**
   *  Apply \f$g^LP_L+g^RP_R\f$
   */
  LorentzSpinorBar
  helicityProjectionOperator(const Complex & gL, const Complex & gR) const {
    LorentzSpinorBar spin;
    spin.setS1(gL*s1());
    spin.setS2(gL*s2());
    spin.setS3(gR*s3());
    spin.setS4(gR*s4());
    return spin;
  }
  //@}

private:
  /**
   * Type of spinor
   */
  SpinorType _type;

  /**
   * Storage of the components.
   */
  std::array<complex<Value>,4> _spin;
};

/// @name Basic mathematical operations
//@{
template <typename Value>
inline LorentzSpinorBar<double>
operator/(const LorentzSpinorBar<Value> & v, Value a) {
  return LorentzSpinorBar<double>(v.s1()/a, v.s2()/a, v.s3()/a, v.s4()/a,v.Type());
}

inline LorentzSpinorBar<double>
operator/(const LorentzSpinorBar<double> & v, Complex a) {
  return LorentzSpinorBar<double>(v.s1()/a, v.s2()/a, v.s3()/a, v.s4()/a,v.Type());
}

template <typename Value>
inline LorentzSpinorBar<Value> operator-(const LorentzSpinorBar<Value> & v) {
  return LorentzSpinorBar<Value>(-v.s1(),-v.s2(),-v.s3(),-v.s4(),v.Type());
}

template <typename ValueA, typename ValueB>
inline LorentzSpinorBar<ValueA>
operator+(LorentzSpinorBar<ValueA> a, const LorentzSpinorBar<ValueB> & b) {
  return a += b;
}

template <typename ValueA, typename ValueB>
inline LorentzSpinorBar<ValueA>
operator-(LorentzSpinorBar<ValueA> a, const LorentzSpinorBar<ValueB> & b) {
  return a -= b;
}

template <typename Value>
inline LorentzSpinorBar<Value>
operator*(const LorentzSpinorBar<Value> & a, double b) {
  return LorentzSpinorBar<Value>(a.s1()*b, a.s2()*b, a.s3()*b, a.s4()*b,a.Type());
}

template <typename Value>
inline LorentzSpinorBar<Value>
operator*(double b, LorentzSpinorBar<Value> a) {
  return a *= b;
}
  
template <typename Value>
inline LorentzSpinorBar<Value>
operator*(const LorentzSpinorBar<Value> & a, Complex b) {
  return LorentzSpinorBar<Value>(a.s1()*b, a.s2()*b, a.s3()*b, a.s4()*b,a.Type());
}

template <typename ValueA, typename ValueB>
inline auto operator*(complex<ValueB> a, const LorentzSpinorBar<ValueA> & v) 
  -> LorentzSpinorBar<decltype(a.real()*v.s1().real())>
{
  return {a*v.s1(), a*v.s2(), a*v.s3(), a*v.s4(),v.Type()};
}

template <typename ValueA, typename ValueB>
inline auto operator*(const LorentzSpinorBar<ValueA> & v, complex<ValueB> b) 
  -> LorentzSpinorBar<decltype(b.real()*v.s1().real())>
{
  return b*v;
}

template <typename ValueA, typename ValueB>
inline auto operator/(const LorentzSpinorBar<ValueA> & v, complex<ValueB> b) 
  -> LorentzSpinorBar<decltype(v.s1().real()/b.real())>
{
  return {v.s1()/b, v.s2()/b, v.s3()/b, v.s4()/b,v.Type()};
}
  
template <typename ValueA, typename ValueB>
inline auto operator*(ValueB a, const LorentzSpinorBar<ValueA> & v) 
  -> LorentzSpinorBar<decltype(a*v.s1().real())>
{
  return {a*v.s1(), a*v.s2(), a*v.s3(), a*v.s4(),v.Type()};
}

template <typename ValueA, typename ValueB>
inline auto operator*(const LorentzSpinorBar<ValueA> & v, ValueB b) 
  -> LorentzSpinorBar<decltype(b*v.s1().real())>
{
  return b*v;
}

template <typename ValueA, typename ValueB>
inline auto operator/(const LorentzSpinorBar<ValueA> & v, ValueB b) 
  -> LorentzSpinorBar<decltype(v.s1().real()/b)>
{
  return {v.s1()/b, v.s2()/b, v.s3()/b, v.s4()/b,v.Type()};
}

}
}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "LorentzSpinorBar.tcc"
#endif 

#endif
