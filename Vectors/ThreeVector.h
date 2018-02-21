// -*- C++ -*-
//
// ThreeVector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2017 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ThreeVector_H
#define ThePEG_ThreeVector_H

/** 
 * @file ThreeVector.h contains the ThreeVector class. ThreeVector can be
 * created with any unit type as template parameter. All basic
 * mathematical operations are supported, as well as a subset of the
 * CLHEP Vector3 functionality.
 */

#include "ThreeVector.fh"
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/UnitIO.h"
#include <cassert>
#include <cmath>

namespace ThePEG {

/** 
 * A 3-component vector. It can be created with any unit type
 * as template parameter.  All basic mathematical operations are
 * supported, as well as a subset of the CLHEP Vector3
 * functionality.
 */
template <typename Value>
class ThreeVector 
{
private:
  /// Value squared
  typedef typename BinaryOpTraits<Value,Value>::MulT Value2;
  /// Value to the 4th power
  typedef typename BinaryOpTraits<Value2,Value2>::MulT Value4;

public:
  /** @name Constructors. */
  //@{
  ThreeVector() 
    : theX(), theY(), theZ() {}

  ThreeVector(Value x, Value y, Value z)
    : theX(x), theY(y), theZ(z) {}

  template<typename ValueB>
  ThreeVector(const ThreeVector<ValueB> & v)
    : theX(v.x()), theY(v.y()), theZ(v.z()) {}
  //@}

public:
  /// @name Component access methods.
  //@{
  Value x() const { return theX; }
  Value y() const { return theY; }
  Value z() const { return theZ; }
  //@}

  /// @name Component set methods.
  //@{
  void setX(Value x)  {  theX = x; }
  void setY(Value y)  {  theY = y; }
  void setZ(Value z)  {  theZ = z; }
  //@}

public:
  /// Squared magnitude \f$x^2+y^2+z^2\f$.
  Value2 mag2() const { return sqr(x()) + sqr(y()) + sqr(z()); }

  /// Magnitude \f$\sqrt{x^2+y^2+z^2}\f$.
  Value  mag() const { return sqrt(mag2()); }

  /// Squared transverse component \f$x^2+y^2\f$.
  Value2 perp2() const { return sqr(x()) + sqr(y()); }

  /// Transverse component \f$\sqrt{x^2+y^2}\f$.
  Value  perp()  const { return sqrt(perp2()); }

  /// Dot product.
  template <typename U>
  typename BinaryOpTraits<Value,U>::MulT
  dot(const ThreeVector<U> & a) const {
    return x()*a.x() + y()*a.y() + z()*a.z();
  }

  /// Squared transverse component with respect to the given axis.
  template <typename U>
  Value2 perp2(const ThreeVector<U> & p) const {
    typedef typename BinaryOpTraits<U,U>::MulT pSqType;
    const pSqType pMag2 = p.mag2();
    assert( pMag2 > pSqType() );
    typename BinaryOpTraits<Value,U>::MulT
      ss = this->dot(p);
    Value2 ret = mag2() - sqr(ss)/pMag2;
    if ( ret <= Value2() )
      ret = Value2();
    return ret;
  }

  /// Transverse component with respect to the given axis.
  template <typename U>
  Value perp(const ThreeVector<U> & p) const {
    return sqrt(perp2(p));
  }

  /// @name Spherical coordinates.
  //@{
  /// Polar angle.
  double theta() const {
    assert(!(x() == Value() && y() == Value() && z() == Value()));
    return atan2(perp(),z());
  }

  /// Azimuthal angle.
  double phi()   const {
    return atan2(y(),x());
  }

  /// Set the polar angle.
  void setTheta(double th) {
    double ma  = mag();
    double ph  = phi();
    setX(ma*sin(th)*cos(ph));
    setY(ma*sin(th)*sin(ph));
    setZ(ma*cos(th));
  }

  /// Set the azimuthal angle.
  void setPhi(double ph) {
    double xy = perp();
    setX(xy*cos(ph));
    setY(xy*sin(ph));
  }
  //@}

  /// Parallel vector with unit length.
  ThreeVector<double> unit() const {
    Value2 mg2 = mag2();
    assert(mg2 > Value2());
    Value mg = sqrt(mg2);
    return ThreeVector<double>(x()/mg, y()/mg, z()/mg);
  }
  
  /// Orthogonal vector.
  ThreeVector<Value> orthogonal() const {
    Value xx = abs(x());
    Value yy = abs(y());
    Value zz = abs(z());
    if (xx < yy) {
      return xx < zz ? ThreeVector<Value>(Value(),z(),-y()) 
	: ThreeVector<Value>(y(),-x(),Value());
    } else {
      return yy < zz ? ThreeVector<Value>(-z(),Value(),x()) 
	: ThreeVector<Value>(y(),-x(),Value());
    }
  }

  /// Azimuthal angle difference, brought into the range \f$(-\pi,\pi]\f$.
  template <typename U>
  double deltaPhi  (const ThreeVector<U> & v2) const {
    double dphi = v2.phi() - phi();
    if ( dphi > Constants::pi ) {
      dphi -= Constants::twopi;
    } else if ( dphi <= -Constants::pi ) {
      dphi += Constants::twopi;
    }
    return dphi;
  } 

  /** 
   * Apply a rotation.
   * @param angle Rotation angle in radians.
   * @param axis Rotation axis.
   */
  template <typename U>
  ThreeVector<Value> & rotate(double angle, const ThreeVector<U> & axis) {
    if (angle == 0.0) 
      return *this;
    const U ll = axis.mag();
    assert( ll > U() );

    const double sa = sin(angle), ca = cos(angle);
    const double dx = axis.x()/ll, dy = axis.y()/ll, dz = axis.z()/ll;
    const Value  xx  = x(), yy = y(), zz = z(); 

    setX((ca+(1-ca)*dx*dx)     * xx
	 +((1-ca)*dx*dy-sa*dz) * yy
	 +((1-ca)*dx*dz+sa*dy) * zz
	 );
    setY(((1-ca)*dy*dx+sa*dz)  * xx
	 +(ca+(1-ca)*dy*dy)    * yy
	 +((1-ca)*dy*dz-sa*dx) * zz
	 );
    setZ(((1-ca)*dz*dx-sa*dy)  * xx
	 +((1-ca)*dz*dy+sa*dx) * yy
	 +(ca+(1-ca)*dz*dz)    * zz
	 );
    return *this;
  }


  /**
   * Rotate the reference frame to a new z-axis.
   */
  ThreeVector<Value> & rotateUz (const Axis & axis) {
    Axis ax = axis.unit();
    double u1 = ax.x();
    double u2 = ax.y();
    double u3 = ax.z();
    double up = u1*u1 + u2*u2;
    if (up>0) {
      up = sqrt(up);
      Value px = x(),  py = y(),  pz = z();
      setX( (u1*u3*px - u2*py)/up + u1*pz );
      setY( (u2*u3*px + u1*py)/up + u2*pz );
      setZ(    -up*px +             u3*pz );
    }
    else if (u3 < 0.) {
      setX(-x());
      setZ(-z()); 
    }
    return *this;
  }

  /**
   * Rotate from a reference frame to the z-axis.
   */
  ThreeVector<Value> & rotateUzBack (const Axis & axis) {
    Axis ax = axis.unit();
    double u1 = ax.x();
    double u2 = ax.y();
    double u3 = ax.z();
    double up = u1*u1 + u2*u2;
    if (up>0) {
      up = sqrt(up);
      Value px = x(),  py = y(),  pz = z();
      setX( ( u1*u3*px + u2*u3*py)/up - up*pz );
      setY( (-u2*px    + u1*py)/up );
      setZ(   u1*px    + u2*py        + u3*pz );
    }
    else if (u3 < 0.) {
      setX(-x());
      setZ(-z()); 
    }
    return *this;
  }

  /// Vector cross-product
  template <typename U>
  ThreeVector<typename BinaryOpTraits<Value,U>::MulT>
  cross(const ThreeVector<U> & a) const {
    typedef ThreeVector<typename BinaryOpTraits<Value,U>::MulT> ResultT;
    return ResultT( y()*a.z()-z()*a.y(),
		   -x()*a.z()+z()*a.x(),
		    x()*a.y()-y()*a.x());
  }
  
  public:  
  /// @name Comparison operators.
  //@{
  bool operator==(const ThreeVector<Value> & a) const {
    return (theX == a.x() && theY == a.y() && theZ == a.z());
  }
  bool operator!=(const ThreeVector<Value> & a) const {
    return !(*this == a);
  }
  bool almostEqual(const ThreeVector<Value> & a, double threshold = 1e-04) const {
    return ((std::abs(theX - a.x()) < threshold) && (std::abs(theY - a.y()) < threshold) && (std::abs(theZ - a.z()) < threshold));
  }
  bool almostUnequal(const ThreeVector<Value> & a, double threshold = 1e-04) const {
    return ! this->almostEqual(a, threshold);
  }
     //@}
  
public:  
  /// @name Mathematical assignment operators.
  //@{
  ThreeVector<Value> & operator+=(const ThreeVector<Value> & a) {
    theX += a.x();
    theY += a.y();
    theZ += a.z();
    return *this;
  }

  ThreeVector<Value> & operator-=(const ThreeVector<Value> & a) {
    theX -= a.x();
    theY -= a.y();
    theZ -= a.z();
    return *this;
  }

  ThreeVector<Value> & operator*=(double a) {
    theX *= a;
    theY *= a;
    theZ *= a;
    return *this;
  }

  ThreeVector<Value> & operator/=(double a) {
    theX /= a;
    theY /= a;
    theZ /= a;
    return *this;
  }
  //@}
  
  /// Cosine of the azimuthal angle between two vectors.
  template <typename U>
  double cosTheta(const ThreeVector<U> & q) const {
    typedef typename BinaryOpTraits<Value,U>::MulT
      ProdType;
    ProdType ptot = mag()*q.mag();
    assert( ptot > ProdType() );
    double arg = dot(q)/ptot;
    if     (arg >  1.0) arg =  1.0;
    else if(arg < -1.0) arg = -1.0;
    return arg;
  }
  
  /// Angle between two vectors.
  template <typename U>
  double angle(const ThreeVector<U> & v) const {
    return acos(cosTheta(v));
  }

private:
  /// @name Vector components
  //@{
  Value theX;
  Value theY;
  Value theZ;
  //@}
};

/// Stream output. Format \f$(x,y,z)\f$.
inline ostream & 
operator<< (ostream & os, const ThreeVector<double> & v)
{
  return os << '(' << v.x() << ',' << v.y() << ',' << v.z() << ')';
}

/// @name Basic mathematical operations
//@{
template <typename Value>
inline ThreeVector<Value>
operator+(ThreeVector<Value> a, 
	  const ThreeVector<Value> & b)
{
  return a += b;
}

template <typename Value>
inline ThreeVector<Value>
operator-(ThreeVector<Value> a, 
	  const ThreeVector<Value> & b)
{
  return a -= b;
}

template <typename Value>
inline ThreeVector<Value> operator-(const ThreeVector<Value> & v) {
  return ThreeVector<Value>(-v.x(),-v.y(),-v.z());
}

template <typename Value>
inline ThreeVector<Value> operator*(ThreeVector<Value> v, double a) {
  return v *= a;
}

template <typename Value>
inline ThreeVector<Value> operator*(double a, ThreeVector<Value> v) {
  return v *= a;
}

template <typename ValueA, typename ValueB>
inline ThreeVector<typename BinaryOpTraits<ValueA,ValueB>::MulT> 
operator*(ValueB a, ThreeVector<ValueA> v) {
  typedef typename BinaryOpTraits<ValueA,ValueB>::MulT ResultT;
  return ThreeVector<ResultT>(a*v.x(), a*v.y(), a*v.z());
}

template <typename ValueA, typename ValueB>
inline ThreeVector<typename BinaryOpTraits<ValueA,ValueB>::MulT> 
operator*(ThreeVector<ValueA> v, ValueB a) {
  return a*v;
}
//@}

/// Vector dot product.
template <typename ValueA, typename ValueB>
inline typename BinaryOpTraits<ValueA,ValueB>::MulT 
operator*(const ThreeVector<ValueA> & a, 
	  const ThreeVector<ValueB> & b)
{
  return a.dot(b);
}

/// A parallel vector with unit length.
template <typename Value>
ThreeVector<double> unitVector(const ThreeVector<Value> & v) {
  return v.unit();
}


/** Output a ThreeVector with units to a stream. */
template <typename OStream, typename UT, typename Value>
void ounitstream(OStream & os, const ThreeVector<Value> & p, UT & u) {
  os << ounit(p.x(), u) << ounit(p.y(), u) << ounit(p.z(), u);
}

/** Input a ThreeVector with units from a stream. */
template <typename IStream, typename UT, typename Value>
void iunitstream(IStream & is, ThreeVector<Value> & p, UT & u) {
  Value x, y, z;
  is >> iunit(x, u) >> iunit(y, u) >> iunit(z, u);
  p = ThreeVector<Value>(x, y, z);
}

}

#endif /* ThePEG_ThreeVector_H */
