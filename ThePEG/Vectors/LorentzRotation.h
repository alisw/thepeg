// -*- C++ -*-
//
// LorentzRotation.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LorentzRotation_H
#define THEPEG_LorentzRotation_H
//
// This is the declaration of the LorentzRotation class.
//
#include "SpinOneLorentzRotation.h"
#include "SpinHalfLorentzRotation.h"
#include "LorentzRotation.fh"

namespace ThePEG {

/**
 * The LorentzRotation class combine a SpinOneLorentzRotation and a 
 * spin SpinHalfLorentzRotation to provide members which can perform the
 * Lorentz transformation of any object. The class ensures that the 
 * two transformations are consistent by only allowing transformations
 * to be made to both the spin-1 and spin-\f$\frac12\f$ members. 
 */
class LorentzRotation {

  /**
   * The external inverseOf needs to be a friend
   */
  friend LorentzRotation inverseOf ( const LorentzRotation & lt );

public:

  /** @name Constructors and destructor. */
  //@{

  /**
   * Default constructor. Gives a unit matrix.
   */
  LorentzRotation() : _half(), _one() {}

  /**
   * Constructor giving the components of a Lorentz boost.
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation (double bx, double by, double bz, double gamma=-1.) 
    : _half(bx,by,bz,gamma), _one(bx,by,bz,gamma) {}

  /**
   * Constructor giving the vector for a Lorentz boost.
   * @param b The boost vector 
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation (const Boost & b, double gamma=-1.)
    : _half(b,gamma), _one(b,gamma) {}
  //@}

  /**
   * Returns true if the Identity matrix.
   */
  bool isIdentity() const { 
    return _half.isIdentity() && _one.isIdentity(); 
  } 

  /**
   * Return the inverse.
   */
  LorentzRotation inverse() const {
    LorentzRotation output;
    output._half = _half.inverse();
    output._one  =  _one.inverse();
    return output;
  }

  /**
   * Inverts the LorentzRotation matrix.
   */
  LorentzRotation & invert() {
    return *this=inverse();
  }

  /**
   *  output operator
   */
  std::ostream & print( std::ostream & os ) const;

  /** @name Set methods for speical cases of simple rotations and boosts */
  //@{

  /**
   * Specify the components of a Lorentz Boost
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation & setBoost (double bx, double by, double bz, double gamma=-1.) {
    _half.setBoost(bx,by,bz,gamma);
    _one.setBoost(bx,by,bz,gamma);
    return *this;
  }

  /**
   * Specify a Lorentz Boost as a vector
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation & setBoost (const Boost & b, double gamma=-1.) {
    _half.setBoost(b,gamma);
    _one.setBoost(b,gamma);
    return *this;
  }

  /**
   * Specify a boost by the given factor along the x-axis
   * @param boost The Lorentz boost  
   */
  LorentzRotation & setBoostX (double boost) {
    _half.setBoostX(boost);
    _one.setBoost(boost,0,0);
    return *this;
  }

  /**
   * Specify a boost by the given factor along the y-axis
   * @param boost The Lorentz boost  
   */
  LorentzRotation & setBoostY (double boost) {
    _half.setBoostY(boost);
    _one.setBoost(0,boost,0);
    return *this;
  }

  /**
   * Specify a boost by the given factor along the z-axis
   * @param boost The Lorentz boost  
   */
  LorentzRotation & setBoostZ (double boost) {
    _half.setBoostZ(boost);
    _one.setBoost(0,0,boost);
    return *this;
  }

  /**
   * Specify a rotation about a general axis by the angle given.
   * @param delta The angle
   * @param axis The axis
   */
  LorentzRotation & setRotate(double delta, const Axis & axis) {
    _half.setRotate(delta,axis);
    _one.setRotate(delta,axis);
    return *this;
  }

  /**
   * Specify a rotation by the given angle about the x-axis
   * @param angle The rotation angle 
   */
  LorentzRotation & setRotateX (double angle) {
    _half.setRotateX(angle);
    _one.setRotateX(angle);
    return *this;
  }

  /**
   * Specify a rotation by the given angle about the y-axis
   * @param angle The rotation angle 
   */
  LorentzRotation & setRotateY (double angle) {
    _half.setRotateZ(angle);
    _one.setRotateZ(angle);
    return *this;
  }

  /**
   * Specify a rotation by the given angle about the z-axis
   * @param angle The rotation angle 
   */
  LorentzRotation & setRotateZ (double angle) {
    _half.setRotateZ(angle);
    _one.setRotateZ(angle);
    return *this;
  }
  //@}

  /** @name Methods to return the spin-\f$\frac12\f$ and spin-1 transformations */
  //@{

  /**
   * The spin-\f$\frac12\f$ transformation
   */
  const SpinHalfLorentzRotation & half() const { return _half; }

  /**
   * The spin-1 transformation
   */
  const SpinOneLorentzRotation & one() const { return _one; }

  /**
   * Automatically cast to the spin-1 transformation
   */
  operator const SpinOneLorentzRotation & () const { return _one; }
  //@}

  /** @name Access methods for the components of the spin-1 rotation */
  //@{

  /**
   *   The xx component
   */
  double xx() const { return _one.xx(); }

  /**
   *   The xy component
   */
  double xy() const { return _one.xy(); }

  /**
   *   The xz component
   */
  double xz() const { return _one.xz(); }

  /**
   *   The xt component
   */
  double xt() const { return _one.xt(); }

  /**
   *   The yx component
   */
  double yx() const { return _one.yx(); }

  /**
   *   The yy component
   */
  double yy() const { return _one.yy(); }

  /**
   *   The yz component
   */
  double yz() const { return _one.yz(); }

  /**
   *   The yt component
   */
  double yt() const { return _one.yt(); }

  /**
   *   The zx component
   */
  double zx() const { return _one.zx(); }

  /**
   *   The zy component
   */
  double zy() const { return _one.zy(); }

  /**
   *   The zz component
   */
  double zz() const { return _one.zz(); }

  /**
   *   The zt component
   */
  double zt() const { return _one.zt(); }

  /**
   *   The tx component
   */
  double tx() const { return _one.tx(); }

  /**
   *   The ty component
   */
  double ty() const { return _one.ty(); }

  /**
   *   The tz component
   */
  double tz() const { return _one.tz(); }

  /**
   *   The tt component
   */
  double tt() const { return _one.tt(); }
  //@}

  /** @name Access methods for the components of the spin-\f$\frac12\f$ rotation */
  //@{
  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s1s1() const { return _half.s1s1(); }

  /**
   *   The \f$(1,2)\f$ component
   */
  Complex s1s2() const { return _half.s1s2(); }

  /**
   *   The \f$(1,3)\f$ component
   */
  Complex s1s3() const { return _half.s1s3(); }

  /**
   *   The \f$(1,4)\f$ component
   */
  Complex s1s4() const { return _half.s1s4(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s1() const { return _half.s2s1(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s2() const { return _half.s2s2(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s3() const { return _half.s2s3(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s4() const { return _half.s2s4(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s1() const { return _half.s3s1(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s2() const { return _half.s3s2(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s3() const { return _half.s3s3(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s4() const { return _half.s3s4(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s1() const { return _half.s4s1(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s2() const { return _half.s4s2(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s3() const { return _half.s4s3(); }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s4() const { return _half.s4s4(); }
  //@}


  /** @name Transformation and product members */
  //@{

  /**
   * Product with a LorentzVector simply returns the rotated vector.
   */
  template <typename Value>
  LorentzVector<Value>
  operator*(const LorentzVector<Value> & lv) const { return one()*lv; }

  /**
   * Product with a Lorentz5Vector simply returns the rotated vector.
   */
  template <typename Value>
  Lorentz5Vector<Value>
  operator*(const Lorentz5Vector<Value> & lv) const { return one()*lv; }

  /**
   * Product of two LorentzRotations (this) * lt - matrix multiplication  
   * @param lt The LorentzRotation we are multiplying
   */
  LorentzRotation operator * (const LorentzRotation & lt) const {
    LorentzRotation output;
    output._half = _half * lt._half;
    output._one  = _one * lt._one; 
    return output;
  }

  /**
   * Multiply by and assign a*=b becomes a= a*b
   */
  LorentzRotation & operator *= (const LorentzRotation & lt) {
    _one *=lt._one;
    _half*=lt._half;
    return *this;
  }

  /**
   *  Transform  (similar to *= but a.transform(b) becomes a = b*a
   */
  LorentzRotation & transform(const LorentzRotation & lt) {
    _half.transform(lt._half);
    _one.transform(lt._one);
    return *this;
  }

  /**
   * Rotation around the x-axis; equivalent to LT = RotationX(delta) * LT
   */
  LorentzRotation & rotateX(double delta) {
    _half.rotateX(delta);
    _one.rotateX(delta);
    return *this;
  }

  /**
   * Rotation around the y-axis; equivalent to LT = RotationY(delta) * LT
   */
  LorentzRotation & rotateY(double delta) {
    _half.rotateY(delta);
    _one.rotateY(delta);
    return *this;
  }

  /**
   * Rotation around the z-axis; equivalent to LT = RotationZ(delta) * LT
   */
  LorentzRotation & rotateZ(double delta) {
    _half.rotateZ(delta);
    _one.rotateZ(delta);
    return *this;
  }
  
  /**
   *  Rotation around specified vector - LT = Rotation(delta,axis)*LT
   */
  LorentzRotation & rotate(double delta, const Axis & axis) {
    _half.rotate(delta,axis);
    _one.rotate(delta,axis);
    return *this;
  }

  /**
   * Pure boost along the x-axis; equivalent to LT = BoostX(beta) * LT
   */
  LorentzRotation & boostX(double beta) {
    _half.boostX(beta);
    _one.boostX(beta);
    return *this;
  }

  /**
   * Pure boost along the y-axis; equivalent to LT = BoostX(beta) * LT
   */
  LorentzRotation & boostY(double beta) {
    _half.boostY(beta);
    _one.boostY(beta);
    return *this;
  }

  /**
   * Pure boost along the z-axis; equivalent to LT = BoostX(beta) * LT
   */
  LorentzRotation & boostZ(double beta) {
    _half.boostZ(beta);
    _one.boostZ(beta);
    return *this;
  }

  /**
   *  boost equivalent to LT = Boost(bx,by,bz) * LT
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation & boost(double bx, double by, double bz, double gamma=-1.) {
    _half.boost(bx,by,bz,gamma);
    _one.boost(bx,by,bz,gamma);
    return *this;
  }

  /**
   *  boost equivalent to LT = Boost(bv) * LT
   * @param bv The boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  LorentzRotation & boost(const Boost & bv, double gamma=-1.) {
    _half.boost(bv,gamma);
    _one.boost(bv,gamma);
    return *this;
  }
  //@}

private:

  /**
   *  The spin-\f$\frac12\f$ rotation
   */
  SpinHalfLorentzRotation _half;

  /**
   *  The spin-1 rotation
   */
  SpinOneLorentzRotation _one;

};

/**
 *  Global method to get the inverse
 */
inline LorentzRotation inverseOf ( const LorentzRotation & lt ) {
  return lt.inverse();
}

/**
 *  output operator
 */
inline std::ostream & operator<< ( std::ostream & os,
				   const  LorentzRotation& lt ) {
  return lt.print(os);
}

}

#endif /* THEPEG_LorentzRotation_H */

