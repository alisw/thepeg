// -*- C++ -*-
//
// SpinOneLorentzRotation.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SpinOneLorentzRotation_H
#define ThePEG_SpinOneLorentzRotation_H

#include "ThePEG/Helicity/HelicityDefinitions.h"
#include "ThePEG/Helicity/LorentzTensor.fh"
#include "ThePEG/Helicity/LorentzRSSpinor.fh"
#include "ThePEG/Helicity/LorentzRSSpinorBar.fh"
#include "ThreeVector.h"
#include <vector>

namespace ThePEG {

/**
 * The SpinOneLorentzRotation class is ... */

class SpinOneLorentzRotation {
public:

  /** @name Constructors and destructor. */
  //@{

  /**
   * Default constructor. Gives a unit matrix.
   */
  SpinOneLorentzRotation() : matrix_(16) {
    xx_() = yy_() = zz_() = tt_() = 1.0;
  }
  
  /**
   * Constructor giving the components of a Lorentz boost.
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinOneLorentzRotation (double bx, double by, double bz, double gamma=-1.) 
    : matrix_(16) {
    setBoost(bx,by,bz,gamma);
  }

  /**
   * Constructor giving the vector for a Lorentz boost.
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  explicit SpinOneLorentzRotation (const Boost & b, double gamma=-1.)
    : matrix_(16) {
    setBoost(b.x(), b.y(), b.z(),gamma);
  }
  //@}

  /**
   * Returns true if the Identity matrix.
   */
  bool isIdentity() const;

  /**
   * Return the inverse.
   */
  SpinOneLorentzRotation inverse() const;

  /**
   * Inverts the SpinOneLorentzRotation matrix.
   */
  SpinOneLorentzRotation & invert() { return *this = inverse(); }

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
  SpinOneLorentzRotation & setBoost (double bx, double by, double bz, double gamma=-1.);

  /**
   * Specify a Lorentz Boost as a vector
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinOneLorentzRotation & setBoost (const Boost & b, double gamma=-1.) {
    return setBoost(b.x(), b.y(), b.z(),gamma); 
  }

  /**
   * Specify a rotation about a general axis by the angle given.
   * @param delta The angle
   * @param axis The axis
   */
  SpinOneLorentzRotation & setRotate(double delta, const Axis & axis);

  /**
   * Specify a rotation by the given angle about the x-axis
   * @param angle The rotation angle 
   */
  SpinOneLorentzRotation & setRotateX (double angle);

  /**
   * Specify a rotation by the given angle about the y-axis
   * @param angle The rotation angle 
   */
  SpinOneLorentzRotation & setRotateY (double angle);

  /**
   * Specify a rotation by the given angle about the z-axis
   * @param angle The rotation angle 
   */
  SpinOneLorentzRotation & setRotateZ (double angle);
  
  //@}

  /** @name Access methods for the components of the spin-1 rotation */
  //@{

  /**
   *   The xx component
   */
  double xx() const { return matrix_[ 0]; }

  /**
   *   The xy component
   */
  double xy() const { return matrix_[ 1]; }

  /**
   *   The xz component
   */
  double xz() const { return matrix_[ 2]; }

  /**
   *   The xt component
   */
  double xt() const { return matrix_[ 3]; }

  /**
   *   The yx component
   */
  double yx() const { return matrix_[ 4]; }

  /**
   *   The yy component
   */
  double yy() const { return matrix_[ 5]; }

  /**
   *   The yz component
   */
  double yz() const { return matrix_[ 6]; }

  /**
   *   The yt component
   */
  double yt() const { return matrix_[ 7]; }

  /**
   *   The zx component
   */
  double zx() const { return matrix_[ 8]; }

  /**
   *   The zy component
   */
  double zy() const { return matrix_[ 9]; }

  /**
   *   The zz component
   */
  double zz() const { return matrix_[10]; }

  /**
   *   The zt component
   */
  double zt() const { return matrix_[11]; }

  /**
   *   The tx component
   */
  double tx() const { return matrix_[12]; }

  /**
   *   The ty component
   */
  double ty() const { return matrix_[13]; }

  /**
   *   The tz component
   */
  double tz() const { return matrix_[14]; }

  /**
   *   The tt component
   */
  double tt() const { return matrix_[15]; }
  //@}

  /** @name Transformation and product members */
  //@{

  /**
   * Product with a LorentzVector simply returns the rotated vector.
   */
  template <typename Value>
  LorentzVector<Value>
  operator*(const LorentzVector<Value> & v) const {
    return LorentzVector<Value>
      (xx()*v.x() + xy()*v.y() + xz()*v.z() + xt()*v.t(),
       yx()*v.x() + yy()*v.y() + yz()*v.z() + yt()*v.t(),
       zx()*v.x() + zy()*v.y() + zz()*v.z() + zt()*v.t(),
       tx()*v.x() + ty()*v.y() + tz()*v.z() + tt()*v.t());
  }

  /**
   * Product with a Lorentz5Vector simply returns the rotated vector.
   */
  template <typename Value>
  Lorentz5Vector<Value>
  operator*(const Lorentz5Vector<Value> & v) const {
    return Lorentz5Vector<Value>
      (xx()*v.x() + xy()*v.y() + xz()*v.z() + xt()*v.t(),
       yx()*v.x() + yy()*v.y() + yz()*v.z() + yt()*v.t(),
       zx()*v.x() + zy()*v.y() + zz()*v.z() + zt()*v.t(),
       tx()*v.x() + ty()*v.y() + tz()*v.z() + tt()*v.t());
  }

  /**
   * Product of two LorentzRotations (this) * lt - matrix multiplication  
   * @param lt The LorentzRotation we are multiplying
   */
  SpinOneLorentzRotation operator * (const SpinOneLorentzRotation & lt) const;

  /**
   * Multiply by and assign a*=b becomes a= a*b
   */
  SpinOneLorentzRotation & operator *= (const SpinOneLorentzRotation & lt) {
    return *this = *this * lt;
  }

  /**
   *  Transform  (similar to *= but a.transform(b) becomes a = b*a
   */
  SpinOneLorentzRotation & transform   (const SpinOneLorentzRotation & lt) {
    return *this = lt * *this;
  }

  /**
   * Rotation around the x-axis; equivalent to LT = RotationX(delta) * LT
   */
  SpinOneLorentzRotation & rotateX(double delta) {
    SpinOneLorentzRotation tmp;
    tmp.setRotateX(delta);
    return *this = tmp * *this;
  }

  /**
   * Rotation around the y-axis; equivalent to LT = RotationY(delta) * LT
   */
  SpinOneLorentzRotation & rotateY(double delta) {
    SpinOneLorentzRotation tmp;
    tmp.setRotateY(delta);
    return *this = tmp * *this;
  }

  /**
   * Rotation around the z-axis; equivalent to LT = RotationZ(delta) * LT
   */
  SpinOneLorentzRotation & rotateZ(double delta) {
    SpinOneLorentzRotation tmp;
    tmp.setRotateZ(delta);
    return *this = tmp * *this;
  }
  
  /**
   *  Rotation around specified vector - LT = Rotation(delta,axis)*LT
   */
  SpinOneLorentzRotation & rotate(double delta, const Axis & axis) {
    SpinOneLorentzRotation tmp;
    tmp.setRotate(delta, axis);
    return *this = tmp * *this;
  }

  /**
   * Pure boost along the x-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinOneLorentzRotation & boostX(double beta) {
    return *this = SpinOneLorentzRotation(beta,0,0) * *this;
  }

  /**
   * Pure boost along the y-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinOneLorentzRotation & boostY(double beta) {
    return *this = SpinOneLorentzRotation(0,beta,0) * *this;
  }

  /**
   * Pure boost along the z-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinOneLorentzRotation & boostZ(double beta) {
    return *this = SpinOneLorentzRotation(0,0,beta) * *this;
  }

  /**
   *  boost equivalent to LT = Boost(bx,by,bz) * LT
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinOneLorentzRotation & boost(double bx, double by, double bz,
				 double gamma=-1.) {
    return *this = SpinOneLorentzRotation(bx,by,bz,gamma) * *this;
  }

  /**
   *  boost equivalent to LT = Boost(bv) * LT
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinOneLorentzRotation & boost(const Boost & b, double gamma=-1.) {
    return *this = SpinOneLorentzRotation(b.x(),b.y(),b.z(),gamma) * *this;
  }
  //@}

private:

  template<typename Value> friend class Helicity::LorentzTensor;
  template<typename Value> friend class Helicity::LorentzRSSpinor;
  template<typename Value> friend class Helicity::LorentzRSSpinorBar;

  /// Matrix components, order: \f$(xx, xy, \ldots, tz, tt)\f$.
  vector<double> matrix_;

  /// Constructor from doubles.
  SpinOneLorentzRotation (double xx, double xy, double xz, double xt,
			  double yx, double yy, double yz, double yt,
			  double zx, double zy, double zz, double zt,
			  double tx, double ty, double tz, double tt);

  /// Component access by index: x=0, t=3.
  double operator()(unsigned int i, unsigned int j) const {
    return matrix_[4*i + j];
  }

  /// @name Component access.
  //@{
  double & xx_() { return matrix_[ 0]; }
  double & xy_() { return matrix_[ 1]; }
  double & xz_() { return matrix_[ 2]; }
  double & xt_() { return matrix_[ 3]; }

  double & yx_() { return matrix_[ 4]; }
  double & yy_() { return matrix_[ 5]; }
  double & yz_() { return matrix_[ 6]; }
  double & yt_() { return matrix_[ 7]; }

  double & zx_() { return matrix_[ 8]; }
  double & zy_() { return matrix_[ 9]; }
  double & zz_() { return matrix_[10]; }
  double & zt_() { return matrix_[11]; }

  double & tx_() { return matrix_[12]; }
  double & ty_() { return matrix_[13]; }
  double & tz_() { return matrix_[14]; }
  double & tt_() { return matrix_[15]; }
  //@}
};

/**
 *  output operator
 */
inline std::ostream & operator<< ( std::ostream & os,
				   const  SpinOneLorentzRotation& lt ) {
  return lt.print(os);
}

}

#endif /* ThePEG_SpinOneLorentzRotation_H */
