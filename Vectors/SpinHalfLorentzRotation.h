// -*- C++ -*-
//
// SpinHalfLorentzRotation.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_SpinHalfLorentzRotation_H
#define THEPEG_SpinHalfLorentzRotation_H
//
// This is the declaration of the SpinHalfLorentzRotation class.
//
#include "ThePEG/Helicity/HelicityDefinitions.h"
#include "ThreeVector.h"

namespace ThePEG {

/**
 * The SpinHalfLorentzRotation class is designed to offer the same
 * features as the HepLorentzRotation class of CLHEP but for the spin-\f$\frac12\f$
 * Lorentz transformation. This is then combined into the general LorentzRotation
 * class of ThePEG to provide the Lorentz transformation for any object as the
 * transformations for higher spin objects can be built from the spin-\f$\frac12\f$
 * and spin-1 transformations.
 *
 * The boost matrix is calculated using the default Dirac matrix representation.
 * Any conversion to other Dirac matrix representations must be handled when the
 * transformation is used.
 */
class SpinHalfLorentzRotation {

  /**
   * The external inverseOf needs to be a friend
   */
  friend SpinHalfLorentzRotation inverseOf ( const SpinHalfLorentzRotation & lt );

public:

  /** @name Constructors and destructor. */
  //@{

  /**
   * Default constructor. Gives a unit matrix.
   */
  SpinHalfLorentzRotation();

  /**
   * Constructor giving the components of a Lorentz boost.
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinHalfLorentzRotation (double bx, double by, double bz, double gamma=-1.);

  /**
   * Constructor giving the vector for a Lorentz boost.
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinHalfLorentzRotation (const Boost & b,double gamma=-1.);
  //@}

  /**
   * Returns true if the Identity matrix.
   */
  bool isIdentity() const;

  /**
   * Return the inverse.
   */
  SpinHalfLorentzRotation inverse() const;

  /**
   * Inverts the SpinHalfLorentzRotation matrix.
   */
  SpinHalfLorentzRotation & invert() { return *this = inverse(); }

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
  SpinHalfLorentzRotation & setBoost (double bx, double by, double bz,double gamma=-1.);

  /**
   * Specify a Lorentz Boost as a vector
   * @param b The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinHalfLorentzRotation & setBoost (const Boost & b,double gamma=-1.);

  /**
   * Specify a boost by the given factor along the x-axis
   * @param boost The Lorentz boost  
   */
  SpinHalfLorentzRotation & setBoostX (double & boost);

  /**
   * Specify a boost by the given factor along the y-axis
   * @param boost The Lorentz boost  
   */
  SpinHalfLorentzRotation & setBoostY (double & boost);

  /**
   * Specify a boost by the given factor along the z-axis
   * @param boost The Lorentz boost  
   */
  SpinHalfLorentzRotation & setBoostZ (double & boost);

  /**
   * Specify a rotation about a general axis by the angle given.
   * @param delta The angle
   * @param axis The axis
   */
  SpinHalfLorentzRotation & setRotate(double delta, const Axis & axis);

  /**
   * Specify a rotation by the given angle about the x-axis
   * @param angle The rotation angle 
   */
  SpinHalfLorentzRotation & setRotateX (double & angle);

  /**
   * Specify a rotation by the given angle about the y-axis
   * @param angle The rotation angle 
   */
  SpinHalfLorentzRotation & setRotateY (double & angle);

  /**
   * Specify a rotation by the given angle about the z-axis
   * @param angle The rotation angle 
   */
  SpinHalfLorentzRotation & setRotateZ (double & angle);
  
  //@}


  /** @name Access methods for the components */
  //@{
  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s1s1() const { return _mx[0][0]; }

  /**
   *   The \f$(1,2)\f$ component
   */
  Complex s1s2() const { return _mx[0][1]; }

  /**
   *   The \f$(1,3)\f$ component
   */
  Complex s1s3() const { return _mx[0][2]; }

  /**
   *   The \f$(1,4)\f$ component
   */
  Complex s1s4() const { return _mx[0][3]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s1() const { return _mx[1][0]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s2() const { return _mx[1][1]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s3() const { return _mx[1][2]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s2s4() const { return _mx[1][3]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s1() const { return _mx[2][0]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s2() const { return _mx[2][1]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s3() const { return _mx[2][2]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s3s4() const { return _mx[2][3]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s1() const { return _mx[3][0]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s2() const { return _mx[3][1]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s3() const { return _mx[3][2]; }

  /**
   *   The \f$(1,1)\f$ component
   */
  Complex s4s4() const { return _mx[3][3]; }

  /**
   *  Fortran style subscript operator
   */
  Complex operator()(unsigned int i, unsigned int j) const {
    assert(i<=3 && j<=3);
    return _mx[i][j];
  }
  //@}


  /** @name Transformation and product members */
  //@{

  /**
   * Product of two SpinHalfLorentzRotations (this) * lt - matrix multiplication  
   * @param lt The SpinHalfLorentzRotation we are multiplying
   */
  SpinHalfLorentzRotation operator * (const SpinHalfLorentzRotation & lt) const;

  /**
   * Multiply by and assign a*=b becomes a= a*b
   */
   SpinHalfLorentzRotation & operator *= (const SpinHalfLorentzRotation & );

  /**
   *  Transform  (similar to *= but a.transform(b) becomes a = b*a
   */
   SpinHalfLorentzRotation & transform   (const SpinHalfLorentzRotation & );

  /**
   * Rotation around the x-axis; equivalent to LT = RotationX(delta) * LT
   */
  SpinHalfLorentzRotation & rotateX(double delta);

  /**
   * Rotation around the y-axis; equivalent to LT = RotationY(delta) * LT
   */
  SpinHalfLorentzRotation & rotateY(double delta);

  /**
   * Rotation around the z-axis; equivalent to LT = RotationZ(delta) * LT
   */
  SpinHalfLorentzRotation & rotateZ(double delta);
  
  /**
   *  Rotation around specified vector - LT = Rotation(delta,axis)*LT
   */
  SpinHalfLorentzRotation & rotate(double delta, const Axis & axis);

  /**
   * Pure boost along the x-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinHalfLorentzRotation & boostX(double beta);

  /**
   * Pure boost along the y-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinHalfLorentzRotation & boostY(double beta);

  /**
   * Pure boost along the z-axis; equivalent to LT = BoostX(beta) * LT
   */
  SpinHalfLorentzRotation & boostZ(double beta);

  /**
   * General boost equivalent to LT = Boost(bx,by,bz) * LT
   * @param bx The x-component of the boost
   * @param by The y-component of the boost
   * @param bz The z-component of the boost
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinHalfLorentzRotation & boost(double bx, double by, double bz, double gamma=-1.);

  /**
   * General boost equivalent to LT = Boost(bv) * LT
   * @param bv The boost vector
   * @param gamma The \f$\gamma\f$ factor (optional)
   */
  SpinHalfLorentzRotation & boost(const Boost & bv, double gamma=-1.);
  //@}

protected:

  /**
   *  Protected constructor giving all the members, no check it is a valid
   *  transformation
   */
  SpinHalfLorentzRotation(Complex s1s1,Complex s1s2,Complex s1s3,Complex s1s4,
				 Complex s2s1,Complex s2s2,Complex s2s3,Complex s2s4,
				 Complex s3s1,Complex s3s2,Complex s3s3,Complex s3s4,
				 Complex s4s1,Complex s4s2,Complex s4s3,Complex s4s4);

private:

  using MatrixT = array<array<Complex,4>,4>;

  SpinHalfLorentzRotation(const MatrixT & m) : _mx(m) {}

  /**
   * The members of the transformation matrix.
   */
  MatrixT _mx;
};

/**
 *  Global method to get the inverse
 */
inline SpinHalfLorentzRotation inverseOf ( const SpinHalfLorentzRotation & lt ) {
  return lt.inverse();
}

/**
 *  output operator
 */
inline std::ostream & operator<< ( std::ostream & os,
				   const  SpinHalfLorentzRotation& lt ) {
  return lt.print(os);
}

}

#endif /* THEPEG_SpinHalfLorentzRotation_H */
