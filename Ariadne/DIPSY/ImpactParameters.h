// -*- C++ -*-
#ifndef DIPSY_ImpactParameters_H
#define DIPSY_ImpactParameters_H
//
// This is the declaration of the ImpactParameters class.
//

#include "Parton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * When two <code>DipoleState</code>s have been generated w.r.t.
 * origo in the transverse, one of them may be displaced and rotated
 * according to an object of the ImpactParameters class before
 * collision. The ImpactParameters are typically generated using an
 * ImpactParameterGenerator object. If not generated according to a
 * flat distribution, a weight may be specified in the
 * ImpactParameters object.
 */
class ImpactParameters {

public:

  /**
   * Use the same Point class as the Parton.
   */
  typedef Parton::Point Point;

public:

  /** @name Standard constructors. */
  //@{
  /**
   * The default constructor, taking the displacement \a b and
   * rotation \a angle as optional arguments. If generated with a
   * non-flat distribution an optional weight \a w can be supplied.
   */
  inline ImpactParameters(const Point & b = Point(), double angle = 0.0,
			  InvEnergy2 w = 1.0/GeV2);

  /**
   * The assignment operator
   */
  inline ImpactParameters & operator=(const ImpactParameters &);
  //@}

  /** @name Simple access functions. */
  //@{
  /**
   * The transverse displacement in impact parameter space of one
   * dipole system w.r.t. the other.
   */
  inline const Point & bVec() const;

  /**
   * The azimuthal rotation of one dipole system w.r.t. the other.
   */
  inline double phi() const;

  /**
   * The cosine of the azimuthal rotation of one dipole system w.r.t. the other.
   */
  inline double cosPhi() const {
    return theCosPhi;
  }

  /**
   * The sine of the azimuthal rotation of one dipole system w.r.t. the other.
   */
  inline double sinPhi() const {
    return theSinPhi;
  }

  /**
   * The weight associated with generating these ImpactParameters.
   */
  inline InvEnergy2 weight() const;
  //@}

  /**
   * Calculate the squared distance between two partons, assuming the second
   * has been rotated and displaced according to these
   * ImpactParameters.
   */
  inline InvEnergy2 dist2(const Parton &, const Parton &) const;

  /**
   * Calculate the squared distance between two partons, assuming the second
   * has been rotated and displaced according to these
   * ImpactParameters.
   */
  inline InvEnergy dist(const Parton &, const Parton &) const;

  /**
   * Rotates the momentum with the azimuthal angle of the second system wrt the first.
   */
  inline TransverseMomentum rotatePT(const TransverseMomentum p) const {
    return TransverseMomentum( p.x()*cosPhi() + p.y()*sinPhi(),
			       -p.x()*sinPhi() + p.y()*cosPhi() );
  }

  /**
   * Rotates the point with the azimuthal angle of the second system wrt the first.
   */
  inline Point rotate(const Point p) const {
    return Point( p.x()*cosPhi() + p.y()*sinPhi(),
		  -p.x()*sinPhi() + p.y()*cosPhi() );
  }

  /**
   * Rotates the momentum with the azimuthal angle of the first system wrt the second.
   */
inline TransverseMomentum invRotatePT(const TransverseMomentum p) const {
  return TransverseMomentum( p.x()*cosPhi() - p.y()*sinPhi(),
			     p.x()*sinPhi() + p.y()*cosPhi() );
}

  /**
   * Rotates the momentum with the azimuthal angle of the first system wrt the second.
   */
inline Point invRotate(const Point p) const {
  return Point( p.x()*cosPhi() - p.y()*sinPhi(),
		p.x()*sinPhi() + p.y()*cosPhi() );
}


  /**
   * Calculate the vector from the first to the second point, assuming the second
   * has been rotated and displaced according to these ImpactParamters.
   */
  inline Point difference(const Point, const Point) const;

private:

  /**
   * The transverse displacement in impact parameter space of one
   * dipole system w.r.t. the other.
   */
  Point theBVec;

  /**
   * The azimuthal rotation of one dipole system w.r.t. the other.
   */
  double thePhi;

  /**
   * The cosine of the azimuthal rotation of one dipole system w.r.t. the other.
   */
  double theCosPhi;

  /**
   * The sine of the azimuthal rotation of one dipole system w.r.t. the other.
   */
  double theSinPhi;

  /**
   * The weight associated with generating these ImpactParameters.
   */
  InvEnergy2 theWeight;

};

}

#include "ImpactParameters.icc"

#endif /* DIPSY_ImpactParameters_H */
