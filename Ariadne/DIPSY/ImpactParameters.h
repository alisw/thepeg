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
			  InvEnergy2 w = 1.0/GeV2)
    : theBVec(b), thePhi(angle), theCosPhi(cos(angle)),
      theSinPhi(sin(angle)), theWeight(w) {}

  /**
   * The assignment operator
   */
  inline ImpactParameters & operator=(const ImpactParameters & x) {
    theBVec = x.theBVec;
    thePhi = x.thePhi;
    theCosPhi = x.theCosPhi;
    theSinPhi = x.theSinPhi;
    theWeight = x.theWeight;
    return *this;
  }
  //@}

  /** @name Simple access functions. */
  //@{
  /**
   * The transverse displacement in impact parameter space of one
   * dipole system w.r.t. the other.
   */
  inline const Point & bVec() const {
    return theBVec;
  }

  /**
   * The azimuthal rotation of one dipole system w.r.t. the other.
   */
  inline double phi() const {
    return thePhi;
  }

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
  inline InvEnergy2 weight() const {
    return theWeight;
  }
  //@}

  /**
   * Calculate the squared distance between two partons, assuming the second
   * has been rotated and displaced according to these
   * ImpactParameters.
   */
  inline InvEnergy2 dist2(const Parton & pi, const Parton & pj) const {
    return sqr(pi.position().x() - bVec().x()
	       - pj.position().x()*cosPhi() - pj.position().y()*sinPhi()) +
      sqr(pi.position().y() - bVec().y()
	  - pj.position().y()*cosPhi() + pj.position().x()*sinPhi());
  }

  /**
   * Calculate the squared distance between two partons, assuming the second
   * has been rotated and displaced according to these
   * ImpactParameters.
   */
  inline InvEnergy dist(const Parton & pi, const Parton & pj) const {
    return sqrt(dist2(pi, pj));
  }

  /**
   * Rotates the momentum with the azimuthal angle of the second
   * system wrt the first.
   */
  inline TransverseMomentum rotatePT(const TransverseMomentum p) const {
    return TransverseMomentum( p.x()*cosPhi() + p.y()*sinPhi(),
			       -p.x()*sinPhi() + p.y()*cosPhi() );
  }

  /**
   * Rotate and translate a point in the second system.
   */
  inline Point translate(const Point & p) const {
    return rotate(p) + bVec();
  }

    /**
     * Rotate and translate the given parton.
     */
  inline void translate(tPartonPtr p) const {
    p->position(translate(p->position()));
    p->pT(rotatePT(p->pT()));
  }

  /**
   * Rotates the point with the azimuthal angle of the second system wrt the first.
   */
  inline Point rotate(const Point & p) const {
    return Point(p.x()*cosPhi() + p.y()*sinPhi(),
		 -p.x()*sinPhi() + p.y()*cosPhi());
  }

  /**
   * Rotates the momentum with the azimuthal angle of the first system wrt the second.
   */
    inline TransverseMomentum invRotatePT(const TransverseMomentum & p) const {
      return TransverseMomentum( p.x()*cosPhi() - p.y()*sinPhi(),
			     p.x()*sinPhi() + p.y()*cosPhi() );
    }

  /**
   * Rotates the momentum with the azimuthal angle of the first system wrt the second.
   */
  inline Point invRotate(const Point & p) const {
    return Point( p.x()*cosPhi() - p.y()*sinPhi(),
		  p.x()*sinPhi() + p.y()*cosPhi() );
  }


  /**
   * Calculate the vector from the first to the second point, assuming the second
   * has been rotated and displaced according to these ImpactParamters.
   */
  inline Point difference(const Point & pi, const Point & pj) const {
    return translate(pj) - pi;
  }

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

#endif /* DIPSY_ImpactParameters_H */
