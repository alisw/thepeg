// -*- C++ -*-
#ifndef ThePEG_Sum20Momentum_H
#define ThePEG_Sum20Momentum_H
//
// This is the declaration of the Sum20Momentum class.
//

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * Here is the documentation of the Sum20Momentum class.
 */
class Sum20Momentum {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Optionally taking a tolerance level as
   * argument.
   */
  Sum20Momentum(double eps = 1.0e-8): theEps(eps) {
    init();
  }

  /**
   * The default constructor. Optionally taking a tolerance level as
   * argument.
   */
  Sum20Momentum(const LorentzMomentum & p, double eps = 1.0e-8): theEps(eps) {
    init();
    add(p);
  }
  //@}

  void init() {
    theSum = LorentzMomentum();
    theX2 = theY2 = theZ2 = theT2 = ZERO;
  }

  Sum20Momentum & add(const LorentzMomentum & p) {
    theSum += p;
    theX2 += sqr(p.x());
    theY2 += sqr(p.y());
    theZ2 += sqr(p.z());
    theT2 += sqr(p.e());
    return *this;
  }

  Sum20Momentum & operator << (const LorentzMomentum & p) {
    return add(p);
  }

  Sum20Momentum & operator >> (const LorentzMomentum & p) {
    return add(-p);
  }

  Sum20Momentum & operator += (const LorentzMomentum & p) {
    return add(p);
  }

  Sum20Momentum & operator -= (const LorentzMomentum & p) {
    return add(-p);
  }

  bool check() const {
    return sqr(theSum.x()/theEps) < theX2 &&  sqr(theSum.y()/theEps) < theY2
       &&  sqr(theSum.z()/theEps) < theZ2 &&  sqr(theSum.e()/theEps) < theT2;
  }

  bool operator () () const {
    return check();
  }

  bool operator ! () const {
    return !check();
  }

private:

  /**
   * The sum of momenta.
   */
  LorentzMomentum theSum;
  
  /**
   * The summed squared of the x-components
   */
  Energy2 theX2;
  
  /**
   * The summed squared of the y-components
   */
  Energy2 theY2;
  
  /**
   * The summed squared of the z-components
   */
  Energy2 theZ2;
  
  /**
   * The summed squared of the t-components
   */
  Energy2 theT2;
  
  /**
   * The tolerance level.
   */
  double theEps;


};

}

#endif /* ThePEG_Sum20Momentum_H */
