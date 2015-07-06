// -*- C++ -*-
#ifndef DIPSY_GaussianImpactGenerator_H
#define DIPSY_GaussianImpactGenerator_H
//
// This is the declaration of the GaussianImpactGenerator class.
//

#include "ImpactParameterGenerator.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The objective of the GaussianImpactGenerator class is to generate
 * GaussianImpacts objects to be used when two
 * <code>DipoleState</code>s have been generated w.r.t.  origo in the
 * transverse to displace and rotate one of them before
 * collision. This base class will generate impact parameters
 * according to a Gaussian distribution, and the weigth in the
 * produced GaussianImpacts objects is set accordingly. Sub-classes
 * may override the generate() function to use any distribution as
 * long as the weight is set accordingly.
 *
 * @see \ref GaussianImpactGeneratorInterfaces "The interfaces"
 * defined for GaussianImpactGenerator.
 */
class GaussianImpactGenerator: public ImpactParameterGenerator {

public:

  /**
   * Use the same Point class as the Parton.
   */
  typedef Parton::Point Point;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GaussianImpactGenerator();

  /**
   * The destructor.
   */
  virtual ~GaussianImpactGenerator();
  //@}

public:

  /**
   * Generate an GaussianImpacts object.
   */
  virtual ImpactParameters generate(double seed = UseRandom::rnd()) const;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GaussianImpactGenerator & operator=(const GaussianImpactGenerator &);

};

}

#endif /* DIPSY_GaussianImpactGenerator_H */
