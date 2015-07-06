// -*- C++ -*-
#ifndef DIPSY_FixedImpactGenerator_H
#define DIPSY_FixedImpactGenerator_H
//
// This is the declaration of the FixedImpactGenerator class.
//

#include "ImpactParameterGenerator.h"
#include "ImpactParameters.h"
#include "ThePEG/Repository/UseRandom.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The objective of the FixedImpactGenerator class is to generate
 * ImpactParameters objects to be used when two
 * <code>DipoleState</code>s have been generated w.r.t.  origo in the
 * transverse to displace and rotate one of them before
 * collision. This class will generate the absolute value of the
 * impact parameter to be within a given interval. Note that the cross
 * sections produced with this objects will not be trustworthy, but
 * the distibution of events will.
 *
 * @see \ref FixedImpactGeneratorInterfaces "The interfaces"
 * defined for FixedImpactGenerator.
 */
class FixedImpactGenerator: public ImpactParameterGenerator {

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
  FixedImpactGenerator();

  /**
   * The copy constructor.
   */
  FixedImpactGenerator(const FixedImpactGenerator &);

  /**
   * The destructor.
   */
  virtual ~FixedImpactGenerator();
  //@}

public:

  /**
   * Generate an ImpactParameters object.
   */
  virtual ImpactParameters generate(double seed = UseRandom::rnd()) const;

  /**
   * The maximum value of the impact parameter.
   */
  InvEnergy maxB() const {
    return theMaxB/Constants::hbarc;
  }

  /**
   * The minimum value of the impact parameter.
   */
  InvEnergy minB() const {
    return theMinB/Constants::hbarc;
  }

  /**
   * The angle of b (-1 is random).
   */
  double bAngle() const {
    return theBAngle;
  }

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

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
   * The minimum value of the impact parameter.
   */
  Length theMinB;

  /**
   * The maximum value of the impact parameter.
   */
  Length theMaxB;

  /**
   * The angle of b (-1 is random).
   */
  double theBAngle;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FixedImpactGenerator & operator=(const FixedImpactGenerator &);

};

}

#endif /* DIPSY_FixedImpactGenerator_H */
