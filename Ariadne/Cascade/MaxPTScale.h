// -*- C++ -*-
#ifndef ARIADNE5_MaxPTScale_H
#define ARIADNE5_MaxPTScale_H
//
// This is the declaration of the MaxPTScale class.
//

#include "ScaleSetter.h"
#include "DipoleState.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The only task of the MaxPTScale class is to calculate a starting
 * scale for the dipole shower from a given DipoleState. This class
 * gives the maximum invariant transverse momentum of all gluons in
 * the sub process as the starting scale. Alternatively It may give
 * the maximum actual transverse momentum of all partons or the scalar
 * sum of all transverse momenta.
 *
 * @see \ref MaxPTScaleInterfaces "The interfaces"
 * defined for MaxPTScale.
 */
class MaxPTScale: public ScaleSetter {

public:

  /**
   * Enumerate the strategies.
   */
  enum PTStrategy {
    invariantPT = 0, /**< Use the invariant transverse momentum of any gluon. */
    actualPT = 1,    /**< Use the actual transverse momentum of any parton. */
    scalarSum = 2    /**< Use the scalar sum of all transverse momenta. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MaxPTScale();

  /**
   * The destructor.
   */
  virtual ~MaxPTScale();
  //@}

public:

  /**
   * Return the starting scale for the dipole shower from the given
   * DipoleState. Returns the maximum invariant transverse momentum of
   * gluons in the event (or the total collision energy if no gluons
   * are found.
   */
  virtual Energy scale(const DipoleState &) const;

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
   * The strategy chosen for finding the highest transverse momentum.
   */
  PTStrategy strategy;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MaxPTScale & operator=(const MaxPTScale &);

};

}

#endif /* ARIADNE5_MaxPTScale_H */
