// -*- C++ -*-
#ifndef ARIADNE5_MinPTScale_H
#define ARIADNE5_MinPTScale_H
//
// This is the declaration of the MinPTScale class.
//

#include "ScaleSetter.h"
#include "DipoleState.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The only task of the MinPTScale class is to calculate a starting
 * scale for the dipole shower from a given DipoleState. This class
 * gives the minimum invariant transverse momentum of all gluons in
 * the sub process as the starting scale
 *
 * @see \ref MinPTScaleInterfaces "The interfaces"
 * defined for MinPTScale.
 */
class MinPTScale: public ScaleSetter {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MinPTScale();

  /**
   * The destructor.
   */
  virtual ~MinPTScale();
  //@}

public:

  /**
   * Return the starting scale for the dipole shower from the given
   * DipoleState. Returns the minimum invariant transverse momentum of
   * gluons in the event (or the total collision energy if no gluons
   * are found.
   */
  virtual Energy scale(const DipoleState &) const;

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
  MinPTScale & operator=(const MinPTScale &);

};

}

#endif /* ARIADNE5_MinPTScale_H */
