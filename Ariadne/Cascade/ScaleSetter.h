// -*- C++ -*-
#ifndef ARIADNE5_ScaleSetter_H
#define ARIADNE5_ScaleSetter_H
//
// This is the declaration of the ScaleSetter class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ScaleSetter.fh"
#include "DipoleState.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The only task of the ScaleSetter class is to calculate a starting
 * scale for the dipole shower from a given DipoleState. This base
 * class only gives the total collision energy, but other strategies
 * can be implemented in sub-classes overriding the virtual scale()
 * function.
 *
 * @see \ref ScaleSetterInterfaces "The interfaces"
 * defined for ScaleSetter.
 */
class ScaleSetter: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ScaleSetter();

  /**
   * The destructor.
   */
  virtual ~ScaleSetter();
  //@}

public:

  /**
   * Return the starting scale for the dipole shower from the given
   * DipoleState. This base class simply return the total energy of
   * the event. More complicated alternatives must be implemented in
   * sub-classes.
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
  ScaleSetter & operator=(const ScaleSetter &);

};

}

#endif /* ARIADNE5_ScaleSetter_H */
