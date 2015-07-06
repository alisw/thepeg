// -*- C++ -*-
#ifndef Ariadne5_StateDipole_H
#define Ariadne5_StateDipole_H
//
// This is the declaration of the StateDipole class.
//

#include "DipoleBase.h"
#include "StateDipole.fh"
#include "Parton.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The StateDipole class is a dummy class to be used by emission
 * models which works on the whole dipole state, rather than on
 * individual dipoles..
 */
class StateDipole: public DipoleBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  StateDipole();

  /**
   * The destructor.
   */
  virtual ~StateDipole() {}
  //@}

protected:

  /** @name Functions relating to the DipoleState to which this belongs. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointer to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  StateDipole & operator=(const StateDipole &);

};

}

#endif /* Ariadne5_StateDipole_H */
