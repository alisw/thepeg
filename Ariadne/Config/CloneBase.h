// -*- C++ -*-
#ifndef ARIADNE5_CloneBase_H
#define ARIADNE5_CloneBase_H
//
// This is the declaration of the CloneBase class.
//

#include "Ariadne5.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "CloneBase.fh"
#include "Ariadne/DIPSY/DipoleState.fh"
#include "Ariadne/Cascade/DipoleState.fh"

namespace Ariadne5 {

/**
 * CloneBase is used as base class for most of the classes in the
 * Ariadne dipole cascade. It defines the basic clone and rebind
 * methods which are used when cloning a whole dipole state. Maybe
 * this class is general enough to be a part of ThePEG base
 * classes.
 */
class CloneBase: public PersistentBase {

public:

  /**
   * A set of pointers to CloneBase objects.
   */
  typedef set<cClonePtr> CloneSet;

  /**
   * The Rebinder class for CloneBase objects.
   */
  typedef Rebinder<CloneBase> TranslationMap;

  /**
   *
   */
  friend class ::Ariadne5::DipoleState;
  friend class ::DIPSY::DipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline CloneBase() {
    //    allocated[uniqueId] = this;
  }

  /**
   * The default constructor.
   */
  inline CloneBase(const CloneBase & x): PersistentBase(x) {
    //    allocated[uniqueId] = this;
  }

  /**
   * The destructor.
   */
  virtual ~CloneBase();
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name The virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const = 0;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

  /**
   * For debugging purposes
   */
  static map<unsigned long, const CloneBase *> allocated;

  /**
   * For debugging purposes
   */
  long allocount() const;

  /**
   * For debugging purposes
   */
  void allocdebug() const;
    

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CloneBase & operator=(const CloneBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

#endif /* ARIADNE5_CloneBase_H */
