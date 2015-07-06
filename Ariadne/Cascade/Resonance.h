// -*- C++ -*-
#ifndef Ariadne5_Resonance_H
#define Ariadne5_Resonance_H
//
// This is the declaration of the Resonance class.
//

#include "Parton.h"
#include "Resonance.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * Here is the documentation of the Resonance class.
 */
class Resonance: public Parton {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Resonance();

  /**
   * The destructor.
   */
  virtual ~Resonance();
  //@}

  /** @name Simple access functions. */
  //@{
  /**
   * The sub system index for which the decay products of this
   * resonance belong.
   */
  int decaySystem() const {
    return theDecaySystem;
  }

  /**
   * Return the parent resonance if it exists.
   */
  tResPtr parentResonance() const {
    return theParentResonance;
  }

  /**
   * Set the parent resonance if it exists.
   */
  void parentResonance(tResPtr r) {
    theParentResonance = r;
    origSystem(r->decaySystem());
  }

  /**
   * Set the sub system index for which the decay products of this
   * resonance belong.
   */
  void decaySystem(int isys) {
    theDecaySystem = isys;
  }
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
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

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

private:

  /**
   * The parent resonance if it exists.
   */
  tResPtr theParentResonance;

  /**
   * The sub system index for which the decay products of this
   * resonance belong.
   */
  int theDecaySystem;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Resonance & operator=(const Resonance &);

};

}

#endif /* Ariadne5_Resonance_H */
