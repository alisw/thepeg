// -*- C++ -*-
#ifndef Ariadne5_ResonanceParton_H
#define Ariadne5_ResonanceParton_H
//
// This is the declaration of the ResonanceParton class.
//

#include "ResonanceParton.fh"
#include "Resonance.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ResonanceParton class inherits from Parton and represent a
 * final state parton which is a decay product of a resonance.
 */
class ResonanceParton: public Parton {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ResonanceParton();

  /**
   * The destructor.
   */
  virtual ~ResonanceParton();
  //@}

public:

  /**
   * Set the parent resonance.
   */
  void setResonance(tResPtr r) {
    theResonance = r;
    origSystem(r->decaySystem());
  }

  /**
   * Add a new final state sibling.
   */
  void addSibling(tParPtr s) {
    theSiblings.insert(s);
  }

  /**
   * Add a new intermediate sibling. If previously present, it will be
   * removed from the final state.
   */
  void addOther(tParPtr o);

  /**
   * Get the parent resonance.
   */
  tResPtr resonance() const {
    return theResonance;
  }

  /**
   * Get the final state siblings.
   */
  const set<tParPtr> & siblings() const {
    return theSiblings;
  }

  /**
   * Get intermediate siblings.
   */
  const set<tParPtr> & others() const {
    return theOthers;
  }

  /**
   * We need to know if an emission has been performed from a sibling.
   */
  virtual void notify(const Emission &);

  /**
   * Return true if this is not a normal parton in a dipole together
   * with an \a other parton.
   */
  virtual bool special(tcParPtr other = tcParPtr()) const {
    return !other || other->system() != system();
  }

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
   * The parent resonance.
   */
  tResPtr theResonance;

  /**
   * The sibling final state partons from the same resonance decay
   * (including their emitted partons).
   */
  set<tParPtr> theSiblings;

  /**
   * Other decay products of the same resonance decay.
   */
  set<tParPtr> theOthers;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ResonanceParton & operator=(const ResonanceParton &);

};

}

#endif /* Ariadne5_ResonanceParton_H */
