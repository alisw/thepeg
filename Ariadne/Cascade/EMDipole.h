// -*- C++ -*-
#ifndef Ariadne5_EMDipole_H
#define Ariadne5_EMDipole_H
//
// This is the declaration of the EMDipole class.
//

#include "DipoleBase.h"
#include "EMDipole.fh"
#include "Parton.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The EMDipole class represents electromagnetic dipoles between partons.
 */
class EMDipole: public DipoleBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EMDipole();

  /**
   * The destructor.
   */
  virtual ~EMDipole() {}
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

public:

  /** @name Simple access functions. */
  //@{
  /**
   * The first parton in this dipole.
   */
  inline tParPtr iPart() const {
    return theIPart;
  }

  /**
   * The second parton in this dipole.
   */
  inline tParPtr oPart() const {
    return theOPart;
  }

  /**
   * Set the first parton in this dipole.
   */
  inline void iPart(tParPtr x) {
    theIPart = x;
  }

  /**
   * Set the second parton in this dipole.
   */
  inline void oPart(tParPtr x) {
    theOPart = x;
  }

  /**
   * Return the squared invariant mass of this dipole.
   */
  Energy2 sdip() const;

  /**
   * The sub system index for which the decay products of this
   * resonance belong.
   */
  int system() const {
    return theSystem;
  }

  /**
   * Set the sub system index for which the decay products of this
   * resonance belong.
   */
  void system(int isys) {
    theSystem = isys;
  }
  //@}

private:

  /**
   * The first parton in this dipole.
   */
  tParPtr theIPart;

  /**
   * The second parton in this dipole.
   */
  tParPtr theOPart;

  /**
   * The sub system index for which the decay products of this
   * resonance belong.
   */
  int theSystem;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EMDipole & operator=(const EMDipole &);

};

}

#endif /* Ariadne5_EMDipole_H */
