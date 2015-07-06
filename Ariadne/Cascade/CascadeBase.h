// -*- C++ -*-
#ifndef Ariadne5_CascadeBase_H
#define Ariadne5_CascadeBase_H
//
// This is the declaration of the CascadeBase class.
//

#include "Ariadne/Config/Ariadne5.h"
#include "Ariadne/Config/CloneBase.h"
#include "CascadeBase.fh"
#include "DipoleState.fh"
#include "AriadneHandler.fh"
#include "ThePEG/Utilities/Current.h"

namespace Ariadne5 {

/**
 * CascadeBase is the base class of all Partons, Dipoles and
 * DipoleState classes in the Ariadne dipole cascade classes. It keeps
 * track of the DioleState to which an object belongs.
 */
class CascadeBase: public CloneBase {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor has an optional pointer to the
   * CascadeHandler in charge as argument.
   */
  inline CascadeBase(): isTouched(true) {}

  /**
   * The destructor.
   */
  inline virtual ~CascadeBase() {}
  //@}

public:

  /** @name Functions relating to the DipoleState and CascadeHandler
   *  to which this belongs. */
  //@{
  /**
   * Get the Ariadne::AriadneHandler in charge of the current generation.
   */
  inline tHandlerPtr handler() const {
    return Current<AriadneHandler>::ptr();
  }

  /**
   * Get the DipoleState to which this Object belongs.
   */
  inline tDipoleStatePtr state() const {
    return theState;
  }

protected:

  /**
   * Set the DipoleState to which this Dipole belongs.
   */
  inline void state(tDipoleStatePtr ds) {
    theState= ds;
  }

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

  /** @name Functions determining if the object has been changed since
      the last generation. */
  //@{
  /**
   * If true, this object has been modified since the last round of
   * generating emissions.
   */
  inline bool touched() const {
    return isTouched;
  }

  /**
   * Signal that this object has been modified.
   */
  inline void touch() {
    isTouched = true;
  }

  /**
   * Signal that all possible emissions involving this object has been
   * generated.
   */
  inline void untouch() {
    isTouched = false;
  }
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
   * The DipoleState to which this Dipole belongs.
   */
  tDipoleStatePtr theState;

  /**
   * If true, this object has been modified since the last round of
   * generating emissions.
   */
  bool isTouched;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CascadeBase & operator=(const CascadeBase &);

};

}

#endif /* Ariadne5_CascadeBase_H */
