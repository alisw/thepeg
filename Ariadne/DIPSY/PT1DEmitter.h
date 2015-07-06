// -*- C++ -*-
#ifndef DIPSY_PT1DEmitter_H
#define DIPSY_PT1DEmitter_H
//
// This is the declaration of the PT1DEmitter class.
//

#include "Emitter.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the PT1DEmitter class.
 *
 * @see \ref PT1DEmitterInterfaces "The interfaces"
 * defined for PT1DEmitter.
 */
class PT1DEmitter: public DIPSY::Emitter {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PT1DEmitter();

  /**
   * The destructor.
   */
  virtual ~PT1DEmitter();
  //@}

public:

  /** @name Virtual functions which may be overridden by sub-classes. */
  //@{
  /**
   * Generate a possible emission or a swing from a given dipole in the
   * given rapidity interval [\a miny,\a maxy].
   */
  virtual void generate(Dipole & dipole, double miny, double maxy) const;

  /**
   * Perform the emission previously generated for the given \a
   * dipole. If no emission has been generated a runtime_error is
   * thrown.
   */
   virtual void emit(Dipole & dipole) const;
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
  PT1DEmitter & operator=(const PT1DEmitter &);

};

}

#endif /* DIPSY_PT1DEmitter_H */
