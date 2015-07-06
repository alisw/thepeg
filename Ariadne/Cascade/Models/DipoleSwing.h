// -*- C++ -*-
#ifndef ARIADNE5_DipoleSwing_H
#define ARIADNE5_DipoleSwing_H
//
// This is the declaration of the DipoleSwing class.
//

#include "Ariadne/Cascade/Emission.h"
#include "Ariadne/Cascade/StateDipole.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The DipoleSwing class contains all information about a generated
 * Swing.
 */
class DipoleSwing: public Emission {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor.
   */
  DipoleSwing(const EmitterBase & inmodel, const DipoleBase & indipole,
	      QCDDipole & d1, QCDDipole & d2, Time dt)
    : Emission(inmodel, indipole), forced(false) {
    setup(d1, d2, dt);
  }

  /**
   * The deault constructor should not normally be used.
   */
  DipoleSwing() {}

  /**
   * The destructor.
   */
  virtual ~DipoleSwing();
  //@}

public:

  void setup(QCDDipole & d1, QCDDipole & d2, Time dt);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual ClonePtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual ClonePtr fullclone() const;
  //@}


public:

  /**
   * The the pair of dipoles that should swing.
   */
  pair<tQCDPtr,tQCDPtr> dipoles;

  /**
   * The colour index of the dipoles that should swing.
   */
  ColourIndex index;

  /**
   * This swing was forced on a tiny string.
   */
  bool forced;

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

  /**
   * Standard debug function to be called from within a debugger.
   */
  void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleSwing & operator=(const DipoleSwing &);

};

}

#endif /* ARIADNE_DipoleSwing_H */
