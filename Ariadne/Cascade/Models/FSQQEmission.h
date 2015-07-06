// -*- C++ -*-
#ifndef ARIADNE5_FSQQEmission_H
#define ARIADNE5_FSQQEmission_H
//
// This is the declaration of the FSQQEmission class.
//

#include "Ariadne/Cascade/Emission.h"
#include "Ariadne/Cascade/QCDDipole.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The FSQQEmission class contains all information about a generated
 * and performed final state splitting of a gluon into a q-qbar pair.
 */
class FSQQEmission: public Emission {

public:

  /** @name Standard constructors and destructors. */
  //@{
   /**
   * The only relevant constructor.
   */
  FSQQEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
	       int infl, Energy inmq)
   : Emission(inmodel, indipole), x1(1.0), x3(1.0), mq(inmq), ifl(infl),
     yo(0.0) {}

  /**
   * The deault constructor should not normally be used.
   */
  FSQQEmission(): x1(1.0), x3(1.0), mq(ZERO), ifl(0), yo(0.0) {}

  /**
   * The destructor.
   */
  virtual ~FSQQEmission();
  //@}

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
   * The energy fraction of the anti-colour emitter (in-parton).
   */
  double x1;

  /**
   * The energy-fraction of the colour emitter (out-parton).
   */
  double x3;

  /**
   * The mass of the quark.
   */
  Energy mq;

  /**
   * The flavour of the quark.
   */
  int ifl;

  /**
   * The rapidity of the anti-quark which is not considered emitted.
   */
  double yo;

  /**
   * The other dipole involved.
   */
  tQCDPtr od;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FSQQEmission & operator=(const FSQQEmission &);

};

}

#endif /* ARIADNE5_FSQQEmission_H */
