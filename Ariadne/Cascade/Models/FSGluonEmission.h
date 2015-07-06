// -*- C++ -*-
#ifndef ARIADNE5_FSGluonEmission_H
#define ARIADNE5_FSGluonEmission_H
//
// This is the declaration of the FSGluonEmission class.
//

#include "Ariadne/Cascade/Emission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The FSGluonEmission class contains all information about a generated
 * and performed final state gluon emission.
 */
class FSGluonEmission: public Emission {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor.
   */
  FSGluonEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
		  double iny1, double iny3)
    : Emission(inmodel, indipole), x1(1.0), x3(1.0), y1(iny1), y3(iny3) {}

  /**
   * The deault constructor should not normally be used.
   */
  FSGluonEmission(): x1(1.0), x3(1.0), y1(0.0), y3(0.0) {}

  /**
   * The destructor.
   */
  virtual ~FSGluonEmission();
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
   * The scaled mass squared of the anti-colour emitter (in-parton).
   */
  double y1;

  /**
   * The scaled mass squared of the colour emitter (out-parton).
   */
  double y3;

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
  FSGluonEmission & operator=(const FSGluonEmission &);

};

}

#endif /* ARIADNE_FSGluonEmission_H */
