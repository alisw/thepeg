// -*- C++ -*-
#ifndef ARIADNE5_ISQtoGEmission_H
#define ARIADNE5_ISQtoGEmission_H
//
// This is the declaration of the ISQtoGEmission class.
//

#include "ISQEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ISQtoGEmission class contains all information about a generated
 * and performed initial-state splitting of a quark into an
 * initial-state gluon and a final-state quark.
 */
class ISQtoGEmission: public ISQEmission {

public:

  /** @name Standard constructors and destructors. */
  //@{
   /**
   * The only relevant constructor.
   */
  ISQtoGEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
		 tRemParPtr rem, const Lorentz5Momentum & ph)
    : ISQEmission(inmodel, indipole, rem, ph) {}

  /**
   * The deault constructor should not normally be used.
   */
  ISQtoGEmission() {}

  /**
   * The destructor.
   */
  virtual ~ISQtoGEmission();
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
   * The neighboring dipole needed to revert an emission.
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
  ISQtoGEmission & operator=(const ISQtoGEmission &);

};

}

#endif /* ARIADNE5_ISQtoGEmission_H */
