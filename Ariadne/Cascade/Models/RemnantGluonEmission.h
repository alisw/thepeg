// -*- C++ -*-
#ifndef ARIADNE5_RemnantGluonEmission_H
#define ARIADNE5_RemnantGluonEmission_H
//
// This is the declaration of the RemnantGluonEmission class.
//

#include "FSGluonEmission.h"
#include "Ariadne/Cascade/EmitterBase.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The RemnantGluonEmission class contains all information about a
 * generated and performed final state gluon emission from a remnant
 * parton.
 */
class RemnantGluonEmission: public FSGluonEmission {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The constructor relevant for real emissions.
   */
  RemnantGluonEmission(const EmitterBase & inmodel,
		       const DipoleBase & indipole,
		       double iny1, double iny3)
    : FSGluonEmission(inmodel, indipole, iny1, iny3),
      isRecoil(false) {}

  /**
   * The constructor relevant for recoil emissions.
   */
  RemnantGluonEmission(const EmitterBase & inmodel,
		       const DipoleBase & indipole)
    : FSGluonEmission(inmodel, indipole, 0.0, 0.0), isRecoil(true) {}

  /**
   * The deault constructor should not normally be used.
   */
  RemnantGluonEmission() {}

  /**
   * The destructor.
   */
  virtual ~RemnantGluonEmission();
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
   * True if this corresponds to the emission of a remnant gluon.
   */
  bool isRecoil;

  /**
   * The rotation of the hard sub-system
   */
  mutable LorentzRotation Rh;

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
  RemnantGluonEmission & operator=(const RemnantGluonEmission &);

};

}

#endif /* ARIADNE_RemnantGluonEmission_H */
