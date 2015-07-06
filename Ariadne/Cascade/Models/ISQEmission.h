// -*- C++ -*-
#ifndef ARIADNE5_ISQEmission_H
#define ARIADNE5_ISQEmission_H
//
// This is the declaration of the ISQEmission class.
//

#include "Ariadne/Cascade/Emission.h"
#include "Ariadne/Cascade/QCDDipole.fh"
#include "Ariadne/Cascade/RemnantParton.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ISQEmission class contains all information about a generated
 * and performed initial-state quark emissions.
 */
class ISQEmission: public Emission {

public:

  /** @name Standard constructors and destructors. */
  //@{
   /**
   * The only relevant constructor.
   */
  ISQEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
		 tRemParPtr rem, const Lorentz5Momentum & ph)
    : Emission(inmodel, indipole), q(tcPDPtr()), mq(ZERO), x(rem->x()),
      exorig(&rem->extractedData()), z(1.0), xi(0.0) {
    pold.first = rem->momentum();
    pold.second = ph;
  }

  /**
   * The deault constructor should not normally be used.
   */
  ISQEmission() {}

  /**
   * The destructor.
   */
  virtual ~ISQEmission();
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
   * The type of the emitted quark.
   */
  tcPDPtr q;

  /**
   * The mass of the emitted quark.
   */
  Energy mq;

  /**
   * Original x-value.
   */
  double x;

  /**
   * Original extracted parton.
   */
  tcPDPtr exorig;

  /**
   * The energy splitting.
   */
  double z;

  /**
   * The lightcone fraction taken by the original hard subsystem-
   */
  double xi;

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
  ISQEmission & operator=(const ISQEmission &);

};

}

#endif /* ARIADNE5_ISQEmission_H */
