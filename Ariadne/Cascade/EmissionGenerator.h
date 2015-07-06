// -*- C++ -*-
#ifndef Ariadne5_EmissionGenerator_H
#define Ariadne5_EmissionGenerator_H
//
// This is the declaration of the EmissionGenerator class.
//

#include "Ariadne/Config/Ariadne5.h"
#include "DipoleBase.fh"
#include "EmitterBase.fh"
#include "Emission.h"
#include <list>

namespace Ariadne5 {

using namespace ThePEG;

/**
 * Here is the documentation of the EmissionGenerator class.
 */
class EmissionGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EmissionGenerator() {}

  /**
   * The main constructor giving the dipole which should be handled.
   */
  EmissionGenerator(tDBPtr d): dipole(d) {}

  /**
   * The destructor.
   */
  ~EmissionGenerator() {}
  //@}

public:

  /**
   * Define ordering (to be used in set<EmissionGenerator>.
   */
  bool operator<(const EmissionGenerator & other) const {
    return dipole < other.dipole;
  }

public:

  /**
   * Initialize by adding all emittors which are able to generate
   * emissions for the dipole. Returns false if no emitters were
   * assigned.
   */
  bool init() const;

  /**
   * Re-initialize if the dipole state has changed. Returns false if
   * no reiitialization was needed or if no emitters were assigned.
   */
  bool reinit() const;

  /**
   * Generate the a phase space point for an emission corresponding to
   * the different emitters.
   */
  Energy generate(Energy rhomin, Energy rhomax) const;

  /**
   * The generated evolution scale.
   */
  Energy rho() const {
    return emission? emission->rho: ZERO;
  }

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
  void persistentInput(PersistentIStream & is);
  //@}

public:

  /**
   * The dipole for which emissions are to be generated.
   */
  tDBPtr dipole;

  /**
   * The emitters which can generate emissions for the dipole.
   */
  mutable list<tEmitterPtr> emitters;

  /**
   * The Emission generated for the dipole.
   */
  mutable EmPtr emission;

  /**
   * Indicate that the dipole has alrady been checked that no further
   * emissions are to be made.
   */
  mutable bool exhausted;

};

}

template <typename OS>
OS & operator<<(OS & os, const Ariadne5::EmissionGenerator & c) {
  c.persistentOutput(os);
  return os;
}

template <typename IS>
IS & operator>>(IS & is, Ariadne5::EmissionGenerator & c) {
  c.persistentInput(is);
  return is;
}

#endif /* Ariadne5_EmissionGenerator_H */
