// -*- C++ -*-
//
// RemnantData.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_RemnantData_H
#define THEPEG_RemnantData_H
//
// This is the declaration of the RemnantData class.
//

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/RemnantData.fh"
#include "ThePEG/PDT/RemnantDecayer.fh"
#include "ThePEG/PDT/DecayMode.h"

namespace ThePEG {

/**
 * The RemnantData class is not a normal ParticleData class. It should
 * never be handled directly by the interface but is automatically
 * created and assigned to an object of the SoftRemnant sub-class of
 * Particle. The SoftRemnant in turn is not a proper Particle, but
 * rather a place holder for what is left of a colliding particle
 * after one or several partons has been extracted. To be able to
 * retrieve properties of the SoftRemnant through its ParticleData the
 * RemnantData is used to dynamically keep track of this.
 *
 * The RemnantData is initialized by the ParticleData corresponding to
 * the colliding particle. For each particle which is extracted the
 * charge charge is changed accordingly. Also the colour charge is
 * changed, but only such that the coloured(), hasColour() and
 * hasAntiColour() returns relevant information. The actual colour can
 * only be singlet, (anti-)triplet or octet.
 *
 * When created the RemnantData object must be given a RemnantDecayer
 * object and a single DecayMode object will be created with this
 * Decayer.
 *
 * @see \ref RemnantDataInterfaces "The interfaces"
 * defined for RemnantData.
 */
class RemnantData: public ParticleData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The standard constructor takes as argument the \a particle type
   * for which this is the remnant and a \a decayer capable of
   * performing the decay.
   */
  RemnantData(tcPDPtr particle, RemDecPtr decayer);
  //@}

public:

  /**
   * The Decayer responsible for for the decay of this remnant.
   */
  const RemnantDecayer & decayer() const {
    return *theDecayer;
  };

  /**
   * Modify the properties to reflect that the given \a parton was
   * extracted.
   */
  bool extract(tcPDPtr parton);

  /**
   * Modify the properties to reflect that the given \a parton which was
   * previously extracted is removed.
   */
  bool remove(tcPDPtr parton);

  /**
   * Modify the properties to reflect that the previously extracted
   * parton, \a oldp, was evolved backwards to the the parton \a newp.
   */
  bool reextract(tcPDPtr oldp, tcPDPtr newp);

protected:
  /**
   * Modify the colour to reflect that the given \a parton was
   * extracted.
   */
  bool fixColour();

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
   * The particle type of the parent.
   */
  tcPDPtr parentPD;

  /**
   * The Decayer responsible for for the decay of this remnant.
   */
  RemDecPtr theDecayer;

  /**
   * The only DecayMode available for this remnant.
   */
  DMPtr decayMode;

  /**
   * The set of extracted particle types.
   */
  multiset<tcPDPtr> extracted;

protected:

  /**
   * The default constructor is protected and must only be used by the
   * PersistentIStream class via the ClassTraits<RemnantData> class.
   */
  RemnantData() {}

  /**
   * The ClassTraits<RemnantData> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<RemnantData>;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<RemnantData> initRemnantData;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RemnantData & operator=(const RemnantData &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RemnantData. */
template <>
struct BaseClassTrait<RemnantData,1> {
  /** Typedef of the first base class of RemnantData. */
  typedef ParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RemnantData class and the shared object where it is defined. */
template <>
struct ClassTraits<RemnantData>
  : public ClassTraitsBase<RemnantData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::RemnantData"; }
  /** Create a Particle object. */
  static TPtr create() { return TPtr::Create(RemnantData()); }
};

/** @endcond */

}

#endif /* THEPEG_RemnantData_H */
