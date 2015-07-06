// -*- C++ -*-
//
// ConstituentParticleData.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ConstituentParticleData_H
#define ThePEG_ConstituentParticleData_H
// This is the declaration of the ConstituentParticleData class.

#include "ThePEG/PDT/ParticleData.h"

namespace ThePEG {

/**
 * ConstituentParticleData inherits from the ParticleData class and is
 * used for quarks, diquarks and gluons to store information about
 * their constituent mass.
 *
 * @see \ref ConstituentParticleDataInterfaces "The interfaces"
 * defined for ConstituentParticleData.
 */
class ConstituentParticleData: public virtual ParticleData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ConstituentParticleData()
    : theConstituentMass(ZERO), theDefaultConstituentMass(ZERO) {}
  //@}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(long newId, string newPDGName);

  /**
   * Create a particle - anti particle pair.
   */
  static PDPair Create(long newId, string newPDGName, string newAntiPDGName);
  //@}

public:

  /**
   * Return the constituent mass of this parton.
   */
  virtual Energy constituentMass() const { return theConstituentMass; }

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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /**
   * Protected constructor only to be used by subclasses or by the
   * Create method.
   */
  ConstituentParticleData(long newId, string newPDGName);

  /**
   * Read setup info from a standard stream. The information that must
   * be supplied is the same as for ParticleData::readSetup with an
   * additional constituent mass (in GeV) added in the end.
   */
  virtual void readSetup(istream & is);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

private:

  /**
   * Utility function for the interface.
   */
  void setConstituentMass(Energy m);

  /**
   * Utility function for the interface.
   */
  Energy defConstituentMass() const;

private:

  /**
   * The constituent mass of this parton.
   */
  Energy theConstituentMass;

  /**
   * The default constituent mass of this parton.
   */
  Energy theDefaultConstituentMass;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ConstituentParticleData> initConstituentParticleData;

  /**
   *  Private and non-existent assignment operator.
   */
  ConstituentParticleData & operator=(const ConstituentParticleData &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ConstituentParticleData. */
template <>
struct BaseClassTrait<ConstituentParticleData,1>: public ClassTraitsType {
  /** Typedef of the first base class of ConstituentParticleData. */
  typedef ParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  ConstituentParticleData class. */
template <>
struct ClassTraits<ConstituentParticleData>:
    public ClassTraitsBase<ConstituentParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ConstituentParticleData"; }
};

/** @endcond */

}

#endif /* ThePEG_ConstituentParticleData_H */
