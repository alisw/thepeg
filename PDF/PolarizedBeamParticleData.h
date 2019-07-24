// -*- C++ -*-
//
// PolarizedBeamParticleData.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PolarizedPolarizedBeamParticleData_H
#define ThePEG_PolarizedPolarizedBeamParticleData_H
// This is the declaration of the PolarizedBeamParticleData class.

#include "BeamParticleData.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "PolarizedBeamParticleData.fh"

namespace ThePEG {

/**
 * PolarizedBeamParticleData inherits from the BeamParticleData class and is used
 * for polarized beam particles 
 *
 * @see \ref PolarizedBeamParticleDataInterfaces "The interfaces"
 * defined for PolarizedBeamParticleData.
 * @see BeamParticleData
 * @see PDFBase
 */
class PolarizedBeamParticleData: public virtual BeamParticleData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PolarizedBeamParticleData() : theLongPolarization(0.0) {}
  //@}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(long newId, string newPDGName);

  /**
   * Create a particle - anti particle pair. Note that setting the
   * parton density object on this particle does not change the parton
   * density of the anti particle iven if synchronized() is true.
   */
  static PDPair Create(long newId, string newPDGName, string newAntiPDGName);
  //@}

  /**
   *  Set-up the spin density matrix
   */
  RhoDMatrix rhoMatrix() const;

  /**
   *  The longitudinal polarization
   */
  double longitudinalPolarization() const {return theLongPolarization;}

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
  PolarizedBeamParticleData(long newId, string newPDGName);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

private:

  /**
   *  The longitudinal polarization
   */
  double theLongPolarization;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<PolarizedBeamParticleData> initPolarizedBeamParticleData;

  /**
   *  Private and non-existent assignment operator.
   */
  PolarizedBeamParticleData & operator=(const PolarizedBeamParticleData &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PolarizedBeamParticleData. */
template <>
struct BaseClassTrait<PolarizedBeamParticleData,1>: public ClassTraitsType {
  /** Typedef of the first base class of PolarizedBeamParticleData. */
  typedef BeamParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  PolarizedBeamParticleData class. */
template <>
struct ClassTraits<PolarizedBeamParticleData>:
    public ClassTraitsBase<PolarizedBeamParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::PolarizedBeamParticleData"; }
};

/** @endcond */

}

#endif /* ThePEG_PolarizedBeamParticleData_H */
