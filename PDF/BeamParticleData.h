// -*- C++ -*-
//
// BeamParticleData.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_BeamParticleData_H
#define ThePEG_BeamParticleData_H
// This is the declaration of the BeamParticleData class.

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/PDFBase.h"
#include "BeamParticleData.xh"

namespace ThePEG {

/**
 * BeamParticleData inherits from the ParticleData class and is used
 * for particles which have information about their sub-structure
 * implemented as a pointer to a PDFBase object.
 *
 * @see \ref BeamParticleDataInterfaces "The interfaces"
 * defined for BeamParticleData.
 * @see ParticleData
 * @see PDFBase
 */
class BeamParticleData: public virtual ParticleData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  BeamParticleData() {}
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

public:

  /**
   * Return a pointer to the parton density object describing the
   * sub-structure of this particle type.
   */
  tcPDFPtr pdf() const { return thePDF; }

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
  BeamParticleData(long newId, string newPDGName);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

private:

  /**
   * Set the parton density object.
   */
  void setPDF(PDFPtr);

private:

  /**
   * The pointer to the parton density object.
   */
  PDFPtr thePDF;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<BeamParticleData> initBeamParticleData;

  /**
   *  Private and non-existent assignment operator.
   */
  BeamParticleData & operator=(const BeamParticleData &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BeamParticleData. */
template <>
struct BaseClassTrait<BeamParticleData,1>: public ClassTraitsType {
  /** Typedef of the first base class of BeamParticleData. */
  typedef ParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  BeamParticleData class. */
template <>
struct ClassTraits<BeamParticleData>:
    public ClassTraitsBase<BeamParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::BeamParticleData"; }
};

/** @endcond */

}

#endif /* ThePEG_BeamParticleData_H */
