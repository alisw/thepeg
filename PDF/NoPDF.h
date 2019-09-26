// -*- C++ -*-
//
// NoPDF.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_NoPDF_H
#define ThePEG_NoPDF_H
// This is the declaration of the NoPDF class.

#include "ThePEG/PDF/PDFBase.h"
// #include "NoPDF.fh"
// #include "NoPDF.xh"

namespace ThePEG {

/**
 * NoPDF inherits from PDFBase and represents particles without
 * sub-structure. The only parton which can be extracted is the
 * incoming particle itself with a momentum distribution which is a
 * delta-function at \f$x=1\f$ (\f$l=0\f$).
 *
 * @see PDFBase
 * @see NoRemnants
 * 
 */
class NoPDF: public PDFBase {

public:

  /** @name Virtual functions mandated by the PDFBase base class. */
  //@{
  /**
   * Return true because we can handle any particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return true if canHandleParticle() and if the corresponding
   * method for remnantHandler() returns true.
   */
  virtual bool canHandle(tcPDPtr particle) const;

  /**
   * Return true if this PDF has a pole at $x=1$ for the given \a
   * particle and \a parton. This default version of the function
   * returns true if \a particle and \a parton is the same.
   */
  virtual bool hasPoleIn1(tcPDPtr particle, tcPDPtr parton) const;

  /**
   * Simply return the particle.
   */
  virtual cPDVector partons(tcPDPtr p) const;

  /**
   * The delta function.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = ZERO) const;
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
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

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<NoPDF> initNoPDF;

  /**
   *  Private and non-existent assignment operator.
   */
  NoPDF & operator=(const NoPDF &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of NoPDF. */
template <>
struct BaseClassTrait<NoPDF,1>: public ClassTraitsType {
  /** Typedef of the first base class of NoPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  NoPDF class. */
template <>
struct ClassTraits<NoPDF>: public ClassTraitsBase<NoPDF> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::NoPDF"; }
};

/** @endcond */

}

#endif /* ThePEG_NoPDF_H */
