// -*- C++ -*-
//
// PDF.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PDF_H
#define ThePEG_PDF_H
// This is the declaration of the PDF class.

#include "ThePEG/PDF/PartonBinInstance.h"

namespace ThePEG {

/**
 * PDF is a simple wrapper class with normal copy-semantics which
 * holds a PDFBase object and a ParticleData object for which to
 * determine parton densities.
 */
class PDF {

public:

  /** @name Standard constructors, assignment and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PDF() {}

  /**
   * Constructor from a given PartonBinInstance.
   */
  PDF(tcPBIPtr pb) {
    if ( !pb ) return;
    thePDF = pb->pdf();
    theParticle = pb->particleData();
  }

  /**
   * Constructor from a given PDFBase and ParticleData object.
   */
  PDF(tcPDFPtr pdf, tcPDPtr pd)
    : thePDF(pdf), theParticle(pd) {}
  //@}

public:

  /** @name Access the parton densities. */
  //@{
  /**
   * Return the density for the given \a parton, for a given \a
   * partonScale and logarithmic momentum fraction \a l assuming the
   * particle has a virtuality \a particleScale.
   */
  double xfl(tcPPtr parton, Energy2 partonScale, double l,
	     Energy2 particleScale = ZERO) const {
    return xfl(parton->dataPtr(), partonScale, l, particleScale);
  }

  /**
   * Return the density for the given \a parton, for a given \a
   * partonScale and momentum fraction \a x assuming the
   * particle has a virtuality \a particleScale.
   */
  double xfx(tcPPtr parton, Energy2 partonScale, double x,
	     double eps = 0.0, Energy2 particleScale = ZERO) const {
    return xfx(parton->dataPtr(), partonScale, x, eps, particleScale);
  }

  /**
   * Return the valence density for the given \a parton, for a given
   * \a partonScale and logarithmic momentum fraction \a l assuming
   * the particle has a virtuality \a particleScale.
   */
  double xfvl(tcPPtr parton, Energy2 partonScale, double l,
	      Energy2 particleScale = ZERO) const {
    return xfvl(parton->dataPtr(), partonScale, l, particleScale);
  }

  /**
   * Return the valence density for the given \a parton, for a given
   * \a partonScale and momentum fraction \a x assuming the particle
   * has a virtuality \a particleScale.
   */
  double xfvx(tcPPtr parton, Energy2 partonScale, double x,
	      double eps = 0.0, Energy2 particleScale = ZERO) const {
    return xfvx(parton->dataPtr(), partonScale, x, eps, particleScale);
  }

  /**
   * Return the density for the given \a parton, for a given \a
   * partonScale and logarithmic momentum fraction \a l assuming the
   * particle has a virtuality \a particleScale.
   */
  double xfl(tcPDPtr parton, Energy2 partonScale, double l,
	     Energy2 particleScale = ZERO) const {
    return thePDF?
      thePDF->xfl(theParticle, parton, partonScale, l, particleScale): 0.0;
  }

  /**
   * Return the density for the given \a parton, for a given \a
   * partonScale and momentum fraction \a x assuming the
   * particle has a virtuality \a particleScale.
   */
  double xfx(tcPDPtr parton, Energy2 partonScale, double x,
	     double eps = 0.0, Energy2 particleScale = ZERO) const {
    return thePDF?
      thePDF->xfx(theParticle, parton, partonScale, x, eps, particleScale): 0.0;
  }

  /**
   * Return the valence density for the given \a parton, for a given
   * \a partonScale and logarithmic momentum fraction \a l assuming
   * the particle has a virtuality \a particleScale.
   */
  double xfvl(tcPDPtr parton, Energy2 partonScale, double l,
	      Energy2 particleScale = ZERO) const {
    return thePDF?
      thePDF->xfvl(theParticle, parton, partonScale, l, particleScale): 0.0;
  }

  /**
   * Return the valence density for the given \a parton, for a given
   * \a partonScale and momentum fraction \a x assuming the particle
   * has a virtuality \a particleScale.
   */
  double xfvx(tcPDPtr parton, Energy2 partonScale, double x,
	      double eps = 0.0, Energy2 particleScale = ZERO) const {
    return thePDF?
      thePDF->xfvx(theParticle, parton, partonScale, x, eps, particleScale): 0.0;
  }
  //@}

  
  /**
   * The parton density object.
   */
  tcPDFPtr pdf() const { return thePDF; }

  /**
   * The particle for which the parton density is used.
   */
  tcPDPtr particle() const { return theParticle; }

  /**
   * Compare for equality.
   */
  bool operator==(const PDF& x) const {
    return
      pdf() == x.pdf() &&
      particle() == x.particle();
  }

  /**
   * Compare for ordering.
   */
  bool operator<(const PDF& x) const {
    return
      pdf() == x.pdf() ?
      particle() < x.particle() :
      pdf() < x.pdf();
  }

private:

  /**
   * The parton density object.
   */
  tcPDFPtr thePDF;

  /**
   * The particle for which the parton density is used.
   */
  tcPDPtr theParticle;

};

}

#endif /* ThePEG_PDF_H */
