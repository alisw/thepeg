// -*- C++ -*-
//
// PDFBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PDFBase_H
#define ThePEG_PDFBase_H
// This is the declaration of the PDFBase class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/PDF/PDFCuts.h"
#include "PDFBase.xh"

namespace ThePEG {

/**
 * PDFBase is the base class for implementing parton density functions
 * for particles with sub-structure. A number of of virtual methods
 * are defined which should be overridden by sub-classes.
 *
 * It is essential that either xfx or xfl is overidden to avoid
 * infinite recursive function calls.
 *
 * A PDFBase object can be assigned to a BeamParticleData object
 * and/or to a PartonExtractor object. A PDFBase has a pointer to a
 * RemnantHandler object which should be capable of generating
 * remnants for all partons which may be extracted by the PDF.
 *
 * @see \ref PDFBaseInterfaces "The interfaces"
 * defined for PDFBase.
 * @see BeamParticleData
 * @see PartonExtractor
 * @see RemnantHandler
 * @see PDFCuts
 */
class PDFBase: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PDFBase();

  /**
   * Copy-constructor.
   */
  PDFBase(const PDFBase &);

  /**
   * Destructor.
   */
  virtual ~PDFBase();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of partons from
   * the given \a particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const = 0;

  /**
   * Return true if canHandleParticle() and if the corresponding
   * method for remnantHandler() returns true for the given \a
   * particle.
   */
  virtual bool canHandle(tcPDPtr particle) const;

  /**
   * Return true if this PDF has a pole at $x=1$ for the given \a
   * particle and \a parton. This default version of the function
   * returns false.
   */
  virtual bool hasPoleIn1(tcPDPtr particle, tcPDPtr parton) const;

  /**
   * Return the partons which this PDF may extract from the given
   * \a particle.
   */
  virtual cPDVector partons(tcPDPtr particle) const = 0;

  /**
   * The density. Return the pdf for the given \a parton inside the
   * given \a particle for the virtuality \a partonScale and
   * logarithmic momentum fraction \a l \f$(l=\log(1/x)\f$. The \a
   * particle is assumed to have a virtuality \a particleScale.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = ZERO) const;

  /**
   * The density. Return the pdf for the given \a parton inside the
   * given \a particle for the virtuality \a partonScale and momentum
   * fraction \a x. The \a particle is assumed to have a virtuality \a
   * particleScale.
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps = 0.0,
		     Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\f$. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * implementation will assume that the difference between a quark
   * and anti-quark distribution is due do valense quarks, but return
   * zero for anything else.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and momentum fraction \a x. The \a particle is
   * assumed to have a virtuality \a particleScale. If not overidden
   * by a sub class this implementation will assume that the
   * difference between a quark and anti-quark distribution is due do
   * valense quarks, but return zero for anything else.
   */
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double x, double eps = 0.0,
		      Energy2 particleScale = ZERO) const;

  /**
   * The sea density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\f$. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * implementation will assume that the difference between a quark
   * and anti-quark distribution is due do valense quarks.
   */
  virtual double xfsl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double l, Energy2 particleScale = ZERO) const;

  /**
   * The sea density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and momentum fraction \a x. The \a particle is
   * assumed to have a virtuality \a particleScale. If not overidden
   * by a sub class this implementation will assume that the
   * difference between a quark and anti-quark distribution is due do
   * valense quarks.
   */
  virtual double xfsx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double x, double eps = 0.0,
		      Energy2 particleScale = ZERO) const;

  /**
   * Generate a momentum fraction. If the PDF contains strange peaks
   * which can be difficult to handle, this function may be
   * overwritten to return an appropriate \f$l=\log(1/x)\f$ for a \a z
   * uniformly distributed in ]0,1[. Also the jacobobian of the
   * \f$l\rightarrow z\f$ variable transformation must in the function
   * multiply the \a jacobian argument. The default version will
   * simply use the function \f$l(z) = l_{\min} +
   * z*(l_{\max}-l_{\min})\f$ (where the limits are set by \a cut).
   */
  virtual double flattenL(tcPDPtr particle, tcPDPtr parton, const PDFCuts &cut,
			  double z, double & jacobian) const;

  /**
   * Generate scale (as a fraction of the maximum scale). If the PDF
   * contains strange peaks which can be difficult to handle, this
   * function may be overwritten to return an appropriate scale
   * \f$Q^2/Q^2_{\max}\f$ for a \a z uniformly distributed in
   * ]0,1[. Also the jacobobian of the \f$Q^2/Q^2_{\max}\rightarrow
   * z\f$ variable transformation must multiply the \a jacobian
   * argument. The default version will simply use the function
   * \f$Q^2/Q^2_{\max} = (Q^2_{\max}/Q^2_{\min})^(z-1)\f$ or, if
   * \f$Q^2_{\min}\f$ is zero, \f$Q^2/Q^2_{\max} = z\f$ (where the
   * limits are set by \a cut).
   */
  virtual double flattenScale(tcPDPtr particle, tcPDPtr parton,
			       const PDFCuts & cut, double l, double z,
			       double & jacobian) const;
  //@}

  /**
   * Pointer to the remnant handler to handle remnant when extracting
   * partons according to these densities.
   */
  tcRemHPtr remnantHandler() const { return theRemnantHandler; }


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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

protected:

  /**
   * A remnant handler which can generate remnants for the parton
   * extracted withfor this PDF
   */
  RemHPtr theRemnantHandler;

protected:

  /**
   * Indicate how to deal with x and Q2 which are out of range.
   */
  enum RangeException {
    rangeFreeze, /**> Freeze the value of the PDF outside the limits. */
    rangeZero,   /**> Set the PDF to zero outside the limits. */
    rangeThrow   /**> Throw an exception if outside the limits. */
  };

  /**
   * Indicate to subclasses how to deal with x and Q2 which are out of
   * range.
   */
  RangeException rangeException;

private:


  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<PDFBase> initPDFBase;

  /**
   *  Private and non-existent assignment operator.
   */
  PDFBase & operator=(const PDFBase &);

};

ThePEG_DECLARE_CLASS_TRAITS(PDFBase,HandlerBase);

}

#endif /* ThePEG_PDFBase_H */
