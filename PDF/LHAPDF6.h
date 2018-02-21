// -*- C++ -*-
//
// LHAPDF6.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2014 Leif Lonnblad, David Grellscheid
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LHAPDF6_H
#define THEPEG_LHAPDF6_H
//
// This is the declaration of the LHAPDF class.
//

#include "ThePEG/PDF/PDFBase.h"

namespace LHAPDF {
	class PDF;
}

namespace ThePEG {

/**
 * The LHAPDF class inherits from PDFBase and implements an interface
 * to the LHAPDF library of parton density function
 * parameterizations. This class is available even if LHAPDF was not
 * properly installed when ThePEG was installed, but will then produce
 * an error in the initialization.
 *
 * Note that the valence densities from the xfvx() and xfvl() function
 * will only work properly for nucleons. All other particles will have
 * zero valence densities.
 *
 * @see \ref LHAPDFInterfaces "The interfaces"
 * defined for LHAPDF.
 */
class LHAPDF: public PDFBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  LHAPDF();

  /**
   * The copy constructor.
   */
  LHAPDF(const LHAPDF &);
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of partons from
   * the given \a particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the partons which this PDF may extract from the given
   * \a particle.
   */
  virtual cPDVector partons(tcPDPtr particle) const;

  /**
   * The density. Return the pdf for the given \a parton inside the
   * given \a particle for the virtuality \a partonScale and momentum
   * fraction \a x (with x = 1-\a eps). The \a particle is assumed to
   * have a virtuality \a particleScale.
   */
  virtual double xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double eps = 0.0,
		     Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l. The \a
   * particle is assumed to have a virtuality \a particleScale. This
   * will only work properly for nucleons. All other particles will
   * have zero valence densities
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double l, Energy2 particleScale = ZERO) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and momentum fraction \a x (with x = 1-\a eps). The
   * \a particle is assumed to have a virtuality \a
   * particleScale. This will only work properly for nucleons. All
   * other particles will have zero valence densities
   */
  virtual double xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		      double x, double eps = 0.0,
		      Energy2 particleScale = ZERO) const;

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
  //@}


  /** @name Simple access function. */
  //@{
  /**
   * The name if the PDF set to be used. The full name including the
   * <code>.LHpdf</code> or <code>.LHgrid</code> suffix.
   */
  const string & PDFName() const { return thePDFName; }
 
  /**
   * The chosen member of the selected PDF set.
   */
  int member() const { return theMember; }

  /**
   * The maximum number of flavours for which non-zero densities are
   * reported. The actual number of flavours may be less depending on
   * the chosen PDF set.
   */
  int maxFlav() const { return theMaxFlav; }
  //@}

protected:

  /** @name Internal helper functions. */
  //@{
  /**
   * Initialize the LHAPDF library for the chosen PDF set if it has
   * not been done before.
   */
  void initPDFptr();

  /**
   * Used by the interface to select a set according to a file name.
   */
  void setPDFName(string name);

  /**
   * Used by the interface to select a member in the current set.
   */
  void setPDFMember(int n);

  /**
   * Interface for simple tests.
   */
  string doTest(string input);
  //@}

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

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used if the LHAPDF library was not installed. */
  class NotInstalled: public InterfaceException {};

  /** @endcond */

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * LHAPDF member object
   */
  ::LHAPDF::PDF * thePDF;

  /**
   * The name of the selected PDF set.
   */
  string thePDFName;

  /**
   * The chosen member of the selected PDF set.
   */
  int theMember;

  /**
   * The maximum number of flavours for which non-zero densities are
   * reported. The actual number of flavours may be less depending on
   * the chosen PDF set.
   */
  int theMaxFlav;

  /**
   * The minimum \f$x\f$-value for the current PDF set.
   */
  double xMin;

  /**
   * The maximum \f$x\f$-value for the current PDF set.
   */
  double xMax;

  /**
   * The minimum \f$Q^2\f$-value for the current PDF set.
   */
  Energy2 Q2Min;

  /**
   * The maximum \f$Q^2\f$-value for the current PDF set.
   */
  Energy2 Q2Max;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHAPDF & operator=(const LHAPDF &);

};

}

#endif /* THEPEG_LHAPDF6_H */
