// -*- C++ -*-
//
// LHAPDF.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LHAPDF_H
#define THEPEG_LHAPDF_H
//
// This is the declaration of the LHAPDF class.
//

#include "ThePEG/PDF/PDFBase.h"

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

  /** Enumerate the allowed particle types. */
  enum PType { nucleonType = 1, /**< (Anti-) proton or neutron. */
	       pionType = 2,    /**< Pion */
	       photonType = 3   /** Photon possible with anomalous component. */
  };

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
   * The particle type. 1=nucleon, 2=pion, 3=photon. No checking is
   * done to see if the selected PDF set in LHAPDF actually can handle
   * this type
   */
  PType ptype() const { return thePType; }

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
   * Read a line from the index file
   */
  bool indexLine(istream & is, int & set, int & mem, string & file,
		 int & pdftyp, int & pdfgup, int & pdfsup,
		 double & xmin, double & xmax,
		 double & q2min, double & q2max) const;

  /**
   * Call the Fortran InitPDFSetM function.
   */
  void initpdfsetm() const;

  /**
   * Call the Fortran InitPDFSetM function.
   */
  void initpdfm() const;

  /**
   * Reset the saved values from the last call to xfx(), etc.
   */
  void lastReset() const;

  /**
   * Acquire a new nset number.
   */
  void setnset() const;

  /**
   * Retrieve the number of members in the chosen PDF.
   */
  int getMaxMember() const;

  /**
   * Get the maximum number of flavours available in the chosen PDF.
   */
  int getMaxFlav() const;

  /**
   * Initialize the LHAPDF library for the chosen PDF set if it has
   * not been done before.
   */
  void checkInit() const;

  /**
   * Retrieve new PDF values for the given parameters if they were
   * changed since the last call.
   */
  void checkUpdate(double x, Energy2 Q2, Energy2 P2) const;

  /**
   * Set the maximum number of simultaneous pdfs that can be used in
   * LHAPDF. Should be set to the parameter <code>nmxset</code> in the
   * <code>parmsetup.inc</code> in the installed LHAPDF library. (By
   * default this is set to 3.)
   */
  void setMaxNSet(int);

  /**
   * Get the maximum number of simultaneous pdfs that can be used in
   * LHAPDF.
   */
  int getMaxNSet() const;

  /**
   * Deduce the min/max values of \f$x\f$ and \f$Q^2\f$ for the
   * selected set.
   */
  void setMinMax();

  /**
   * Used by the interface to select a set and member according to a number.
   */
  void setPDFNumber(int n);

  /**
   * Used by the interface to select a get the index number of the
   * currently chosen set and member.
   */
  int getPDFNumber() const;

  /**
   * Used by the interface to select a set and member according to the
   * old PDFLIB numbers.
   */
  void setPDFLIBNumbers(int group, int num);

  /**
   * Used by the interface to select a set and member according to the
   * old PDFLIB numbers.
   */
  string setPDFLIBNumbers(string);

  /**
   * Used by the interface to select a get the old PDFLIB numbers of the
   * currently chosen set and member.
   */
  pair<int,int> getPDFLIBNumbers() const;

  /**
   * Used by the interface to select a set according to a file name.
   */
  void setPDFName(string name);

  /**
   * Used by the interface to select a member in the current set.
   */
  void setPDFMember(int n);

  /**
   * Try to determine the path to where the LHAPDF index file is located.
   */
  static std::string getIndexPath();

  /**
   * Try to find a LHAPDF index file, open it in the given file stream.
   */
  static bool openLHAIndex(ifstream & is);

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
  class NotInstalled: public Exception {};

  /** Function to throw a NotInstalled exception. */
  static void throwNotInstalled();

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
   * The particle type. 1=nucleon, 2=pion, 3=photon. No checking is
   * done to see if the selected PDF set in LHAPDF actually can handle
   * this type
   */
  PType thePType;

  /**
   * The name if the PDF set to be used. Should be the full name
   * including the <code>.LHpdf</code> or <code>.LHgrid</code> suffix.
   */
  string thePDFName;

  /**
   * The chosen member of the selected PDF set.
   */
  int theMember;

  /**
   * If this is a photon PDF, this describes the option for how to
   * treat the anomalous component.
   */
  int thePhotonOption;

  /**
   * If this PDF allows partonic photons inside a hadron, enable this
   */
  bool enablePartonicGamma;

  /**
   * The verbosity of the output from the LHAPDF library.
   */
  int theVerboseLevel;

  /**
   * The maximum number of flavours for which non-zero densities are
   * reported. The actual number of flavours may be less depending on
   * the chosen PDF set.
   */
  int theMaxFlav;

  /**
   * The LHAPDF nset number (minus one) to be used by this object. If
   * below zero, the object has not been initialized.
   */
  mutable int nset;

  /**
   * Save the last \f$Q^2\f$ value used to avoid recalculation.
   */
  mutable Energy2 lastQ2;

  /**
   * Save the last \f$x\f$ value used to avoid recalculation.
   */
  mutable double lastX;

  /**
   * Save the last \f$P^2\f$ value used for off-shell photon to avoid
   * recalculation.
   */
  mutable Energy2 lastP2;

  /**
   * Save the last function values returned from the LHAPDF library.
   */
  mutable vector<double> lastXF;

  /**
   * The maximum number of simultaneous pdfs that can be used in
   * LHAPDF. Should be set to the parameter <code>nmxset</code> in the
   * <code>parmsetup.inc</code> in the installed LHAPDF library. (By
   * default this is set to 3.)
   */
  static int MaxNSet;

  /**
   * The last nset number used by an LHAPDF object.
   */
  static int lastNSet;

  /**
   * The last names used in the initialization of a given nset number.
   */
  static vector<string> lastNames;

  /**
   * The last mem used in the initialization of a given nset number.
   */
  static vector<int> lastMem;

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LHAPDF> initLHAPDF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHAPDF & operator=(const LHAPDF &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHAPDF. */
template <>
struct BaseClassTrait<LHAPDF,1> {
  /** Typedef of the first base class of LHAPDF. */
  typedef PDFBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHAPDF class and the shared object where it is defined. */
template <>
struct ClassTraits<LHAPDF>
  : public ClassTraitsBase<LHAPDF> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::LHAPDF"; }
  /** Return the name of the shared library be loaded to get access to
   *  the LeptonLeptonPDF class and every other class it uses (except
   *  the base class). */
  static string library() { return "ThePEGLHAPDF.so"; }
};

/** @endcond */

}

#endif /* THEPEG_LHAPDF_H */
