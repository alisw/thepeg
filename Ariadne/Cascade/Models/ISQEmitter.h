// -*- C++ -*-
#ifndef ARIADNE5_ISQEmitter_H
#define ARIADNE5_ISQEmitter_H
//
// This is the declaration of the ISQEmitter class.
//

#include "Ariadne/Cascade/EmitterBase.h"
#include "Ariadne/Cascade/RemnantParton.h"
#include "ISQEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ISQEmitter class implements common functions for initial-state
 * quark emissions.
 *
 * @see \ref ISQEmitterInterfaces "The interfaces"
 * defined for ISQEmitter.
 */
class ISQEmitter: public EmitterBase {

public:

  /**
   * Convenient typedef
   */
  typedef vector< pair<Energy2,double> > PDFLimits;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ISQEmitter();

  /**
   * The destructor.
   */
  virtual ~ISQEmitter();
  //@}

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * If the given EmissionModel overlaps with this model for the given
   * Emitter, return true if this model should take precedence. Must
   * only be called for Emitters for which canHandle() is true.
   */
  virtual bool overrides(const EmitterBase &, DipoleBase &) const;
  //@}

public:

  /**
   * Called by perform in sub-classes to set momenta.
   */
  virtual void setMom(QCDDipole & d, const ISQEmission & e,
		      tRemParPtr rem, tcPDPtr qex) const;

  /**
   * Called by revert in sub-classes to reset momenta.
   */
  virtual void revertMom(QCDDipole & d, const ISQEmission & e, tRemParPtr rem) const;

  /**
   * Helper function to define maximum PDF ratios in suitable mt2
   * intervals.
   */
  PDFLimits
  maxPDFRatios(tRemParPtr rem, Energy2 mt2max, Energy2 mt2min, tcPDPtr g,
	       Energy2 s, Energy2 mh2, Energy2 mq2, Energy2 Q2 = ZERO) const;

public:

  /** Exception class in case there is not enough energy to emit any quark. */
  struct KinematicsException: public Exception {};

  /** Exception class in case the PDF ratio wwas not properly overestimated. */
  struct WeightException: public Exception {};

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * Maximum x-value above which PDF-ratios above the maximum is
   * silently ignored.
   */
  double maxPDFX;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ISQEmitter & operator=(const ISQEmitter &);

};

}

#endif /* ARIADNE5_ISQEmitter_H */
