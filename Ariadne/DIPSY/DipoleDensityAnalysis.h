// -*- C++ -*-
#ifndef DIPSY_DipoleDensityAnalysis_H
#define DIPSY_DipoleDensityAnalysis_H
//
// This is the declaration of the DipoleDensityAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the DipoleDensityAnalysis class.
 *
 * @see \ref DipoleDensityAnalysisInterfaces "The interfaces"
 * defined for DipoleDensityAnalysis.
 */
class DipoleDensityAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleDensityAnalysis();

  /**
   * The destructor.
   */
  virtual ~DipoleDensityAnalysis();
  //@}

public:

  /** @name Standard virtual functions inherited from the classes. */
  //@{
  /**
   * Initialize the analysis object.
   */
  virtual void initialize();

  /**
   * Analyze a given collision. Given left- and right-moving dipole
   * states, \a dl and \a dr, an ImpactParameters object, \a b, a
   * cross section object, \a xsec, and the total summed dipole-dipole
   * scattering probabilities, \a fsum, and the total \a weight
   * associated with the generated states.
   */
  virtual void analyze(const DipoleState & dl, const DipoleState & dr,
		       const ImpactParameters & b, const DipoleXSec & xsec,
		       double fsum, CrossSection weight);

  /**
   * Finalize the analysis, (compute statistics etc.). \a neve is the
   * number of times analyze() has been called since last
   * initialize().
   */
  virtual void finalize(long neve);
  //@}

  /**
   * Fill dipoles from a state.
   */
  void fill(const DipoleState & ds);

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

private:

  /**
   * Sum of weights.
   */
  double sumw;

  /**
   * The bare distribution in dipole sizes.
   */
  FactoryBase::tH1DPtr dipsizes;

  /**
   * The bare distribution in gluon transverse momenta.
   */
  FactoryBase::tH1DPtr gluonpts;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleDensityAnalysis & operator=(const DipoleDensityAnalysis &);

};

}

#endif /* DIPSY_DipoleDensityAnalysis_H */
