// -*- C++ -*-
#ifndef DIPSY_GlauberAnalysis_H
#define DIPSY_GlauberAnalysis_H
//
// This is the declaration of the GlauberAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the GlauberAnalysis class.
 *
 * @see \ref GlauberAnalysisInterfaces "The interfaces"
 * defined for GlauberAnalysis.
 */
class GlauberAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /**
   * A convenient typedef.
   */
  typedef QTY<4,0,0>::Type CrossSection2;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GlauberAnalysis();

  /**
   * The destructor.
   */
  virtual ~GlauberAnalysis();
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
   * The assumed total nucleon - nucleon cross section
   */
  CrossSection nnSigTot;

  /**
   * The assumed elastic nucleon - nucleon cross section
   */
  CrossSection nnSigEl;

  /**
   * The assumed inelastic, non-diffracitve nucleon - nucleon cross section
   */
  CrossSection nnSigInND;

  /**
   * Number of strategies.
   */
  static const int nstr = 7;

  /**
   * The total number of points
   */
  double sumNTot;

  /**
   * The total number of points with collisions
   */
  vector<double> sumNColl;

  /**
   * The sum of individual weights
   */
  vector<CrossSection> sumWeights, sumsig, sumsigel, sumsignd, sumsigdt,
  sumcoll, sumcoll2, sumpart, sumpart2;

  /**
   * The sum of individual squared weights
   */
  vector<CrossSection2> sumsig2, sumsigel2, sumsignd2, sumsigdt2;

  /**
   * The sum of the cross sections seen so far.
   */
  CrossSection sum, sumND;

  /**
   * The sum of the squared cross sections seen so far.
   */
  CrossSection2 sum2, sumND2;

  /**
   * Helper function for printing.
   */
  void printstrat(int) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GlauberAnalysis & operator=(const GlauberAnalysis &);

};

}

#endif /* DIPSY_GlauberAnalysis_H */
