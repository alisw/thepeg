// -*- C++ -*-
#ifndef DIPSY_NoFluctXSecAnalysis_H
#define DIPSY_NoFluctXSecAnalysis_H
//
// This is the declaration of the NoFluctXSecAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the NoFluctXSecAnalysis class.
 *
 * @see \ref NoFluctXSecAnalysisInterfaces "The interfaces"
 * defined for NoFluctXSecAnalysis.
 */
class NoFluctXSecAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NoFluctXSecAnalysis();

  /**
   * The destructor.
   */
  virtual ~NoFluctXSecAnalysis();
  //@}

public:

  /** @name Standard virtual functions inherited from the classes. */
  //@{
  /**
   * Initialize the analysis object.
   */
  virtual void initialize();

  /**
   * Analyze a given set of collision given a set of left- and
   * right-moving dipole states in \a vl and \a vr, an a set of
   * ImpactParameters object in \a vb,  and a cross section object, \a
   * xsec. Also a jacobian, \a jac, is supplied in case of varying
   * total energy.
   */
  virtual void analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
		       const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
		       const Vec3D & probs, double jac);


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
   * The sum over interactions.
   */
  CrossSection sum00;

  /**
   * The sum over interactions with left-going cascades averaged.
   */
  CrossSection sumL0, sum2L0;

  /**
   * The sum over interactions with right-going cascades averaged.
   */
  CrossSection sum0R;


  /**
   * The sum over squared interactions with both cascades averaged.
   */
  CrossSection sumLR;

  /**
   * Convenient typedef.
   */
  typedef QTY<4,0,0>::Type CrossSection2;

  /**
   * The sum over squared interactions.
   */
  CrossSection2 sum002;

  /**
   * The sum over squared interactions with left-going cascades averaged.
   */
  CrossSection2 sumL02, sum2L02;

  /**
   * The sum over squared interactions with right-going cascades averaged.
   */
  CrossSection2 sum0R2;


  /**
   * The sum over squared interactions with both cascades averaged.
   */
  CrossSection2 sumLR2;

  /**
   * The number of interactions.
   */
  long n00;

  /**
   * The number interactions with left-going cascades averaged.
   */
  long nL0, n2L0;

  /**
   * The bumber interactions with right-going cascades averaged.
   */
  long n0R;


  /**
   * The number interactions with both cascades averaged.
   */
  long nLR;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NoFluctXSecAnalysis & operator=(const NoFluctXSecAnalysis &);

};

}

#endif /* DIPSY_NoFluctXSecAnalysis_H */
