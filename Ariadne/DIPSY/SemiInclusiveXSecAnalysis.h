// -*- C++ -*-
#ifndef DIPSY_SemiInclusiveXSecAnalysis_H
#define DIPSY_SemiInclusiveXSecAnalysis_H
//
// This is the declaration of the SemiInclusiveXSecAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the SemiInclusiveXSecAnalysis class.
 *
 * @see \ref SemiInclusiveXSecAnalysisInterfaces "The interfaces"
 * defined for SemiInclusiveXSecAnalysis.
 */
class SemiInclusiveXSecAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SemiInclusiveXSecAnalysis();

  /**
   * The destructor.
   */
  virtual ~SemiInclusiveXSecAnalysis();
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
   * Convenient typedef.
   */
  typedef QTY<4,0,0>::Type CrossSection2;

  /**
   * The sums.
   */
  CrossSection sigtot, signd, sigel, sigdt, sigdr, sigdl, sigdd;

  /**
   * The sum of squares.
   */
  CrossSection2 sigtot2, signd2, sigel2, sigdt2, sigdr2, sigdl2, sigdd2;

  /**
   * The number of b-values.
   */
  long ntot;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SemiInclusiveXSecAnalysis & operator=(const SemiInclusiveXSecAnalysis &);

};

}

#endif /* DIPSY_SemiInclusiveXSecAnalysis_H */
