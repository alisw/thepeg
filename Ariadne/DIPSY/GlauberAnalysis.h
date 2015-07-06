// -*- C++ -*-
#ifndef DIPSY_GlauberAnalysis_H
#define DIPSY_GlauberAnalysis_H
//
// This is the declaration of the GlauberAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#include <valarray>

namespace DIPSY {

using namespace ThePEG;
using std::valarray;


/**
 * Here is the documentation of the GlauberAnalysis class.
 *
 * @see \ref GlauberAnalysisInterfaces "The interfaces"
 * defined for GlauberAnalysis.
 */
class GlauberAnalysis: public DIPSY::DipoleAnalysisHandler {

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
  static const int nstr = 8;

  /**
   * The sums.
   */
  valarray<double> sigtot, signd, sigel, sigdt, sigdr, sigdl, sigdd;
  valarray<double> sumlr, sum2lr, sumlr2, suml2r, sumr2l;

  /**
   * The sum of squares.
   */
  valarray<double> sigtot2, signd2, sigel2, sigdt2, sigdr2, sigdl2, sigdd2;
  valarray<double> sTwLR, sTw2LR, sTwLR2, sTwL2R, sTwR2L;
  valarray<double> sT2wLR, sT2w2LR, sT2wLR2, sT2wL2R, sT2wR2L;
  double swLR, sw2LR, swLR2, swL2R, swR2L;

  /**
   * The number of b-values.
   */
  long ntot;

  /**
   * The number of t-bins for the elastic cross section
   */
  int Nt;

  /**
   * The maximum t-value for the elastic cross section.
   */
  Energy2 tMax;

  /**
   * The t-values, first and second strategy
   */
  vector<FactoryBase::tH1DPtr> hists;

  /**
   * Calulate all transition probabilities
   */
  valarray<double> getT(const DipoleState & dr, const DipoleState & dl,
			const ImpactParameters & b, const DipoleXSec & xsec,
			double fsum) const;

  /**
   * helper function for printing
   */
  void print(valarray<double> sig, valarray<double> sig2,
	     int ntot, string xstype) const;
  /**
   * helper function for histogram booking;
   */
  void bookHistos();

  /**
   * helper function for printing
   */
  string getStrat(int i) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GlauberAnalysis & operator=(const GlauberAnalysis &);

};

}

#endif /* DIPSY_GlauberAnalysis_H */
