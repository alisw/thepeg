// -*- C++ -*-
#ifndef DIPSY_ElasticXSecAnalysis_H
#define DIPSY_ElasticXSecAnalysis_H
//
// This is the declaration of the ElasticXSecAnalysis class.
//

#include "DipoleAnalysisHandler.h"

#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the ElasticXSecAnalysis class.
 *
 * @see \ref ElasticXSecAnalysisInterfaces "The interfaces"
 * defined for ElasticXSecAnalysis.
 */
class ElasticXSecAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ElasticXSecAnalysis();

  /**
   * The destructor.
   */
  virtual ~ElasticXSecAnalysis();
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
   * Number of bins in impact parameter.
   */
  int nb;

  /**
   * Size of interval in impact parameter.
   */
  InvEnergy db;

  /**
   * The weighted sum of the interaction probability, binned in b.
   */
  vector<CrossSection> sumTw;

  /**
   * The weighted sum of the square of
   * the interaction probability, binned in b.
   */
  vector<CrossSection> sumT2w;

  /**
   * The weighted sum of the 4:th power of
   * the interaction probability, binned in b.
   */
  vector<CrossSection> sumT4w;

  /**
   * The weighted sum of the elastic interaction amplitude, binned in b.
   */
  vector<CrossSection> sumElAw;

  /**
   * The weighted sum of the square of
   * the elastic interaction amplitude, binned in b.
   */
  vector<CrossSection> sumElA2w;

  /**
   * The sum of the real weights, binned in b.
   */
  vector<CrossSection> sumw;

  /**
   * The sum of the real elastic weights, binned in b.
   */
  vector<CrossSection> sumElw;

  /**
   * The sum of the number of collisions seen so far.
   */
  vector<long> sumn;

  /**
   * Number of bins in q.
   */
  int nq;

  /**
   * Size of interval in q.
   */
  Energy dq;

  /**
   * The weighted sum of the differential elastic amplitude in q.
   */
  vector<CrossSection> sumDEw;

  /**
   * The weighted sum of the square of the differential elastic amplitude in q.
   */
  vector<CrossSection> sumDE2w;

  /**
   * The differential elastic cross section in q.
   */
  FactoryBase::tH1DPtr dSigmadq;

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

  map<InvEnergy, pair<double,CrossSection> > bmap;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ElasticXSecAnalysis & operator=(const ElasticXSecAnalysis &);

};

}

#endif /* DIPSY_ElasticXSecAnalysis_H */
