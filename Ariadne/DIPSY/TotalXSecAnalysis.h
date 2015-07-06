// -*- C++ -*-
#ifndef DIPSY_TotalXSecAnalysis_H
#define DIPSY_TotalXSecAnalysis_H
//
// This is the declaration of the TotalXSecAnalysis class.
//

#include "Ariadne/DIPSY/DipoleAnalysisHandler.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * TotalXSecAnalysis class calculates the total cross section in the
 * presampling phase of the DipoleEventHandler.
 *
 * @see \ref TotalXSecAnalysisInterfaces "The interfaces"
 * defined for TotalXSecAnalysis.
 */
class TotalXSecAnalysis: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TotalXSecAnalysis();

  /**
   * The destructor.
   */
  virtual ~TotalXSecAnalysis();
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
   * The sum of the cross sections seen so far.
   */
  CrossSection sum;

  /**
   * The sum of the squared cross sections seen so far.
   */
  QTY<4,0,0>::Type sum2;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TotalXSecAnalysis & operator=(const TotalXSecAnalysis &);

};

}

#endif /* DIPSY_TotalXSecAnalysis_H */
