// -*- C++ -*-
#ifndef DIPSY_DipoleAnalysisHandler_H
#define DIPSY_DipoleAnalysisHandler_H
//
// This is the declaration of the DipoleAnalysisHandler class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleAnalysisHandler.fh"
#include "Ariadne/DIPSY/DipoleState.fh"
#include "Ariadne/DIPSY/ImpactParameters.h"
#include "Ariadne/DIPSY/DipoleXSec.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * The DipoleAnalysisHandler class can be used as a base class for
 * sub-classes which can be used in the initialization of a
 * DipoleEventHandler, where statistics on, eg. total and elastic
 * cross sections, can be collected during the presampling phase.
 *
 * @see \ref DipoleAnalysisHandlerInterfaces "The interfaces"
 * defined for DipoleAnalysisHandler.
 */
class DipoleAnalysisHandler: public HandlerBase {

public:

  /**
   * Convenient typedef
   */
  typedef vector< vector< vector<double> > > Vec3D;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleAnalysisHandler();

  /**
   * The destructor.
   */
  virtual ~DipoleAnalysisHandler();
  //@}

public:


  /** @name Standard virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Initialize the analysis object.
   */
  virtual void initialize() = 0;

  /**
   * Analyze a given collision given left- and right-moving dipole
   * states, \a dl and \a dr, an ImpactParameters object, \a b, a
   * cross section object, \a xsec, the total summed dipole-dipole
   * scattering probabilities, \a fsum, and the total \a weight
   * associated with the generated states.
   */
  virtual void analyze(const DipoleState & dr, const DipoleState & dl,
		       const ImpactParameters & b, const DipoleXSec & xsec,
		       double fsum, CrossSection weight);

  /**
   * Analyze a given set of collision given a set of left- and
   * right-moving dipole states in \a vl and \a vr, an a set of
   * ImpactParameters object in \a vb, and a cross section object, \a
   * xsec, and the total summed dipole-dipole scattering
   * probabilities, \a probs. Also a jacobian, \a jac, is supplied in
   * case of varying total energy.
   */
  virtual void analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
		       const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
		       const Vec3D & probs, double jac);

  /**
   * Finalize the analysis, (compute statistics etc.). \a neve is the
   * number of times analyze() has been called since last
   * initialize().
   */
  virtual void finalize(long neve) = 0;
  //@}

  /**
   * Write out stub for output line.
   */
  ostream & stub(string) const;

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleAnalysisHandler & operator=(const DipoleAnalysisHandler &);

};

}

#endif /* DIPSY_DipoleAnalysisHandler_H */
