// -*- C++ -*-
#ifndef PYTHIA7_SpaceShowerHandler_H
#define PYTHIA7_SpaceShowerHandler_H
// This is the declaration of the SpaceShowerHandler class.

#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "Pythia7/Shower/SpaceShower.h"
#include "SpaceShowerHandler.fh"
// #include "SpaceShowerHandler.xh"

#include "TimeShowerHandler.fh"

namespace Pythia7 {

/**
 * SpaceShowerHandler is a wrapper around the internal
 * Shower::SpaceShower class to perform a space-like parton shower.
 *
 * @see \ref SpaceShowerHandlerInterfaces "The interfaces"
 * defined for SpaceShowerHandler.
 */
class SpaceShowerHandler: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SpaceShowerHandler();

  /**
   * The copy constructor.
   */
  inline SpaceShowerHandler(const SpaceShowerHandler &);

  /**
   * The destructor.
   */
  virtual ~SpaceShowerHandler();
  //@}

public:

  /**
   * Set the static parameters of the underlying model object and
   * return an instance.
   */
  Shower::SpaceShower * getModel();

  /**
   * Allowed ISR showers for incoming hadron: = 0: none; = 1: QCD
   * ones; = 2 : also allow photon emission; = 3 : allow spacelike
   * photon branchings, if supported by the PDF used (not the case
   * currently) (does not address VMD part of resolved photons).
   */
  inline int hadronShower() const;

  /**
   * Allowed ISR showers for incoming lepton: = 0: none; = 1: photon
   * emission; = 2 : allow spacelike photon branchings, if supported
   * by the PDF used (not the case currently) (does not address VMD
   * part of resolved photons).
   */
  inline int leptonShower() const;

  /**
   * Number of allowed quark flavours in \f$g\rightarrow q \bar{q}\f$
   * branching.
   */
  inline int nQuarks() const;

  /**
   * Running of alpha_strong in evolution: = 0: fixed; = 1: scale
   * \f$Q^2\f$; = 2: scale \f$p_\perp^2\f$.  But note that PDF's are
   * always evaluated at \f$Q^2\f$.
   */
  inline int alphaSMode() const;

  /**
   * Maximum virtuality setting the starting point of the evolution: =
   * 0: \f$s\f$; = 1: \f$\hat{s}\f$; = 2: average \f$m_\perp^2\f$; =
   * 3: smallest \f$m_\perp^2\f$.
   */
  inline int maxVirtuality() const;

  /**
   * Use of matrix element corrections:
   * = 0: no; = 1: yes.
   */
  inline int MEMode() const;

  /**
   * Resum the effect of multiple soft gluon emissions: = 0: no; = 1: yes.
   */
  inline int softGluonResum() const;

  /**
   * Restrict first emission within cone given by colour flow in hard
   * process.  = 0: no; = 1: yes, isotropic phi angle inside cone; = 2:
   * yes, also with anisotropic phi angle inside cone.
   */
  inline int finalCone() const;

  /**
   * Q2 ordering is normal, but could be relaxed (to be developed in
   * future).  = 0: yes; = 1: no, allow maximum kinematically possible
   * (toy scenario).
   */
  inline int Q2Order() const;

  /**
   * Use angular ordering?
   */
  inline bool angularOrdering() const;

  /**
   * Azimuthal asymmetry induced by gluon polarization.
   * = 0: no; = 1: yes.
   */
  inline int phiPolAsym() const;

  /**
   * Azimuthal asymmetry induced by colour coherence.
   * = 0: no; = 1: yes.
   */
  inline int phiCoherAsym() const;

  /**
   * Use the scale variable of original partons to restrict
   * branchings.  = 0: no; = 1: yes, \f$Q^2 <\f$ scale; = 2: yes,
   * \f$p_\perp^2 <\f$ scale, = 3: yes, \f$(E\theta)^2 <\f$ scale; =
   * 4: yes, \f$\theta^2 <\f$ scale.  Here theta is the full opening
   * angle of a branching, defined in the rest frame of the event. (In
   * all cases relations are approximate.)
   */
  inline int respectScale() const;

  /**
   * Parton shower cut-off mass for QCD emissions.
   */
  inline Energy Q0() const;

  /**
   * Parton shower cut-off mass for photon coupling to coloured particle.
   */
  inline Energy Q0ChgQ() const;

  /**
   * Parton shower cut-off mass for pure QED branchings. Assumed <= Q0CHGQ.
   */
  inline Energy Q0ChgL() const;

  /**
   * Fixed alpha_strong value for AlphaSMode == 0.
   */
  inline double alphaSFix() const;

  /**
   * Lambda_QCD(five flavours) in alpha_strong for AlphaSMode >= 1.
   */
  inline Energy Lambda5() const;

  /**
   * Fixed alpha_em value.
   */
  inline double alphaEMFix() const;

  /**
   * Minimum energy of emitted QCD parton in rest frame of subprocess.
   */
  inline Energy EMinEmitted() const;

  /**
   * Minimum fraction \f$1 - z\f$ of emitted QCD parton, in addition
   * to other limits.
   */
  inline double zMinEmitted() const;

  /**
   * Minimum \f$x\f$ fraction of emitted photon - matched to treatment of
   * photon PDF.
   */
  inline double xMinEmittedChg() const;

  /**
   * Smallest particle mass for QED evolution (= electron mass).
   */
  inline Energy tinyQChg() const;

  /**
   * Vanishingly small parton density.
   */
  inline double tinyPDF() const;

  /**
   * Vanishingly small product of splitting kernels and parton density
   * ratios.
   */
  inline double tinyKernelPDF() const;

  /**
   * Vanishingly small recoil mass in branching kinematics reconstruction.
   */
  inline double tinyKinPrec() const;

  /**
   * Safety margin in \f$x\f$ that heavy flavour evolution is at all
   * possible.
   */
  inline double heavyEvol() const;

  /**
   * Extra preweight in QED shower evolution, to avoid maximum violation.
   */
  inline double extraPreweight() const;

  /**
   * Maximum allowed \f$x\f$ when reconstructing back to heavy flavour
   * from gluon or photon.
   */
  inline double heavyMax() const;

  /**
   * Mimimum gap in \f$Q^2\f$ values to allow iteration when parton
   * density vanishes.
   */
  inline double Q2StartFrac() const;

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The object implementing the actual model.
   */
  Shower::SpaceShower * theShowerModel;

  /**
   * Allowed ISR showers for incoming hadron: = 0: none; = 1: QCD
   * ones; = 2 : also allow photon emission; = 3 : allow spacelike
   * photon branchings, if supported by the PDF used (not the case
   * currently) (does not address VMD part of resolved photons).
   */
  int theHadronShower;

  /**
   * Allowed ISR showers for incoming lepton: = 0: none; = 1: photon
   * emission; = 2 : allow spacelike photon branchings, if supported
   * by the PDF used (not the case currently) (does not address VMD
   * part of resolved photons).
   */
  int theLeptonShower;

  /**
   * Number of allowed quark flavours in g -> q qbar branching.
   */
  int theNQuarks;

  /**
   * Running of alpha_strong in evolution:
   * = 0: fixed; = 1: scale Q^2; = 2: scale pT^2. 
   * But note that PDF's are always evaluated at Q2.
   */
  int theAlphaSMode;

  /**
   * Maximum virtuality setting the starting point of the evolution: =
   * 0: s; = 1: sHat; = 2: average mT^2; = 3: smallest mT^2.
   */
  int theMaxVirtuality;

  /**
   * Use of matrix element corrections:
   * = 0: no; = 1: yes.
   */
  int theMEMode;

  /**
   * Resum the effect of multiple soft gluon emissions: = 0: no; = 1: yes.
   */
  int theSoftGluonResum;

  /**
   * Restrict first emission within cone given by colour flow in hard
   * process.  = 0: no; = 1: yes, isotropic phi angle inside cone; = 2:
   * yes, also with anisotropic phi angle inside cone.
   */
  int theFinalCone;

  /**
   * Q2 ordering is normal, but could be relaxed (to be developed in
   * future).  = 0: yes; = 1: no, allow maximum kinematically possible
   * (toy scenario).
   */
  int theQ2Order;

  /**
   * Use angular ordering if non-zero.
   */
  int useAngularOrdering;

  /**
   * Azimuthal asymmetry induced by gluon polarization.
   * = 0: no; = 1: yes.
   */
  int thePhiPolAsym;

  /**
   * Azimuthal asymmetry induced by colour coherence.
   * = 0: no; = 1: yes.
   */
  int thePhiCoherAsym;

  /**
   * Use the scale variable of original partons to restrict
   * branchings.  = 0: no; = 1: yes, the Q2 < scale; = 2: yes, the pT2
   * < scale, = 3: yes, the (E*theta)^2 < scale; = 4: yes, the theta^2
   * < scale.  Here theta is the full opening angle of a branching,
   * defined in the rest frame of the event. (In all cases relations
   * are approximate.)
   */
  int theRespectScale;

  /**
   * Parton shower cut-off mass for QCD emissions.
   */
  Energy theQ0;

  /**
   * Parton shower cut-off mass for photon coupling to coloured particle.
   */
  Energy theQ0ChgQ;

  /**
   * Parton shower cut-off mass for pure QED branchings. Assumed <= Q0CHGQ.
   */
  Energy theQ0ChgL;

  /**
   * Fixed alpha_strong value for AlphaSMode == 0.
   */
  double theAlphaSFix;

  /**
   * Lambda_QCD(five flavours) in alpha_strong for AlphaSMode >= 1.
   */
  Energy theLambda5;

  /**
   * Fixed alpha_em value. 
   */
  double theAlphaEMFix;

  /**
   * Minimum energy of emitted QCD parton in rest frame of subprocess.
   */
  Energy theEMinEmitted;

  /**
   * Minimum fraction 1 - z of emitted QCD parton, in addition to
   * other limits.
   */
  double theZMinEmitted;

  /**
   * Minimum x fraction of emitted photon - matched to treatment of
   * photon PDF.
   */
  double theXMinEmittedChg;

  /**
   * Smallest particle mass for QED evolution (= electron mass).
   */
  Energy theTinyQChg;

  /**
   * Vanishingly small parton density.
   */
  double theTinyPDF;

  /**
   * Vanishingly small product of splitting kernels and parton density
   * ratios.
   */
  double theTinyKernelPDF;

  /**
   * Vanishingly small recoil mass in branching kinematics reconstruction. 
   */
  double theTinyKinPrec;

  /**
   * Safety margin in x that heavy flavour evolution is at all possible.
   */
  double theHeavyEvol;

  /**
   * Extra preweight in QED shower evolution, to avoid maximum violation.
   */
  double theExtraPreweight;

  /**
   * Maximum allowed x when reconstructing back to heavy flavour from
   * gluon or photon.
   */
  double theHeavyMax;

  /**
   * Mimimum gap in Q2 values to allow iteration when parton density
   * vanishes.
   */
  double theQ2StartFrac;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SpaceShowerHandler> initSpaceShowerHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  SpaceShowerHandler & operator=(const SpaceShowerHandler &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::SpaceShowerHandler.
 */
template <>
struct BaseClassTrait<Pythia7::SpaceShowerHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::SpaceShowerHandler. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::SpaceShowerHandler class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Pythia7::SpaceShowerHandler>
  : public ClassTraitsBase<Pythia7::SpaceShowerHandler> {
  /** Return the class name. */
  static string className() { return "Pythia7::SpaceShowerHandler"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Pythia7::SpaceShowerHandler class and every other class it uses
   *  (except the base class). */
  static string library() { return "libP7Shower.so"; }

};

/** @endcond */

}

#include "SpaceShowerHandler.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "SpaceShowerHandler.tcc"
#endif

#endif /* PYTHIA7_SpaceShowerHandler_H */
