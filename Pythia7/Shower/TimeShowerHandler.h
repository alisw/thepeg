// -*- C++ -*-
#ifndef PYTHIA7_TimeShowerHandler_H
#define PYTHIA7_TimeShowerHandler_H
// This is the declaration of the TimeShowerHandler class.

#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "Pythia7/Shower/TimeShower.h"

#include "TimeShowerHandler.fh"
// #include "TimeShowerHandler.xh"

namespace Pythia7 {

/**
 * TimeShowerHandler is a wrapper around the internal
 * Shower::TimeShower class to perform a space-like parton shower.
 *
 * @see \ref TimeShowerHandlerInterfaces "The interfaces"
 * defined for TimeShowerHandler.
 */
class TimeShowerHandler: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline TimeShowerHandler();

  /**
   * The copy constructor.
   */
  inline TimeShowerHandler(const TimeShowerHandler &);

  /**
   * The destructor.
   */
  virtual ~TimeShowerHandler();
  //@}

public:

  /**
   * Set the static parameters of the underlying model object and
   * return an instance.
   */
  Shower::TimeShower * getModel();

  /**
   * Angular ordering; = 0: off, = 1: on for g emission, = 2: on also
   * for \f$g\rightarrow q \bar{q}\f$ splitting.
   */
  inline int angularOrder() const;

  /**
   * Number of allowed quark flavours in \f$g\rightarrow q \bar{q}\f$
   * branching.
   */
  inline int nQuark() const;

  /**
   * Running of alpha_strong in evolution: = 0: fixed; = 1: scale
   * \f$Q^2/4\f$; = 2: scale \f$p_\perp^2\f$; = 3: scale
   * \f$p_\perp^2\f$, except for \f$g\rightarrow q \bar{q}\f$, where
   * it is \f$Q^2/4\f$.
   */
  inline int alphaSMode() const;

  /**
   * Also allow a QED shower together with QCD ones: = 0: no; = 1:
   * radiate on-shell photons; = 2 : also allow photons to branch.
   */
  inline int MEMode() const;

  /**
   * Also allow a QED shower together with QCD ones: = 0: no;
   * = 1: radiate on-shell photons; = 2 : also allow photons to branch.
   */
  inline int QEDShower() const;

  /**
   * Restrict first emission within cone given by colour flow in hard
   * process.  = 0: no; = 1: yes, isotropic phi angle inside cone; =
   * 2: yes, also with anisotropic phi angle inside cone.
   */
  inline int initialCone() const;

  /**
   * Azimuthal asymmetry induced by gluon polarization.
   * = 0: no; = 1: yes.
   */
  inline int phiPolAsym() const;

  /**
   * Azimuthal asymmetry induced by colour coherence.
   * = 0: no; = 1: yes.
   */
  inline int phoCoherAsym() const;

  /**
   * Use the scale variable of original partons to restrict
   * branchings.  = 0: no; = 1: yes, \f$Q^2 <\f$ scale; = 2: yes,
   * \f$p_\perp^2 <\f$ scale, = 3: yes, \f$(E\theta)^2 <\f$ scale; =
   * 4: yes, \f$\theta^2 <\f$ scale.  (In all cases relations are
   * approximate.)
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
   * \f$\Lambda_{\mbox{QCD}}\f$ (five flavours) in alpha_strong for
   * AlphaSMode >= 1.
   */
  inline Energy Lambda5() const;

  /**
   * Fixed \f$\alpha_{\mbox{EM}}\f$ value.
   */
  inline double alphaEMFix() const;

  /**
   * Fraction of Q0 cut-off mass used as safety margin in 
   * daughter mass sum. Relevant for total parton multiplicity.
   */
  inline double Q0FracPS() const;

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
  Shower::TimeShower * theShowerModel;

  /**
   * Angular ordering; = 0: off, = 1: on for g emission, = 2: on also
   * for g -> q qbar splitting.
   */
  int theAngularOrdering; 


  /**
   * Number of allowed quark flavours in g -> q qbar branching.
   */
  int theNQuark; 


  /**
   * Running of alpha_strong in evolution: = 0: fixed; = 1: scale
   * Q^2/4; = 2: scale pT^2; = 3: scale pT2, except for g -> q qbar,
   * where it is Q^2/4 .
   */
  int theAlphaSMode; 


  /**
   * Also allow a QED shower together with QCD ones: = 0: no; = 1:
   * radiate on-shell photons; = 2 : also allow photons to branch.
   */
  int theMEMode; 


  /**
   * Also allow a QED shower together with QCD ones: = 0: no;
   * = 1: radiate on-shell photons; = 2 : also allow photons to branch.
   */
  int theQEDShower; 


  /**
   * Restrict first emission within cone given by colour flow in hard
   * process.  = 0: no; = 1: yes, isotropic phi angle inside cone; =
   * 2: yes, also with anisotropic phi angle inside cone.
   */
  int theInitialCone; 


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
   * < scale.  (In all cases relations are approximate.)
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
   * Fixed alpha_EM value.
   */
  double theAlphaEMFix; 


  /**
   * Fraction of Q0 cut-off mass used as safety margin in 
   * daughter mass sum. Relevant for total parton multiplicity.
   */
  double theQ0FracPS; 


private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<TimeShowerHandler> initTimeShowerHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  TimeShowerHandler & operator=(const TimeShowerHandler &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::TimeShowerHandler.
 */
template <>
struct BaseClassTrait<Pythia7::TimeShowerHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::TimeShowerHandler. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::TimeShowerHandler class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Pythia7::TimeShowerHandler>
  : public ClassTraitsBase<Pythia7::TimeShowerHandler> {
  /** Return the class name. */
  static string className() { return "Pythia7::TimeShowerHandler"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Pythia7::TimeShowerHandler class and every other class it uses
   *  (except the base class). */
  static string library() { return "libP7Shower.so"; }

};

/** @endcond */

}

#include "TimeShowerHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TimeShowerHandler.tcc"
#endif

#endif /* PYTHIA7_TimeShowerHandler_H */
