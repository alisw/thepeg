// -*- C++ -*-
//
// Cuts.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Cuts_H
#define THEPEG_Cuts_H
//
// This is the declaration of the Cuts class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Cuts.fh"
#include "OneCutBase.h"
#include "TwoCutBase.h"
#include "MultiCutBase.h"
#include "JetFinder.h"
#include "FuzzyTheta.h"

namespace ThePEG {

/**
 * Cuts is a class for implementing kinematical cuts in ThePEG. The
 * class itself only implements cuts on the total momentum of the hard
 * sub-process, implemented as minimum and maximum values of \f$x_1\f$
 * and \f$x_2\f$ (or \f$\hat{s}=x_1x_2S_{tot}\f$ and
 * \f$\hat{y}=\log(x_1/x_2)/2\f$. Further cuts can be implemented
 * either by inheriting from this base class, in which the virtual
 * cut() function should be overridden, or by assigning objects of
 * class OneCutBase, TwoCutBase and MultiCutBase defining cuts on
 * single particles, pairs of particles and groups of particles in the
 * hard sub-process respectively.
 *
 * The Cuts object must be initialized specifying the overall
 * laboratory frame, giving the total squared invariant mass, \f$S\f$,
 * and the rapidity, \f$Y\f$, of the colliding particles in this
 * frame. The colliding particles are thus assumed to be directed
 * along the \f$z\f$-axis.
 *
 * For each event, the Cuts object must also be initialized giving the
 * squared invarint mass, \f$\hat{s}\f$, and the total rapidity,
 * \f$\hat{y}\f$, of the hard sub-process in the center-of-mass frame
 * of the colliding particles. Note that this means that the
 * transformation between the lab frame and the rest frame of the hard
 * sub-process is assumed to be a simple boost along the z-axis.
 *
 * @see \ref CutsInterfaces "The interfaces"
 * defined for Cuts.
 */
class Cuts: public Interfaced {

public:

  /**
   * A vector of OneCutBase pointers.
   */
  typedef vector<OneCutPtr> OneCutVector;

  /**
   * A vector of TwoCutBase pointers.
   */
  typedef vector<TwoCutPtr> TwoCutVector;

  /**
   * A vector of MultiCutBase pointers.
   */
  typedef vector<MultiCutPtr> MultiCutVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Cuts(Energy MhatMin=2*GeV);

  /**
   * The destructor.
   */
  virtual ~Cuts();
  //@}

public:

  /** @name Initialization functions. */
  //@{
  /**
   * Initialize this object specifying the maximum total invariant
   * mass squared, \a smax, and the total rapidity, \a Y, of the
   * colliding particles (for the maximum invariant mass). A sub-class
   * overriding this function must make sure the base-class function
   * is called. This function should be called once in the beginning
   * of a run.
   */
  virtual void initialize(Energy2 smax, double Y);

  /**
   * Initialize this object for a new event. A sub-class overriding
   * this function must make sure the base-class function is called.
   * This function is called before the generation of a new
   * sub-process, before the incoming partons have been generated.
   */
  virtual void initEvent();

  /**
   * Set information about the invariant mass squared, \a shat, and
   * rapidity, \a yhat, of the hard sub-process. The rapidity should
   * be given wrt. the center of mass of the colliding particles. A
   * sub-class overriding this function must make sure the base-class
   * function is called. This function is called before the generation
   * of a new sub-process, after the incoming partons have been
   * generated. If \a mirror is true any questions regarding cuts on
   * the sub-process in the functions minYStar(tcPDPtr),
   * maxYStar(tcPDPtr p), passCuts(const tcPDVector &, const
   * vector<LorentzMomentum> &, tcPDPtr, tcPDPtr) and passCuts(const
   * tcPVector &, tcPDPtr t1, tcPDPtr) will assume that the z-axis is
   * reversed in the sub-process rest frame. Returns false if the
   * given values were outside of the cuts.
   */
  virtual bool
  initSubProcess(Energy2 shat, double yhat, bool mirror = false) const;
  //@}

  /** @name Check functions to see if a state has passed the cuts or not. */
  //@{
  /**
   * Check if the outgoing particles, with the given types and
   * momenta, from a sub-process passes the cuts. The particles must
   * be given in the rest frame of tha hard sub-process, and the
   * initSubProcess must have been called before. Also the types of
   * the incoming partons, \a t1 and \a t2, may be given if availible.
   */
  virtual bool passCuts(const tcPDVector & ptype, const vector<LorentzMomentum> & p,
			tcPDPtr t1 = tcPDPtr(), tcPDPtr t2 = tcPDPtr()) const;

  /**
   * Check if the outgoing particles from a sub-process passes the
   * cuts. The particles must be given in the rest frame of tha hard
   * sub-process, and the initSubProcess must have been called
   * before. Also the types of the incoming partons, \a t1 and \a t2,
   * may be given if availible.
   */
  bool passCuts(const tcPVector & p,
		tcPDPtr t1 = tcPDPtr(), tcPDPtr t2 = tcPDPtr()) const;

  /**
   * Check if the incoming and outgoing particles in the given
   * sub-process passes the cuts. The sub-process must be given in its
   * rest frame, and the initSubProcess must have been called before.
   */
  bool passCuts(const SubProcess & sub) const;

  /**
   * Check if the given collision passes the cuts. The collision must
   * be given in its rest frame.
   */
  bool passCuts(const Collision & coll) const;
  //@}

  /** @name Access to cuts of the underlying cut objects. */
  //@{
  /**
   * Return the minimum allowed squared invariant mass of two outgoing
   * partons of type \a pi and \a pj. This function first determines
   * the minimum from the corresponding function from in TwoCutBase
   * objects. If no minimum was found, one is derived from
   * minKTClus(), minDurham(), minKT() and minDeltaR(), if possible.
   */
  Energy2 minSij(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the negative of the squared
   * invariant mass of an incoming parton of type \a pi and an
   * outgoing parton of type \a po. This function first determines the
   * minimum from the corresponding function from in TwoCutBase
   * objects. If no minimum was found, one is derived from minKT(), if
   * possible.
   */
  Energy2 minTij(tcPDPtr pi, tcPDPtr po) const;

  /**
   * Return the minimum allowed value of \f$\Delta
   * R_{ij}=\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ of two
   * outgoing partons of type \a pi and \a pj. Simply returns the
   * maximum of the results from calling the corresponding function in
   * the TwoCutBase objects.
   */
  double minDeltaR(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the longitudinally invariant
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$\min(p_{\perp i}, p_{\perp
   * j})\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ for two outgoing
   * partons, or simply \f$p_{\perp i}\f$ or \f$p_{\perp j}\f$ for a
   * single outgoing parton. Returns 0 if both partons are incoming. A
   * null pointer indicates an incoming parton, hence the type of the
   * incoming parton is irrelevant. Simply returns the maximum of the
   * results from calling the corresponding function in the TwoCutBase
   * objects.
   */
  Energy minKTClus(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the Durham
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$2\min(E_j^2, E_j^2)(1-\cos\theta_{ij})/\hat{s}\f$ for two
   * outgoing partons. Simply returns the maximum of the results from
   * calling the corresponding function in the TwoCutBase objects.
   */
  double minDurham(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the transverse momentum of an
   * outgoing parton. This function first determines the minimum from
   * the corresponding function from in OneCutBase objects. If no
   * minimum was found, one is derived from minKTClus(), if possible.
   */
  Energy minKT(tcPDPtr p) const;

  /**
   * Return the minimum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system. Simply returns the maximum of the results from calling
   * the corresponding function in the OneCutBase objects.
   */
  double minEta(tcPDPtr p) const;

  /**
   * Return the maximum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system. Simply returns the minimum of the results from calling
   * the corresponding function in the OneCutBase objects.
   */
  double maxEta(tcPDPtr p) const;

  /**
   * Return the minimum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system. Simply returns the maximum of the results from calling
   * the corresponding function in the OneCutBase objects.
   */
  double minRapidityMax(tcPDPtr p) const;

  /**
   * Return the maximum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system. Simply returns the minimum of the results from calling
   * the corresponding function in the OneCutBase objects.
   */
  double maxRapidityMin(tcPDPtr p) const;

  /**
   * Return the minimum allowed rapidity of an outgoing parton of the
   * given type in the center-of-mass system of the hard sub-process.
   * Only available after initSubProcess() has been called.
   */
  double minYStar(tcPDPtr p) const;

  /**
   * Return the minimum allowed rapidity of an outgoing parton of the
   * given type in the center-of-mass system of the hard sub-process.
   * Only available after initSubProcess() has been called.
   */
  double maxYStar(tcPDPtr p) const;

  /**
   * Return the minimum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for
   * efficiency. Simply returns the maximum of the results from
   * calling the corresponding function in the MultiCutBase objects.
   */
  Energy2 minS(const tcPDVector & pv) const;

  /**
   * Return the maximum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for
   * efficiency. Simply returns the minimum of the results from
   * calling the corresponding function in the MultiCutBase objects.
   */
  Energy2 maxS(const tcPDVector & pv) const;
  //@}

  /** @name Direct access to underlying cut objects. */
  //@{
  /**
   * Return a vector of pointers to objects of the given class (with
   * base class OneCutBase).
   */
  template <typename T>
  vector<typename Ptr<T>::transient_const_pointer>
  oneCutObjects() const;

  /**
   * Return a vector of pointers to objects of the given class (with
   * base class TwoCutBase).
   */
  template <typename T>
  vector<typename Ptr<T>::transient_const_pointer>
  twoCutObjects() const;

  /**
   * Return a vector of pointers to objects of the given class (with
   * base class MultiCutBase).
   */
  template <typename T>
  vector<typename Ptr<T>::transient_const_pointer>
  multiCutObjects() const;

  /**
   * Return the objects defining cuts on single outgoing partons from the
   * hard sub-process.
   */
  const OneCutVector& oneCuts() const { return theOneCuts; }

  /**
   * Return the objects defining cuts on pairs of particles in the hard
   * sub-process.
   */
  const TwoCutVector& twoCuts() const { return theTwoCuts; }

  /**
   * Return the objects defining cuts on sets of outgoing particles from the
   * hard sub-process.
   */
  const MultiCutVector& multiCuts() const { return theMultiCuts; }

  /**
   * Return the jet finder
   */
  Ptr<JetFinder>::tptr jetFinder() const { return theJetFinder; }

  /**
   * Add a OneCutBase object.
   */
  void add(tOneCutPtr c) { theOneCuts.push_back(c); }

  /**
   * Add a TwoCutBase object.
   */
  void add(tTwoCutPtr c) { theTwoCuts.push_back(c); }

  /**
   * Add a MultiCutBase object.
   */
  void add(tMultiCutPtr c) { theMultiCuts.push_back(c); }
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * The maximum allowed total invariant mass squared allowed for
   * events to be considered.
   */
  Energy2 SMax() const { return theSMax; }


  /**
   * The total rapidity of the colliding particles corresponding to
   * the maximum invariant mass squared, SMax().
   */
  double Y() const { return theY; }

  /**
   * The invariant mass squared of the hard sub-process of the event
   * being considered.
   */
  Energy2 currentSHat() const { return theCurrentSHat; }

  /**
   * The total rapidity of hard sub-process (wrt. the rest system of
   * the colliding particles so that currentYHat() + Y() gives the
   * true rapidity) of the event being considered.
   */
  double currentYHat() const { return theCurrentYHat; }

  //@}

  /** @name Functions to inquire about specific cuts. */
  //@{
  /**
   * The minimum allowed value of \f$\hat{s}\f$.
   */
  Energy2 sHatMin() const { return max(sqr(theMHatMin), theX1Min*theX2Min*SMax()); }

  /**
   * The maximum allowed value of \f$\hat{s}\f$.
   */
  Energy2 sHatMax() const { return min(sqr(theMHatMax), theX1Max*theX2Max*SMax()); }

  /**
   * Check if the given \f$\hat{s}\f$ is within the cuts.
   */
  bool sHat(Energy2 sh) const { 
    return sh > sHatMin() && sh <= sHatMax()*(1.0 + 1000.0*Constants::epsilon);
  }

  /**
   * The minimum allowed value of \f$\sqrt{\hat{s}}\f$.
   */
  Energy mHatMin() const { return max(theMHatMin, sqrt(theX1Min*theX2Min*SMax())); }

  /**
   * The maximum allowed value of \f$\sqrt{\hat{s}}\f$.
   */
  Energy mHatMax() const { return min(theMHatMax, sqrt(theX1Max*theX2Max*SMax())); }

  /**
   * The minimum value of the rapidity of the hard sub-process
   * (wrt. the rest system of the colliding particles).
   */
  double yHatMin() const;

  /**
   * The maximum value of the rapidity of the hard sub-process
   * (wrt. the rest system of the colliding particles).
   */
  double yHatMax() const;

  /**
   * Check if the given \f$\hat{y}\f$ is within the cuts.
   */
  bool yHat(double y) const;

  /**
   * The minimum value of the positive light-cone fraction of the hard
   * sub-process.
   */
  double x1Min() const;

  /**
   * The maximum value of the positive light-cone fraction of the hard
   * sub-process.
   */
  double x1Max() const;

  /**
   * Check if the given \f$x_1\f$ is within the cuts.
   */
  bool x1(double x) const;

  /**
   * The minimum value of the negative light-cone fraction of the hard
   * sub-process.
   */
  double x2Min() const;

  /**
   * The maximum value of the negative light-cone fraction of the hard
   * sub-process.
   */
  double x2Max() const;

  /**
   * Check if the given \f$x_2\f$ is within the cuts.
   */
  bool x2(double x) const;

  /**
   * The minimum allowed value of the scale to be used in PDF's and
   * coupling constants.
   */
  Energy2 scaleMin() const { return theScaleMin; }

  /**
   * The maximum allowed value of the scale to be used in PDF's and
   * coupling constants.
   */
  Energy2 scaleMax() const { return theScaleMax; }

  /**
   * Check if the given scale is within the cuts.
   */
  bool scale(Energy2 Q2) const { return Q2 > scaleMin() && Q2 < scaleMax(); }

  /**
   * Set true if a matrix element is should be using this cut and is
   * mirrored along the z-axis .
   */
  bool subMirror() const { return theSubMirror; }

  /**
   * Return the overall cut weight
   */
  double cutWeight() const { return theCutWeight; }

  /**
   * Set the cut weight as appropriate from the call to the last n-cut
   * object.
   */
  void lastCutWeight(double w) const { theLastCutWeight = w; }

  /**
   * Return the fuzziness object
   */
  Ptr<FuzzyTheta>::tcptr fuzzy() const { return theFuzzyTheta; }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isInside(const Value& v, const Value& lower, const Value& upper, double& weight) const {
    if ( !fuzzy() ) {
      if ( v >= lower && v <= upper )
	return true;
      weight = 0.0;
      return false;
    }
    return fuzzy()->isInside<CutType>(v,lower,upper,weight);
  }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isLessThan(const Value& v, const Value& upper, double& weight) const {
    if ( !fuzzy() ) {
      if ( v <= upper )
	return true;
      weight = 0.0;
      return false;
    }
    return fuzzy()->isLessThan<CutType>(v,upper,weight);
  }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isLargerThan(const Value& v, const Value& lower, double& weight) const {
    if ( !fuzzy() ) {
      if ( v >= lower )
	return true;
      weight = 0.0;
      return false;
    }
    return fuzzy()->isLargerThan<CutType>(v,lower,weight);
  }
  //@}

public:

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
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
   * Helper function used by the interface.
   */
  Energy maxMHatMin() const;

  /**
   * Helper function used by the interface.
   */
  Energy minMHatMax() const;

  /**
   * Helper function used by the interface.
   */
  double maxYHatMin() const;

  /**
   * Helper function used by the interface.
   */
  double minYHatMax() const;

  /**
   * Helper function used by the interface.
   */
  double maxX1Min() const;

  /**
   * Helper function used by the interface.
   */
  double minX1Max() const;

  /**
   * Helper function used by the interface.
   */
  double maxX2Min() const;

  /**
   * Helper function used by the interface.
   */
  double minX2Max() const;

  /**
   * Helper function used by the interface.
   */
  Energy2 maxScaleMin() const;

  /**
   * Helper function used by the interface.
   */
  Energy2 minScaleMax() const;

private:

  /**
   * The maximum allowed total invariant mass squared allowed for
   * events to be considered.
   */
  Energy2 theSMax;

  /**
   * The total rapidity of the colliding particles corresponding to
   * the maximum invariant mass squared, SMax().
   */
  double theY;

  /**
   * The invariant mass squared of the hard sub-process of the event
   * being considered.
   */
  mutable Energy2 theCurrentSHat;

  /**
   * The total rapidity of hard sub-process (wrt. the rest system of
   * the colliding particles so that currentYHat() + Y() gives the
   * true rapidity) of the event being considered.
   */
  mutable double theCurrentYHat;

  /**
   * The minimum allowed value of \f$\sqrt{\hat{s}}\f$.
   */
  Energy theMHatMin;

  /**
   * The maximum allowed value of \f$\sqrt{\hat{s}}\f$.
   */
  Energy theMHatMax;

  /**
   * The minimum value of the rapidity of the hard sub-process
   * (wrt. the rest system of the colliding particles).
   */
  double theYHatMin;

  /**
   * The maximum value of the rapidity of the hard sub-process
   * (wrt. the rest system of the colliding particles).
   */
  double theYHatMax;

  /**
   * The minimum value of the positive light-cone fraction of the hard
   * sub-process.
   */
  double theX1Min;

  /**
   * The maximum value of the positive light-cone fraction of the hard
   * sub-process.
   */
  double theX1Max;

  /**
   * The minimum value of the negative light-cone fraction of the hard
   * sub-process.
   */
  double theX2Min;

  /**
   * The maximum value of the negative light-cone fraction of the hard
   * sub-process.
   */
  double theX2Max;

  /**
   * The minimum allowed value of the scale to be used in PDF's and
   * coupling constants.
   */
  Energy2 theScaleMin;

  /**
   * The maximum allowed value of the scale to be used in PDF's and
   * coupling constants.
   */
  Energy2 theScaleMax;

  /**
   * The objects defining cuts on single outgoing partons from the
   * hard sub-process.
   */
  OneCutVector theOneCuts;

  /**
   * The objects defining cuts on pairs of particles in the hard
   * sub-process.
   */
  TwoCutVector theTwoCuts;

  /**
   * The objects defining cuts on sets of outgoing particles from the
   * hard sub-process.
   */
  MultiCutVector theMultiCuts;

  /**
   * An optional jet finder used to define cuts on the level of
   * reconstructed jets.
   */
  Ptr<JetFinder>::ptr theJetFinder;

  /**
   * Set to true if a matrix element is should be using this cut and is
   * mirrored along the z-axis .
   */
  mutable bool theSubMirror;

  /**
   * The overall cut weight
   */
  mutable double theCutWeight;

  /**
   * The cut weight as appropriate from the call to the last n-cut
   * object.
   */
  mutable double theLastCutWeight;

  /**
   * The fuzziness object
   */
  Ptr<FuzzyTheta>::ptr theFuzzyTheta;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Cuts> initCuts;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Cuts & operator=(const Cuts &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Cuts. */
template <>
struct BaseClassTrait<Cuts,1> {
  /** Typedef of the first base class of Cuts. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Cuts class and the shared object where it is defined. */
template <>
struct ClassTraits<Cuts>
  : public ClassTraitsBase<Cuts> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Cuts"; }
};

/** @endcond */

}

#endif /* THEPEG_Cuts_H */
