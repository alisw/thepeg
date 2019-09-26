// -*- C++ -*-
//
// MEBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MEBase_H
#define ThePEG_MEBase_H
// This is the declaration of the MEBase class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/MatrixElement/DiagramBase.h"
#include "ThePEG/MatrixElement/ColourLines.h"
#include "ThePEG/MatrixElement/Amplitude.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/Handlers/StandardXComb.fh"
#include "ReweightBase.h"

#include "ThePEG/Handlers/EventHandler.fh"
#include "ThePEG/Handlers/StandardEventHandler.fh"
#include "ThePEG/Handlers/SubProcessHandler.fh"
#include "ThePEG/PDF/PartonBin.fh"

#include "MEBase.fh"

namespace ThePEG {

/**
 * The MEBase class is the base class of all objects
 * representing hard matrix elements in ThePEG. There are three
 * methods which must be overridden by a concrete subclass:<BR>
 *
 * includedDiagrams(tcPDPair) should return a vector of DiagramBase
 * objects describing the diagrams used for this matrix element for
 * the given pair of incoming parton types. These DiagramBases are
 * used to identify the incoming and outgoing partons which can be
 * handled by the process generation scheme, and is also used to
 * cnstruct a corresponding SubProcess object.
 *
 * scale() should return the scale associated with the phase space
 * point set with the last call to setKinematics(...) or
 * generateKinematics(...).
 *
 * me() should return the the matrix element squared using the the
 * type and momentum of the incoming and outgoing partons, previously
 * set by the setKinematics(...) or generateKinematics(...) member
 * functions, accessible through the methods meMomenta() and
 * mePartonData() inherited from LastXCombInfo, and/or from
 * information stored by sub classes. The returned value should be
 * dimensionless suitable scaled by the total invariant mass squared
 * (accessible through the sHat() member function). Any user of this
 * method must make sure that the setKinematics(...) member function
 * has been appropriately called before.
 *
 * colourGeometries() should return a Selector with the possible
 * ColourLines objects weighted by their relative probabilities given
 * the information set by the last call to setKinematics(...) or
 * generateKinematics(...).
 *
 * There are other virtula functions which may be overridden as listed
 * below.
 *
 * @see \ref MEBaseInterfaces "The interfaces"
 * defined for MEBase.
 * @see DiagramBase
 * @see ColourLines
 * 
 */
class MEBase: public HandlerBase, public LastXCombInfo<StandardXComb> {

public:

  /** A vector of pointers to DiagramBase objects. */
  typedef vector<DiagPtr> DiagramVector;
  /** The size_type used in the DiagramVector. */
  typedef DiagramVector::size_type DiagramIndex;
  /** A vector of pointers to ReweightBase objects. */
  typedef vector<ReweightPtr> ReweightVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MEBase();

  /**
   * Destructor.
   */
  virtual ~MEBase();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes.. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const = 0;

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const = 0;

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const = 0;

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const = 0;

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const;

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries().
   */
  void setKinematics(tPPair in, const PVector & out);

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics() {}

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub);

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub, const ColourLines* cl);

  /**
   * The number of internal degreed of freedom used in the matrix
   * element. This default version returns 0;
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector.
   */
  virtual bool generateKinematics(const double * r) = 0;

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { return true; }

  /**
   * If this is a dependent matrix element in a ME group, return true,
   * if cuts should be inherited from the head matrix element, i.e. no
   * cut is being applied to the dependent matrix element if the head
   * configuration has passed the cuts.
   */
  virtual bool headCuts() const { return false; }

  /**
   * If this is a dependent matrix element in a ME group, return true,
   * if cuts should be ignored.
   */
  virtual bool ignoreCuts() const { return false; }

  /**
   * If this is a dependent matrix element in a ME group, return true,
   * if it applies to the process set in lastXComb()
   */
  virtual bool apply() const { return true; }

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const = 0;

  /**
   * Return true, if this matrix element will generate momenta for the
   * incoming partons itself.  The matrix element is required to store
   * the incoming parton momenta in meMomenta()[0,1]. No mapping in
   * tau and y is performed by the PartonExtractor object, if a
   * derived class returns true here. The phase space jacobian is to
   * include a factor 1/(x1 x2).
   */
  virtual bool haveX1X2() const { return false; }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return false; }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { return false; }

  /**
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return false; }

  /**
   * Comlete a SubProcess object using the internal degrees of freedom
   * generated in the last generateKinematics() (and possible other
   * degrees of freedom which was intergated over in dSigHatDR(). This
   * default version does nothing. Will be made purely virtual in the
   * future.
   */
  virtual void generateSubCollision(SubProcess &);

  /**
   * Clear the information previously provided by a call to
   * setKinematics(...).
   */
  virtual void clearKinematics();

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const = 0;

  /**
   * Return true, if this matrix element does not want to
   * make use of mirroring processes; in this case all
   * possible partonic subprocesses with a fixed assignment
   * of incoming particles need to be provided through the diagrams
   * added with the add(...) method.
   */
  virtual bool noMirror () const { return false; }

  /**
   * Return all possible diagrams.
   */
  const DiagramVector & diagrams() const {
    if ( theDiagrams.empty() ) getDiagrams();
    return theDiagrams;
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const = 0;

  /**
   * Select a ColpurLines geometry. The default version returns a
   * colour geometry selected among the ones returned from
   * colourGeometries(tcDiagPtr).
   */
  virtual const ColourLines &
  selectColourGeometry(tcDiagPtr diag) const;

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const {
    return Selector<DiagramIndex>();
  }


  /**
   * Select a diagram. Default version uses diagrams(const
   * DiagramVector &) to select a diagram according to the
   * weights. This is the only method used that should be outside of
   * MEBase.
   */
  virtual DiagramIndex diagram(const DiagramVector &) const;

  /**
   * Return true if this matrix element has associated (p)reWeight
   * objects assigned.
   */
  inline bool reweighted() const {
    return reweights.size() > 0 || preweights.size() > 0;
  }

  /**
   * With the information previously supplied with the
   * setKinematics(...) methods, return the combined effects of the
   * reweighters.
   */
  double reWeight() const;

  /**
   * With the information previously supplied with the
   * setKinematics(...) methods, return the comined effects of the
   * peweighters.
   */
  double preWeight() const;

  /**
   * Add objects to the list of reweighters.
   */
  void addReweighter(tReweightPtr rw);

  /**
   * Add objects to the list of preweighters.
   */
  void addPreweighter(tReweightPtr rw);

  /**
   * Return the amplitude associated with this matrix element. This
   * function is allowed to return the null pointer if the amplitude
   * is not available.
   */
  Ptr<Amplitude>::pointer amplitude() const { return theAmplitude; }

  /**
   * Set the amplitude associated with this matrix element.
   */
  void amplitude(Ptr<Amplitude>::pointer amp) { theAmplitude = amp; }
  //@}

public:

  /** @name Acces information about the last generated phase space point. */
  //@{
  /**
   * Return the last set invariant mass squared.
   */
  Energy2 sHat() const { return lastSHat(); }

  /**
   * Return the factor with which this matrix element was last
   * pre-weighted.
   */
  double preweight() const { return lastPreweight(); }

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches() {}

  /**
   * For the given event generation setup return a xcomb object
   * appropriate to this matrix element.
   */
  virtual StdXCombPtr makeXComb(Energy newMaxEnergy, const cPDPair & inc,
				tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				const PBPair & newPartonBins, tCutsPtr newCuts,
				const DiagramVector & newDiagrams, bool mir,
				const PartonPairVec& allPBins,
				tStdXCombPtr newHead = tStdXCombPtr(),
				tMEPtr newME = tMEPtr());

  /**
   * For the given event generation setup return a dependent xcomb object
   * appropriate to this matrix element.
   */
  virtual StdXCombPtr makeXComb(tStdXCombPtr newHead,
				const PBPair & newPartonBins,
				const DiagramVector & newDiagrams,
				tMEPtr newME = tMEPtr());

  /**
   * Fill the projectors object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  virtual void fillProjectors() { }

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr);

  /**
   * Retrieve information obtained in the calculation of the cross
   * section to be used later when selecting diagrams and colour flow.
   */
  const DVector & meInfo() const;

  /**
   * Save information obtained in the calculation of the cross
   * section to be used later when selecting diagrams and colour flow.
   */
  void meInfo(const DVector & info) const;

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  virtual int maxMultCKKW() const { return theMaxMultCKKW; }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  virtual int minMultCKKW() const { return theMinMultCKKW; }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this will set the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  virtual void maxMultCKKW(int mult) { theMaxMultCKKW = mult; }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this will set the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  virtual void minMultCKKW(int mult) { theMinMultCKKW = mult; }

  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const {}
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * To be used by sub classes in the getDiagrams() method to add
   * included diagrams.
   */
  void add(DiagPtr dp) const { theDiagrams.push_back(dp); }

  /**
   * Access the momenta set by the last call to generateKinematics().
   */
  vector<Lorentz5Momentum> & meMomenta();
  using LastXCombInfo<StandardXComb>::meMomenta;

  /**
   * Set the matrix element squared as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  void lastME2(double v) const;
  using LastXCombInfo<StandardXComb>::lastME2;

  /**
   * Set the last preweight factor
   */
  void lastPreweight(double w) const;
  using LastXCombInfo<StandardXComb>::lastPreweight;

  /**
   * Set the partonic cross section as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  void lastMECrossSection(CrossSection v) const;
  using LastXCombInfo<StandardXComb>::lastMECrossSection;

  /**
   * Set the PDF weight as calculated
   * for the last phase space point, if the matrix
   * element does supply PDF weights. This may optionally
   * be used by a matrix element for caching.
   */
  void lastMEPDFWeight(double v) const;
  using LastXCombInfo<StandardXComb>::lastMEPDFWeight;

  /**
   * Set the coupling weight as calculated
   * for the last phase space point
   */
  void lastMECouplings(double v) const;
  using LastXCombInfo<StandardXComb>::lastMECouplings;

  /**
   * Set the last jacobian obtained when generating the kinematics for
   * the call to dSigHatDR.
   */
  void jacobian(double j);
  using LastXCombInfo<StandardXComb>::jacobian;

  /**
   * Initialize all member variables from another
   * MEBase object.
   *
   * @TODO remove?
   */
  void use(tcMEPtr other);

  /**
   * Initialize the diagrams from another MEBase object.
   */
  void useDiagrams(tcMEPtr other) const;

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The diagrams included for this matrix element.
   */
  mutable DiagramVector theDiagrams;

  /**
   * The reweight objects modifying this matrix element.
   */
  ReweightVector reweights;

  /**
   * The preweight objects modifying this matrix element.
   */
  ReweightVector preweights;

  /**
   * The amplitude associated with this matrix element.
   */
  Ptr<Amplitude>::pointer theAmplitude;

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  int theMaxMultCKKW;

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  int theMinMultCKKW;

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<MEBase> initMEBase;

  /**
   *  Private and non-existent assignment operator.
   */
  MEBase & operator=(const MEBase &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * MEBase.
 */
template <>
struct BaseClassTrait<MEBase,1>: public ClassTraitsType {
  /** Typedef of the base class of MEBase. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * MEBase class.
 */
template <>
struct ClassTraits<MEBase>: public ClassTraitsBase<MEBase> {
  /** Return the class name. */
  static string className() { return "ThePEG::MEBase"; }
};

/** @endcond */

}

#include "ThePEG/Handlers/StandardXComb.h"

#endif /* ThePEG_MEBase_H */
