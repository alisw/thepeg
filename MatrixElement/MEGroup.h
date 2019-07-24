// -*- C++ -*-
//
// MEGroup.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MEGroup_H
#define ThePEG_MEGroup_H
// This is the declaration of the MEGroup class.

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Cuts/Cuts.fh"
#include "MEGroup.fh"

namespace ThePEG {

/**
 * The MEGroup class represents a 'head' matrix element
 * in association with a group of dependent matrix elements.
 * It basically acts as a wrapper around its head matrix element
 * however supplying additional information to the corresponding
 * StdXCombGroup object.
 *
 * @see StdXCombGroup
 * 
 */
class MEGroup: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MEGroup();

  /**
   * Destructor.
   */
  virtual ~MEGroup();
  //@}

public:

  /** @name Virtual functions from MEBase. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return head()->orderInAlphaS(); }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return head()->orderInAlphaEW(); }

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const { return head()->me2(); }

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const { return head()->scale(); }

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const { return head()->alphaS(); }

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const { return head()->alphaEM(); }

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics() {
    MEBase::setKinematics();
    head()->setKinematics();
  }

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub) { head()->constructVertex(sub); }

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub, const ColourLines* cl) {
    head()->constructVertex(sub,cl);
  }

  /**
   * The number of internal degreed of freedom used in the matrix
   * element. This default version returns 0;
   */
  virtual int nDim() const { return theNDim; }

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS () const { return head()->wantCMS(); }

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const { return head()->dSigHatDR(); }

  /**
   * Return true, if this matrix element will generate momenta for the
   * incoming partons itself.  The matrix element is required to store
   * the incoming parton momenta in meMomenta()[0,1]. No mapping in
   * tau and y is performed by the PartonExtractor object, if a
   * derived class returns true here. The phase space jacobian is to
   * include a factor 1/(x1 x2).
   */
  virtual bool haveX1X2() const { return head()->haveX1X2(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1 () const { return head()->havePDFWeight1(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2 () const { return head()->havePDFWeight2(); }

  /**
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return head()->keepRandomNumbers(); }

  /**
   * Comlete a SubProcess object using the internal degrees of freedom
   * generated in the last generateKinematics() (and possible other
   * degrees of freedom which was intergated over in dSigHatDR(). This
   * default version does nothing. Will be made purely virtual in the
   * future.
   */
  virtual void generateSubCollision(SubProcess & sub) { head()->generateSubCollision(sub); }

  /**
   * Clear the information previously provided by a call to
   * setKinematics(...).
   */
  virtual void clearKinematics();

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const {
    head()->diagrams();
    useDiagrams(head());
  }

  /**
   * Return true, if this matrix element does not want to
   * make use of mirroring processes; in this case all
   * possible partonic subprocesses with a fixed assignment
   * of incoming particles need to be provided through the diagrams
   * added with the add(...) method.
   */
  virtual bool noMirror () const { return head()->noMirror(); }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const { return head()->colourGeometries(diag); }

  /**
   * Select a ColpurLines geometry. The default version returns a
   * colour geometry selected among the ones returned from
   * colourGeometries(tcDiagPtr).
   */
  virtual const ColourLines &
  selectColourGeometry(tcDiagPtr diag) const { return head()->selectColourGeometry(diag); }

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const  { 
    return head()->diagrams(dv); 
  }

  /**
   * Select a diagram. Default version uses diagrams(const
   * DiagramVector &) to select a diagram according to the
   * weights. This is the only method used that should be outside of
   * MEBase.
   */
  virtual DiagramIndex diagram(const DiagramVector & dv) const {
    DiagramIndex res = head()->diagram(dv); 
    return res;
  }

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr xc) {
    MEBase::setXComb(xc);
    head()->setXComb(xc);
  }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  virtual int maxMultCKKW() const { return head()->maxMultCKKW(); }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  virtual int minMultCKKW() const { return head()->minMultCKKW(); }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this will set the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  virtual void maxMultCKKW(int mult) { head()->maxMultCKKW(mult); }

  /**
   * If this matrix element is to be used together with others for
   * CKKW reweighting and veto, this will set the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  virtual void minMultCKKW(int mult) { head()->minMultCKKW(mult); }

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches() { head()->flushCaches(); }

  /**
   * Collect information on the last evaluated phasespace
   * point for verification or debugging purposes. This
   * only called, if the StdXCombGroup did accumulate
   * a non-zero cross section from this ME group.
   */
  virtual void lastEventStatistics() {}
  //@}

public:

  /**
   * Return the head matrix element.
   */
  tMEPtr head() const { return theHead; }

  /**
   * Visit the dependent matrix elements
   */
  const MEVector& dependent() const { return theDependent; }

  /**
   * Set the head matrix element.
   */
  void head(tMEPtr me) { theHead = me; }

  /**
   * Access the dependent matrix elements
   */
  MEVector& dependent() { return theDependent; }

  /**
   * Return the random number offset to access the random
   * numbers provided for the given matrix element to generate
   * dependent kinematics.
   */
  int dependentOffset(tMEPtr dep) const;

  /**
   * For the given event generation setup return an xcomb object
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
   * Create a dependent xcomb object to be used
   * for the given process steered bythe head object and 
   * dependent matrix element.
   */
  virtual vector<StdXCombPtr> makeDependentXCombs(tStdXCombPtr xcHead,
						  const cPDVector& proc,
						  tMEPtr depME,
						  const PartonPairVec& allPBins) const;

  /**
   * Fill the projectors object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  virtual void fillProjectors() { head()->fillProjectors(); }

  /**
   * Return true, if projectors will be used
   */
  virtual bool willProject() const { return false; }

  /**
   * Return true, if this MEGroup will reweight the contributing cross
   * sections.
   */
  virtual bool groupReweighted() const { return false; }

  /**
   * Reweight the head cross section
   */
  virtual double reweightHead(const vector<tStdXCombPtr>&) { return 1.; }

  /**
   * Reweight the dependent cross section
   */
  virtual double reweightDependent(tStdXCombPtr, const vector<tStdXCombPtr>&) { return 1.; }

  /**
   * Return true, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  virtual bool subProcessGroups() const { return true; }

public:

  /**
   * Return true, if the same additional random numbers
   * should be presented to any of the dependent
   * matrix elements.
   */
  virtual bool uniformAdditional() const = 0;

  /**
   * Given a process from the head matrix element,
   * return a list of diagrams which should be considered for
   * the given dependent matrix element.
   */
  virtual MEBase::DiagramVector dependentDiagrams(const cPDVector& proc,
						  tMEPtr depME) const = 0;

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

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   * The head matrix element.
   */
  MEPtr theHead;

  /**
   * The dependent matrix elements.
   */
  MEVector theDependent;

  /**
   * Offsets to access additional random numbers
   * required by the dependent matrix elements.
   */
  map<tMEPtr,int> theNDimMap;

  /**
   * The total number of random numbers required.
   */
  int theNDim;

private:

  /**
   * Describe a class with persistent data.
   */
  static AbstractClassDescription<MEGroup> initMEGroup;

  /**
   *  Private and non-existent assignment operator.
   */
  MEGroup & operator=(const MEGroup &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * MEGroup.
 */
template <>
struct BaseClassTrait<MEGroup,1> {
  /** Typedef of the base class of MEGroup. */
  typedef MEBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * MEGroup class.
 */
template <>
struct ClassTraits<MEGroup>: public ClassTraitsBase<MEGroup> {
  /** Return the class name. */
  static string className() { return "ThePEG::MEGroup"; }
};

/** @endcond */

}

#endif /* ThePEG_MEGroup_H */
