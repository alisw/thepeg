// -*- C++ -*-
//
// StandardXComb.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
// Copyright (C) 2009-2011 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardXComb_H
#define ThePEG_StandardXComb_H
// This is the declaration of the StandardXComb class.

#include "ThePEG/Config/ThePEG.h"
#include "SubProcessHandler.fh"
#include "ThePEG/PDF/PartonExtractor.fh"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/Utilities/VSelector.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Utilities/XSecStat.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/SubProcessHandler.fh"
#include "StandardXComb.fh"

namespace ThePEG {

/**
 * The StandardXComb class inherits from the more general XComb class
 * which stores all information about the generation of a hard
 * sub-proces for a given pair of incoming particles, a pair of
 * extracted partons, etc. This class stores more information related
 * to thestandard process generation scheme in ThePEG, such as the
 * PartonExtractor and MEBase object used. It also does some of the
 * administration of the process generation.
 *
 * The main function is dSigDR() which returns the differential cross
 * section w.r.t. a given vector of random numbers in the interval
 * ]0,1[. In the initialization this is used to pre-sample the phase
 * space. In the generation phase it is used to give the cross section
 * for a phase space point, and if this StandardXComb is chosen the
 * construct() function is called to generate the actual sub-process.
 *
 * @see ParonExtractor
 * @see MEBase
 * @see Cuts
 * @see StdXCombGroup
 */
class StandardXComb: public XComb {

public:

  /** A vector of DiagramBase objects. */
  typedef MEBase::DiagramVector DiagramVector;

  /** A vector of indices. */
  typedef MEBase::DiagramIndex DiagramIndex;

  /** MEBase needs to be a friend. */
  friend class MEBase;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Standard constructor.
   */
  StandardXComb(Energy newMaxEnergy, const cPDPair & inc,
		tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
		tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
		const PBPair & newPartonBins, tCutsPtr newCuts, tMEPtr newME,
		const DiagramVector & newDiagrams, bool mir,
		tStdXCombPtr newHead = tStdXCombPtr());

  /**
   * Constructor given a head xcomb.
   */
  StandardXComb(tStdXCombPtr newHead,
		const PBPair & newPartonBins, tMEPtr newME,
		const DiagramVector & newDiagrams);

  /**
   * Default constructor.
   */
  StandardXComb();

  /**
   * Destructor.
   */
  virtual ~StandardXComb();

  /**
   * Constructor used by MEBase to create a temporary object to store info.
   */
  StandardXComb(tMEPtr me, const tPVector & parts, DiagramIndex i);

  //@}

  /** @name Utilities for incoming partons. */
  //@{

  /**
   * Properly setup the PartonBinInstance objects provided a sub
   * process has been constructed using this XComb.
   */
  void recreatePartonBinInstances(Energy2 scale);

  /**
   * Fill the variables needed to generate remnants; momenta will be
   * used from the partons set in this xcomb, but random numbers need
   * to be provided to (re)generate variables not fixed by the
   * incoming partons.
   */
  void refillPartonBinInstances(const double* r);

  /**
   * Setup information on incoming partons depending
   * on the information previously supplied through the
   * choice of diagram and incoming momenta in the first
   * two entries of meMomenta(). Partons are not actually
   * extracted from the incoming particles, though a subprocess
   * detached from the current Event may be created.
   */
  bool setIncomingPartons(tStdXCombPtr labHead = tStdXCombPtr());

  /**
   * Fill phase space information as far as possible
   */
  void fill(const PPair& newParticles,
	    const PPair& newPartons,
	    const vector<Lorentz5Momentum>& newMEMomenta,
	    const DVector& newLastRandomNumbers = DVector());

  //@}

  /** @name Access the assigned objects used in the generation. */
  //@{
  /**
   * Return a pointer to the corresponding sub-process handler. May be
   * null if the standard process generation in ThePEG was not used.
   */
  tcSubHdlPtr subProcessHandler() const { return theSubProcessHandler; }

  /**
   * The matrix element to be used.
   */
  tMEPtr matrixElement() const { return theME; }

  /**
   * Return a pointer to the head XComb this XComb
   * depends on. May return NULL, if this is not a
   * member of a XComb group.
   */
  tStdXCombPtr head() const { return theHead; }

  /**
   * Set the head XComb pointer.
   */
  void head(tStdXCombPtr headXC) { theHead = headXC; }

  /**
   * Return a selector object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  Selector<tStdXCombPtr>& projectors() { return theProjectors; }

  /**
   * Return a selector object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  const Selector<tStdXCombPtr>& projectors() const { return theProjectors; }

  /**
   * Return a pointer to a projector xcomb which will generate a subprocess
   * different from the one just integrated.
   */
  tStdXCombPtr lastProjector() const { return theProjector; }

  /**
   * Set a pointer to a projector xcomb which will generate a subprocess
   * different from the one just integrated.
   */
  void lastProjector(tStdXCombPtr pxc) { theProjector = pxc; }

  //@}

  /** @name Main functions used for the generation. */
  //@{
  /**
   * Try to determine if this subprocess is at all possible.
   */
  virtual bool checkInit();

  /**
   * The number of dimensions of the phase space used to generate this
   * process.
   */
  virtual int nDim() const { return theNDim; }

  /**
   * Return the parton extraction dimensions
   */
  const pair<int,int>& partonDimensions() const { return partonDims; }

  /**
   * Return true, if the current configuration will pass the cuts
   */
  bool willPassCuts();

  /**
   * Return the cut weight encountered from fuzzy cuts
   */
  double cutWeight() const { return theCutWeight; }

  /**
   * Reset all saved data about last generated phasespace point;
   */
  virtual void clean();

  /**
   * Return true, if kinematics have already been generated
   */
  bool kinematicsGenerated() const  { return theKinematicsGenerated; }

  /**
   * Indicate that kinematics have been generated
   */
  void didGenerateKinematics() { theKinematicsGenerated = true; }

  /**
   * Generate a phase space point from a vector \a r of \a nr numbers
   * in the interval ]0,1[ and return the corresponding differential
   * cross section.
   */
  virtual CrossSection dSigDR(const pair<double,double> ll, int nr, const double * r);

  /**
   * If this XComb has a head XComb, return the cross section
   * differential in the variables previously supplied. The PDF weight
   * is taken from the lastPDFWeight supplied by the head XComb
   * object.
   */
  CrossSection dSigDR(const double * r);

  /**
   * Return the PDF weight used in the last call to dSigDR
   */
  double lastPDFWeight() const { return theLastPDFWeight; }

  /**
   * Return the cross section calculated in the last call to dSigDR
   */
  CrossSection lastCrossSection() const { return theLastCrossSection; }

  /**
   * Check if a reshuffling is required when constructing the hard
   * subprocess.
   */
  void checkReshufflingNeeds();

  /**
   * Return true if a reshuffling is required when constructing the hard
   * subprocess. 
   */
  bool needsReshuffling() const { return theNeedsReshuffling; }

  /**
   * Perform the reshuffling from hardProcessMass to mass values,
   * given outgoing momenta
   */
  void reshuffle(vector<Lorentz5Momentum>&) const;

  /**
   * Construct a sub-process object from the information available.
   */
  virtual tSubProPtr construct();
  //@}

  /** @name Functions used for collecting statistics. */
  //@{
  /**
   * The statistics object for this XComb.
   */
  virtual const XSecStat & stats() const { return theStats; }

  /**
   * Select the current event. It will later be rejected with a
   * probability given by \a weight.
   */
  virtual void select(double weight) { theStats.select(weight); }

  /**
   * Accept the current event assuming it was previously selcted.
   */
  virtual void accept() { theStats.accept(); }

  /**
   * Reweight a selected and accepted event.
   */
  void reweight(double oldWeight, double newWeight) {
    theStats.reweight(oldWeight,newWeight);
  }

  /**
   * Reject the current event assuming it was previously accepted. If
   * weighted events are produced, the \a weight should be the same as
   * the previous call to select(double).
   */
  virtual void reject(double weight = 1.0) { theStats.reject(weight); }

  /**
   * Reset statistics.
   */
  virtual void reset() { theStats.reset(); }
  //@}

  /** @name Access information used by the MEBase object. */
  //@{
  /**
   * The diagrams used by the matrix element.
   */
  const DiagramVector & diagrams() const { return theDiagrams; }

  /**
   * True if the TreeDiagram's for this matrix element should in fact
   * be mirrored before used to create an actual sub-rocess.
   */
  bool mirror() const { return isMirror; }

  /**
   * Return the momenta of the partons to be used by the matrix
   * element object, in the order specified by the TreeDiagram objects
   * given by the matrix element.
   */
  const vector<Lorentz5Momentum> & meMomenta() const { return theMEMomenta; }

  /**
   * Return the last selected diagram.
   */
  tcDiagPtr lastDiagram() const { 
    if ( !theExternalDiagram )
      return diagrams()[lastDiagramIndex()];
    return theExternalDiagram;
  }

  /**
   * Return the parton types to be used by the matrix element object,
   * in the order specified by the TreeDiagram objects given by the
   * matrix element.
   */
  const cPDVector & mePartonData() const { return theMEPartonData; }

  /**
   * Return the index of the last selected diagram.
   */
  DiagramIndex lastDiagramIndex() const { return theLastDiagramIndex; }

  /**
   * Get information saved by the matrix element in the calculation of
   * the cross section to be used later when selecting diagrams and
   * colour flow.
   */
  const DVector & meInfo() const { return theMEInfo; }

  /**
   * Set information saved by the matrix element in the calculation of
   * the cross section to be used later when selecting diagrams and
   * colour flow.
   */
  void meInfo(const DVector & info) { theMEInfo = info; }

  /**
   * Return the random numbers used to generate the
   * last phase space point, if the matrix element
   * requested so.
   */
  const DVector& lastRandomNumbers() const { return theLastRandomNumbers; }

  /**
   * Get the last jacobian obtained when generating the kinematics
   * for the call to dSigHatDR.
   */
  double jacobian() const { return theLastJacobian; }

  /**
   * Return the matrix element squared as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  double lastME2() const { return theLastME2; }

  /**
   * Return the last preweight factor
   */
  double lastPreweight() const { return theLastPreweight; }

  /**
   * Return the partonic cross section as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  CrossSection lastMECrossSection() const { return theLastMECrossSection; }

  /**
   * Return the PDF weight as calculated
   * for the last phase space point, if the matrix
   * element does supply PDF weights. This may optionally
   * be used by a matrix element for caching.
   */
  double lastMEPDFWeight() const { return theLastMEPDFWeight; }

  /**
   * Return the coupling factor as calculated for the lats phase space
   * point.
   */
  double lastMECouplings() const { return theLastMECouplings; }
  //@}

  /**
   * Construct the corresponding SubProcess object if it hasn't been
   * done before.
   */
  virtual void newSubProcess(bool group = false);

  /**
   * Return the momenta of the partons to be used by the matrix
   * element object, in the order specified by the TreeDiagram objects
   * given by the matrix element.
   */
  vector<Lorentz5Momentum> & meMomenta() { return theMEMomenta; }

  /**
   * Access the random numbers used to generate the
   * last phase space point, if the matrix element
   * requested so.
   */
  DVector& lastRandomNumbers() { return theLastRandomNumbers; }

  /**
   * Return the parton types to be used by the matrix element object,
   * in the order specified by the TreeDiagram objects given by the
   * matrix element.
   */
  cPDVector & mePartonData() { return theMEPartonData; }

  /**
   * Set a diagram to be used instead of the one selected by the matrix
   * element.
   */
  void externalDiagram(tcDiagPtr diag) { theExternalDiagram = diag; }

  /**
   * Set the last selected diagram.
   */
  void lastDiagramIndex(DiagramIndex i) { theLastDiagramIndex = i; }

  /**
   * Set the PDF weight used in the last call to dSigDR
   */
  void lastPDFWeight(double w) { theLastPDFWeight = w; }

  /**
   * Set the cross section calculated in the last call to dSigDR
   */
  void lastCrossSection(CrossSection s) { theLastCrossSection = s; }

  /**
   * Set the last jacobian obtained when generating the kinematics for
   * the call to dSigHatDR.
   */
  void jacobian(double j) { theLastJacobian = j; }

  /**
   * Set the matrix element squared as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  void lastME2(double v) { theLastME2 = v; }

  /**
   * Set the last preweight factor
   */
  void lastPreweight(double w) { theLastPreweight = w; }

  /**
   * Set the partonic cross section as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  void lastMECrossSection(CrossSection v) { theLastMECrossSection = v; }

  /**
   * Set the PDF weight as calculated
   * for the last phase space point, if the matrix
   * element does supply PDF weights. This may optionally
   * be used by a matrix element for caching.
   */
  void lastMEPDFWeight(double v) { theLastMEPDFWeight = v; }

  /**
   * Set the coupling factor
   */
  void lastMECouplings(double v) { theLastMECouplings = v; }

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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * The corresponding sub-process handler
   */
  tSubHdlPtr theSubProcessHandler;

  /**
   * The matrix element to be used.
   */
  tMEPtr theME;

  /**
   * Statistics gathering for this XComb.
   */
  XSecStat theStats;

  /**
   * The diagrams used by the matrix element.
   */
  DiagramVector theDiagrams;

  /**
   * True if the TreeDiagram's for this matrix element should in fact
   * be mirrored before used to create an actual sub-rocess.
   */
  bool isMirror;

  /**
   * The number of dimensions of the phase space used to generate this
   * process.
   */
  int theNDim;

protected:

  /**
   * The number of dimensions of the phase space used for each of the
   * incoming partons.
   */
  pair<int,int> partonDims;

private:

  /**
   * True, if kinematics have already been generated
   */
  bool theKinematicsGenerated;

  /**
   * The momenta of the partons to be used by the matrix element
   * object, in the order specified by the TreeDiagram objects given
   * by the matrix element.
   */
  vector<Lorentz5Momentum> theMEMomenta;

  /**
   * The parton types to be used by the matrix element object, in the
   * order specified by the TreeDiagram objects given by the matrix
   * element.
   */
  cPDVector theMEPartonData;

  /**
   * A diagram to be used instead of the one selected by the matrix element.
   */
  tcDiagPtr theExternalDiagram;

  /**
   * The last selected tree diagram.
   */
  DiagramIndex theLastDiagramIndex;

  /**
   * Information saved by the matrix element in the calculation of the
   * cross section to be used later when selecting diagrams and colour
   * flow.
   */
  DVector theMEInfo;

  /**
   * The random numbers used to generate the
   * last phase space point, if the matrix element
   * requested so.
   */
  DVector theLastRandomNumbers;

  /**
   * The PDF weight used in the last call to dSigDR
   */
  double theLastPDFWeight;

  /**
   * The cross section calculated in the last call to dSigDR
   */
  CrossSection theLastCrossSection;

  /**
   * Save the last jacobian obtained when generating the kinematics for
   * the call to dSigHatDR.
   */
  double theLastJacobian;

  /**
   * The matrix element squared as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  double theLastME2;

  /**
   * The last preweight factor
   */
  double theLastPreweight;

  /**
   * The partonic cross section as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  CrossSection theLastMECrossSection;

  /**
   * The PDF weight as calculated
   * for the last phase space point, if the matrix
   * element does supply PDF weights. This may optionally
   * be used by a matrix element for caching.
   */
  double theLastMEPDFWeight;

  /**
   * The coupling factor
   */
  double theLastMECouplings;

  /**
   * A pointer to the head XComb this XComb
   * depends on. May return NULL, if this is not a
   * member of a XComb group.
   */
  tStdXCombPtr theHead;

  /**
   * A selector object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  Selector<tStdXCombPtr> theProjectors;

  /**
   * A pointer to a projector xcomb which will generate a subprocess
   * different from the one just integrated.
   */
  tStdXCombPtr theProjector;

  /**
   * True, if cuts have already been checked
   */
  bool checkedCuts;

  /**
   * The result of the last call to willPassCuts
   */
  bool passedCuts;

  /**
   * The cut weight encountered from fuzzy cuts
   */
  double theCutWeight;

  /**
   * True if a reshuffling is required when constructing the hard
   * subprocess.
   */
  bool theNeedsReshuffling;

  /**
   * Calculate the reshuffling equation given the coefficients
   */
  double reshuffleEquation(double, const vector<pair<Energy2,Energy2> >&, Energy2) const;

  /**
   * Solve the reshuffling equation given the coefficients
   */
  double solveReshuffleEquation(const vector<pair<Energy2,Energy2> >&, Energy2) const;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardXComb> initStandardXComb;
 
  /**
   * Private and non-existent assignment operator.
   */
  StandardXComb & operator=(const StandardXComb &);

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * StandardXComb.
 */
template <>
struct BaseClassTrait<StandardXComb,1>: public ClassTraitsType {
  /** Typedef of the base class of StandardXComb. */
  typedef XComb NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * StandardXComb class.
 */
template <>
struct ClassTraits<StandardXComb>:
    public ClassTraitsBase<StandardXComb> {
  /** Return the class name. */
  static string className() { return "ThePEG::StandardXComb"; }
};

/** @endcond */

}

#endif /* ThePEG_StandardXComb_H */
