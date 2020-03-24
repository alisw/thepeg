// -*- C++ -*-
//
// StdXCombGroup.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StdXCombGroup_H
#define ThePEG_StdXCombGroup_H
// This is the declaration of the StdXCombGroup class.

#include "StandardXComb.h"
#include "StdXCombGroup.fh"
#include "ThePEG/MatrixElement/MEGroup.fh"

namespace ThePEG {

/**
 * The StdXCombGroup class represents a 'head' XComb object
 * in association with a group of dependent XComb objects.
 *
 * @see MEGroup
 */
class StdXCombGroup: public StandardXComb {

  /** MEBase needs to be a friend. */
  friend class MEBase;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Standard constructor.
   */
  StdXCombGroup(Energy newMaxEnergy, const cPDPair & inc,
		tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
		tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
		const PBPair & newPartonBins, tCutsPtr newCuts, tMEGroupPtr newME,
		const DiagramVector & newDiagrams, bool mir,
		tStdXCombPtr newHead = tStdXCombPtr());

  /**
   * Default constructor.
   */
  StdXCombGroup();

  /**
   * Destructor.
   */
  virtual ~StdXCombGroup();

public:

  /**
   * Reset all saved data about last generated phasespace point;
   */
  virtual void clean();

  /**
   * The number of dimensions of the phase space used to generate this
   * process.
   */
  virtual int nDim() const;

  /**
   * Generate a phase space point from a vector \a r of \a nr numbers
   * in the interval ]0,1[ and return the corresponding differential
   * cross section.
   */
  virtual CrossSection dSigDR(const pair<double,double> ll, int nr, const double * r);

  /**
   * Return the cross section calculated from the head matrix element
   */
  CrossSection lastHeadCrossSection() const { return theLastHeadCrossSection; }

  /**
   * Visit the dependent XComb objects
   */
  const vector<StdXCombPtr>& dependent() const { return theDependent; }

  /**
   * Return the matrix element group steered by this
   * XComb group.
   */
  tcMEGroupPtr meGroup() const { return theMEGroup; }

  /**
   * Initialize this XComb group
   */
  void build(const PartonPairVec& allPBins);

  /**
   * Construct a sub-process object from the information available.
   */
  virtual tSubProPtr construct();

  /**
   * Construct the corresponding SubProcess object if it hasn't been
   * done before.
   */
  virtual void newSubProcess(bool);

  /**
   * Set the cross section calculated from the head matrix element
   */
  void lastHeadCrossSection(CrossSection xs) { theLastHeadCrossSection = xs; }

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
   * The MEGroup object
   */
  MEGroupPtr theMEGroup;

  /**
   * The dependent XComb objects
   */
  vector<StdXCombPtr> theDependent;

  /**
   * The cross section calculated from the head matrix element
   */
  CrossSection theLastHeadCrossSection;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StdXCombGroup> initStdXCombGroup;
 
  /**
   * Private and non-existent assignment operator.
   */
  StdXCombGroup & operator=(const StdXCombGroup &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * StdXCombGroup.
 */
template <>
struct BaseClassTrait<StdXCombGroup,1> {
  /** Typedef of the base class of StdXCombGroup. */
  typedef StandardXComb NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * StdXCombGroup class.
 */
template <>
struct ClassTraits<StdXCombGroup>:
    public ClassTraitsBase<StdXCombGroup> {
  /** Return the class name. */
  static string className() { return "ThePEG::StdXCombGroup"; }
};

/** @endcond */

}

#endif /* ThePEG_StdXCombGroup_H */
