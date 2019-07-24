// -*- C++ -*-
//
// PartonExtractor.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PartonExtractor_H
#define ThePEG_PartonExtractor_H
// This is the declaration of the PartonExtractor class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "PartonExtractor.xh"

namespace ThePEG {

/**
 * The PartonExtractor is a base class defining the interface to
 * objects responsible for extracting partons from particles. It is
 * used by a SubProcessHandler which combines one PartonExtractor with
 * a number of MEBase objects which are the used in an XComb in a
 * StandardEventHandler to generate hard sub-processes.
 *
 * PartonExtractor inherits from the general HandlerBase class and
 * from the LastXCombInfo class to have easy acces to information
 * about the currently chosen hard sub-process.
 *
 * @see \ref PartonExtractorInterfaces "The interfaces"
 * defined for PartonExtractor.
 * @see SubProcessHandler
 * @see MEBase
 * @see EventHandler
 * @see StandardEventHandler
 * @see XComb
 * @see HandlerBase
 *
 */
class PartonExtractor: public HandlerBase, public LastXCombInfo<> {

  /** XComb is a friend. */
  friend class XComb;

public:

  /** A map of PartonBinInstance objects indexed by the extracted parton. */
  typedef map<cPPtr,PBIPtr> PartonBinInstanceMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PartonExtractor();

  /**
   * Destructor.
   */
  virtual ~PartonExtractor();
  //@}

public:

  /** @name Virtual functions which may be overridden in sub-classes. */
  //@{
  /**
   * Return true if this parton extractor can handle the given types
   * of incoming particles.
   */
  virtual bool canHandle(const cPDPair &) { return true; }

  /**
   * Return a vector of possible pairs of parton bins which can be
   * produced within a given maximum total particle-particle
   * invariant mass squared, \a maxEnergy sBin.
   */
  virtual PartonPairVec getPartons(Energy maxEnergy, const cPDPair &,
				   const Cuts &) const;

  /**
   * May be overriden by sub-classes which have their own oppinion
   * about what scale to use in a hard subprocess. The default version
   * simply returns the previously selected scale.
   */
  virtual Energy2 newScale();

  /**
   * Connect the remnants with the colour lines of the extracted
   * parton.
   */
  virtual void colourConnect(tPPtr particle, tPPtr parton,
			     const tPVector & remnants) const;

  /**
   * If remnants has already been created for the given parton, remove
   * them from the given step and generate new remnants corresponding
   * to the parton newp and add them to the step. The new parton bins
   * are returned.
   * @throws Veto if remnant generation failed for whatever reason.
   */
  virtual PBIPair newRemnants(tPPair oldp, tPPair newp, tStepPtr step);

  /**
   * Determine the number of random numbers needed to calculate
   * \f$\hat{s}\f$ and the product of all densitiy functions.
   */
  virtual pair<int,int> nDims(const PBPair & pbins);

  /**
   * Prepare the given parton bin instances for generating a new
   * event.
   */
  virtual void prepare(const PBIPair & pbins);

  /**
   * Update information on the given parton bin instances
   */
  virtual void updatePartonBinInstances(const PBIPair & pbins);

  /**
   * Generate \f$l=\log(1/x)\f$ for all parton extractions.
   */
  virtual bool generateL(const PBIPair & pbins,
			 const double * r1, const double * r2);

  /**
   * Used by generateL() for each of the final parton
   * bins. Direction<0> is set to positive(negative) for the
   * first(second) final bin.
   */
  virtual void generateL(PartonBinInstance & pb, const double * r);

  /**
   * Generate the rest of the degrees of freedom to calculate
   * \f$\hat{s}\f$ and the product of all densitiy functions.
   */
  virtual Energy2 generateSHat(Energy2 s, const PBIPair & pbins,
			       const double * r1, const double * r2,
			       bool mepartons = false);

  /**
   * Return the product of all density functions. If noLastPDF.first
   * (.second) is true, then the PDF value multiplied by the momentum 
   * fraction for the last extracted parton from the first (second) 
   * incoming particle will be excluded.
   */
  virtual double fullFn(const PBIPair & pbins, Energy2 scale,
			pair<bool,bool> noLastPDF = make_pair(false,false));

  /**
   * Construct remnants and add them to the step.
   */
  virtual void construct(const PBIPair & pbins, tStepPtr step) const;

  /**
   * Construct remnants for partons created outside of this
   * extractor. Information about the incoming partons should be set
   * in \a pbins and the hard subprocess should be present in \a
   * sub. Generated remnants will be added to the \a step.
   * @throws Veto if remnant generation failed for whatever reason.
   */
  virtual void constructRemnants(const PBIPair & pbins, tSubProPtr sub,
			 tStepPtr step) const;

  /**
   * Get boost for hard subsystem and boost remnants. To be called
   * after re-constructing remnants and obtaining new momenta of the
   * partons entering the hard subsystem, but ignoring detailed energy
   * and momentum conservation. Perform boosts of the remnants to
   * conserve energy and momentum and return the boost needed for the
   * hard subsystem. \a bins contains the current state of the
   * incoming partons, including the momenta obtained after the
   * remnant generation. \a k1 and \a k2 contains the momenta of the
   * incoming partons before the remnant generation. If either side
   * has no new remnants, \a side1 and/or \a side2 should be false.
   */
  virtual LorentzRotation
  boostRemnants(PBIPair & bins, LorentzMomentum k1, LorentzMomentum k2,
	   bool side1, bool side2) const;
  //@}

  /** @name Access information about the current paron extraction. */
  //@{
  /**
   * Return the corresponding parton bin instance for a given
   * extracted parton.
   */
  tPBIPtr partonBinInstance(tcPPtr) const;

  /**
   * Set the XComb object describing the current hard sub-process.
   */
  void select(tXCombPtr newXComb);

  //@}

  /**
   * The maximum number of attempts allowed when trying to generate
   * remnants.
   */
  int maxTries() const { return theMaxTries; }

  /**
   * Return the PDFBase object to be used for the incoming particle
   * type. If one of theSpecialDensities matches the particle type it
   * is returned, otherwise if particle is a BeamParticleData use the
   * PDFBase object specified there. If also this fails, return a
   * NoPDF object.
   */
  tcPDFPtr getPDF(tcPDPtr particle) const;

protected:

  /** @name Functions used by the main virtual functions. Some of
      these may be overridden in sub-classes. */
  //@{

  /**
   * Used by generateSHat() for each of the final parton
   * bins. Direction<0> is set to positive(negative) for the
   * first(second) final bin, \a pb. Should ask the remnant handler to
   * generate what is needed to construct the extracted parton
   * momentum. \a r is a pointer to an array of random numbers to be
   * used and \a shat is the approximate invariant mass squared of the
   * hard system produced by the extracted parton and the primary
   * parton from the other side. \a first is the momentum of the
   * original incoming particle.
   *
   * @return false if no remnants could be generated.
   */
  virtual bool generate(PartonBinInstance & pb, const double * r,
			Energy2 shat, const Lorentz5Momentum & first,
			bool haveMEPartons = false);

  /**
   * Used by the public fullFn() for each of the final parton bins.
   */
  virtual double fullFn(const PartonBinInstance & pb,
			bool noLastPDF = false);

  /**
   * Used by the public construct() for each of the final parton
   * bins. If boost is false, no boost is necessary to give the
   * remnants proper momenta. 
   */
  virtual void construct(PartonBinInstance & pb,
			 tStepPtr step, bool boost = true) const;

  /**
   * Used by the public newRemnants() for each of the parton bins.
   * @throws Veto if remnant generation failed for whatever reason.
   */
  PBIPtr newRemnants(tPBIPtr oldpb, tPPtr newp, const LorentzMomentum & k);

  /**
   * Used by the public newRemnants() for each of the parton bins.
   */
  void addNewRemnants(tPBIPtr oldpb, tPBIPtr newpb, tStepPtr step);

  /**
   * Transform remnant momentum. Assuming remnants have been generated
   * with momentum \a Pr without considering energy-momentum
   * conservation, shift the momentum, possibly compensating with the
   * momentum of the hard subsystem, \a Ph. For information the
   * momentum of the parton entering the hard subsystem from the other
   * side, \a k, and the momentum of the remnants parent particle , \a
   * P is given. Note that Direction<0> must be set to determine if
   * the parent particle is to be assumed to go in the positive or
   * negative direction.
   */
  virtual void transformRemnants(LorentzMomentum & Ph, LorentzMomentum & Pr,
				 const LorentzMomentum & k,
				 const LorentzMomentum & P) const;

  /**
   * Construct remnants recursively for the parton represented by \a
   * pb. Used by constructRemnants(const PBIPair &, tSubProPtr, tStepPtr).
   * Shift the momentum, \a Ph, of the hard subsystem to conserve
   * energy and momentum if necessary. The momentum, \a k, of the
   * parton coming into the hard subsystem from the other side is
   * given for information. Note that Direction<0> must be set to
   * determine if the parent particle is to be assumed to go in the
   * positive or negative direction.
   * @throws Veto if remnant generation failed for whatever reason.
   */
  virtual void
  constructRemnants(PartonBinInstance & pb, LorentzMomentum & Ph,
		    const LorentzMomentum & k) const;
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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /**
   * Add parton bins to pbins for the given incoming particle and the
   * specified cuts.
   */
  virtual void addPartons(tPBPtr incoming ,const PDFCuts & cuts,
			  tcPDFPtr pdf ,PartonVector & pbins) const;

  /**
   * The NoPDF object.
   */
  tcPDFPtr noPDF() const { return theNoPDF; }

  /**
   * Connect the first (\a anti) coloured particle in the given range
   * (not equal to \a parton) and connect it to the colour \a line.
   */
  template <typename Iterator>
  void findConnect(tColinePtr line, tPPtr parton, bool anti,
		   Iterator first, Iterator last) const {
    for ( ; first != last; ++first ) {
      if ( *first != parton &&  (**first).hasColour(anti) &&
	   !(**first).colourLine(anti) ) {
	line->addColoured(*first, anti);
	return;
      }
    }
    throw RemColException(*this);
  }

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

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The PartonBinInstance's used mapped to the respective partons.
   */
  PartonBinInstanceMap& partonBinInstances() const { 
    assert(lastXCombPtr());
    return lastXCombPtr()->partonBinInstanceMap();
  }

  /**
   * A list of special PDFBase objects to be used.
   */
  vector<PDFPtr> theSpecialDensities;

  /**
   *  PDFBase object to override  first PDF
   */
  PDFPtr theFirstPDF;

  /**
   *  PDFBase object to override second PDF
   */
  PDFPtr theSecondPDF;

  /**
   * The NoPDF object.
   */
  PDFPtr theNoPDF;

  /**
   * The maximum number of tries allowed when trying to produce
   * remnants.
   */
  int theMaxTries;

  /**
   * True if this extractor should override the \f$l\f$-generation in
   * the PDFs and generate a flat distribution in \f$\log(\hat{s})\f$
   * and y.
   */
  bool flatSHatY;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<PartonExtractor> initPartonExtractor;

  /**
   *  Private and non-existent assignment operator.
   */
  PartonExtractor & operator=(const PartonExtractor &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of PartonExtractor. */
template <>
struct BaseClassTrait<PartonExtractor,1>: public ClassTraitsType {
  /** Typedef of the first base class of PartonExtractor. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  PartonExtractor class. */
template <>
  /** Return a platform-independent class name */
struct ClassTraits<PartonExtractor>: public ClassTraitsBase<PartonExtractor> {
  static string className() { return "ThePEG::PartonExtractor"; }
};

/** @endcond */

}

#endif /* ThePEG_PartonExtractor_H */
