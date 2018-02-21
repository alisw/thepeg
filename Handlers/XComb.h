// -*- C++ -*-
//
// XComb.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_XComb_H
#define ThePEG_XComb_H
// This is the declaration of the XComb class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDF/PartonExtractor.fh"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/Utilities/AnyReference.h"
#include "ThePEG/Utilities/VSelector.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Handlers/EventHandler.fh"
#include "ThePEG/Cuts/Cuts.fh"

namespace ThePEG {

/**
 * The XComb class stores all information about the generation of a
 * hard sub-proces for a given pair of incoming particles, a pair of
 * extracted partons, total parton-parton energy squared and a
 * PartonExtractor object.
 *
 * When an event is generated, the objects used in the generation can
 * be assigned an XComb object for easy acces to the corresponding
 * information. To facilitate this, the corresponding classes inherits
 * from the LastXCombInfo class which provides the relefant access
 * functions.
 *
 * @see PartonExtractor
 * @see Cuts
 * @see LastXCombInfo
 */
class XComb: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Standard constructor.
   */
  XComb(Energy newMaxEnergy, const cPDPair & inc,
	tEHPtr newEventHandler,	tPExtrPtr newExtractor, tCascHdlPtr newCKKW,
	const PBPair & newPartonBins,	tCutsPtr newCuts);

  /**
   * Default constructor.
   */
  XComb();

  /**
   * Destructor.
   */
  virtual ~XComb();
  //@}



  /** @name Access the assigned objects used in the generation. */
  //@{
  /**
   * Return a reference to the corresponding collision handler
   */
  const EventHandler & eventHandler() const { return *theEventHandler; }

  /**
   * Return a pointer to the corresponding collision handler
   */
  tEHPtr eventHandlerPtr() const { return theEventHandler; }

  /**
   * A pointer to the parton extractor.
   */
  tPExtrPtr pExtractor() const { return thePartonExtractor; }

  /**
   * A pointer to the kinematical cuts.
   */
  tCutsPtr cuts() const { return theCuts; }

  /**
   * Return a possibly null pointer to a CascadeHandler to be used for
   * CKKW-reweighting.
   */
  tCascHdlPtr CKKWHandler() const { return theCKKW; }
  //@}

  /** @name Access information about incoming particles and partons. */
  //@{
  /**
   * The incoming particle types.
   */
  const cPDPair & particles() const { return theParticles; }

  /**
   * The incoming parton types.
   */
  const cPDPair & partons() const { return thePartons; }

  /**
   * Additional information about the incoming partons.
   */
  const PBPair & partonBins() const { return thePartonBins; }
  
  /**
   * The maximum cm energy for this process.
   */
  Energy maxEnergy() const { return theMaxEnergy; }

  /**
   * Returns true if this XComb does not correspond to a proper
   * subprocess generation. I.e. if we are only generating a partial
   * event and the incoming particles and partons are not used
   * explicitly.
   */
  bool empty() const { return !theEventHandler; }
  //@}

  /** @name Manipulate and acces information about the last selected
      phase space point. */
  //@{

  /**
   * Reset all saved data about last generated phasespace point;
   */
  virtual void clean();

  /**
   * Set information about currently generated partons.
   */
  void setPartonBinInstances(PBIPair pbis, Energy2 scale);

  /**
   * Prepare this XComb for producing a sub-process.
   */
  void prepare(const PPair &);

  /**
   * Return the pair of incoming particle instances.
   */
  const PPair & lastParticles() const { return theLastParticles; }

  /**
   * Return the pair of incoming parton instances.
   */
  const PPair & lastPartons() const { return theLastPartons; }

  /**
   * Set the pair of incoming parton instances.
   */
  void lastPartons(PPair pp) { theLastPartons = pp; }

  /**
   * Return the SubProcess object corresponding to the last generated
   * sub-process.
   */
  tSubProPtr subProcess() const { return theSub; }

  /**
   * Set the SubProcess object corresponding to the last generated
   * sub-process.
   */
  void subProcess(tSubProPtr);

  /** A map of PartonBinInstance objects indexed by the extracted parton. */
  typedef map<cPPtr,PBIPtr> PartonBinInstanceMap;

  /**
   * Access the parton bin instance map (used by the parton extractor)
   */
  PartonBinInstanceMap& partonBinInstanceMap() { return thePartonBinInstanceMap; }

  /**
   * Return the parton bin instance map (used by the parton extractor)
   */
  const PartonBinInstanceMap& partonBinInstanceMap() const { return thePartonBinInstanceMap; }

  /**
   * Additional information about the incoming partons.
   */
  const PBIPair & partonBinInstances() const { return thePartonBinInstances; }

  /**
   * Additional information about the incoming partons.
   */
  PBIPair & partonBinInstances() { return thePartonBinInstances; }

  /**
   * Return the corresponding parton bin instance for a given
   * extracted parton.
   */
  tPBIPtr partonBinInstance(tcPPtr) const;

  /**
   * The last generated total energy squared of the incoming particles.
   */
  Energy2 lastS() const { return theLastS; }

  /**
   * Set the last generated total energy squared of the incoming
   * particles.
   */
  void lastS(Energy2 s) { theLastS = s; }

  /**
   * The last generated total energy squared of the incoming prtons.
   */
  Energy2 lastSHat() const { return theLastSHat; }

  /**
   * Set the last generated total energy squared of the incoming
   * prtons.
   */
  void lastSHat(Energy2 sh) { theLastSHat = sh; }

  /**
   * lastSHat()/lastS().
   */
  double lastTau() const { return lastSHat()/lastS(); }

  /**
   * The last generated rapidity of the hard scattering sub-system.
   */
  double lastY() const { return theLastY; }

  /**
   * Set the last generated rapidity of the hard scattering sub-system.
   */
  void lastY(double y) { theLastY = y; }

  /**
   * Log of one over the momentum fraction of the first incoming
   * particle w.r.t. the maximum allowed energy.
   */
  double lastP1() const { return theLastP1P2.first; }

  /**
   * Log of one over the momentum fraction of the second incoming
   * particle w.r.t. the maximum allowed energy.
   */
  double lastP2() const { return theLastP1P2.second; }

  /**
   * Set log of one over the momentum fraction of the incoming
   * particles w.r.t. the maximum allowed energy.
   */
  void lastP1P2(pair<double,double> pp) { theLastP1P2 = pp; }

  /**
   * Log of one over the first incoming parton momentum fraction
   * w.r.t. the first incoming particle.
   */
  double lastL1() const { return theLastL1L2.first; }

  /**
   * Log of one over the second incoming parton momentum fraction
   * w.r.t. the second incoming particle.
   */
  double lastL2() const { return theLastL1L2.second; }

  /**
   * Set log of one over the incoming parton momentum fractions
   * w.r.t. the incoming particles.
   */
  void lastL1L2(pair<double,double>);

  /**
   * The first incoming parton momentum fraction w.r.t. the
   * first incoming particle.
   */
  double lastX1() const { return theLastX1X2.first; }

  /**
   * The second incoming parton momentum fraction
   * w.r.t. the second incoming particle.
   */
  double lastX2() const { return theLastX1X2.second; }

  /**
   * Set the incoming parton momentum fractions w.r.t. the incoming
   * particles.
   */
  void lastX1X2(pair<double,double>);

  /**
   * Return 1-lastX1() to highest possible precision for
   * x\f$\rightarrow\f$ 1.
   */
  double lastE1() const { return theLastE1E2.first; }

  /**
   * Return 1-lastX2() to highest possible precision for
   * x\f$\rightarrow\f$ 1.
   */
  double lastE2() const { return theLastE1E2.second; }

  /**
   * Set one minus the incoming parton momentum fractions w.r.t. the
   * incoming particles.
   */
  void lastE1E2(pair<double,double>);

  /**
   * Get the last chosen scale of the hard scattering.
   */
  Energy2 lastScale() const { return theLastScale; }

  /**
   * Set the last chosen scale of the hard scattering.
   */
  void lastScale(Energy2 Q2) { theLastScale = Q2; }

  /**
   * Get the last chosen central scale of the hard scattering.
   */
  Energy2 lastCentralScale() const { 
    return 
      theLastCentralScale != ZERO ?
      theLastCentralScale :
      lastScale(); 
  }

  /**
   * Set the last chosen central scale of the hard scattering.
   */
  void lastCentralScale(Energy2 Q2) { theLastCentralScale = Q2; }

  /**
   * Get the last chosen shower scale.
   */
  Energy2 lastShowerScale() const { 
    return 
      theLastShowerScale != ZERO ?
      theLastShowerScale :
      lastCentralScale(); 
  }

  /**
   * Set the last chosen showr scale.
   */
  void lastShowerScale(Energy2 Q2) { theLastShowerScale = Q2; }

  /**
   * Get the \f$\alpha_S\f$ used in the hard scattering. Is negative
   * if no value has been set.
   */
  double lastAlphaS() const { return theLastAlphaS; }

  /**
   * Set the \f$\alpha_S\f$ used in the hard scattering.
   */
  void lastAlphaS(double a) { theLastAlphaS = a; }

  /**
   * Get the \f$\alpha_{EM}\f$ used in the hard scattering. Is negative
   * if no value has been set.
   */
  double lastAlphaEM() const { return theLastAlphaEM; }

  /**
   * Set the \f$\alpha_{EM}\f$ used in the hard scattering.
   */
  void lastAlphaEM(double a) { theLastAlphaEM = a; }
  //@}

public:

  /**
   * Check for meta information
   */
  bool hasMeta(int id) const {
    return theMeta.find(id) != theMeta.end();
  }

  /**
   * Set meta information.
   */
  template<class T>
  void meta(int id, T& ref) {
    theMeta[id] = AnyReference(ref);
  }

  /**
   * Erase meta information.
   */
  void eraseMeta(int id) {
    theMeta.erase(id);
  }

  /**
   * Retrieve meta information.
   */
  template<class T>
  T& meta(int id) const {
    return theMeta.find(id)->second.cast<T>();
  }

  /**
   * Set the local parton bin info objects for this XComb.
   */
  void setPartonBinInfo();

  /**
   * Create PartonBinInstance objects for this XComb.
   */
  void createPartonBinInstances();

  /**
   * Set the pair of incoming particle instances.
   */
  void lastParticles(const PPair & p) { theLastParticles = p; }

  /**
   * Set information about currently generated partons.
   */
  void resetPartonBinInstances(const PBIPair & newBins) { thePartonBinInstances = newBins; }

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
   * The corresponding collision handler
   */
  tEHPtr theEventHandler;

  /**
   * A pointer to the parton extractor.
   */
  tPExtrPtr thePartonExtractor;

  /**
   * A pointer to a CascadeHandler to be used for CKKW-reweighting.
   */
  tCascHdlPtr theCKKW;

  /**
   * A pointer to the kinematical cuts used.
   */
  tCutsPtr theCuts;

  /**
   * The incoming particle types.
   */
  cPDPair theParticles;

  /**
   * The incoming parton types.
   */
  cPDPair thePartons;

  /**
   * The parton bin instance map (used by the parton extractor)
   */
  PartonBinInstanceMap thePartonBinInstanceMap;

  /**
   * Additional information about the incoming partons.
   */
  PBPair thePartonBins;

  /**
   * Additional information about the origins of the incoming partons.
   */
  PBPair theParticleBins;

  /**
   * Additional information about the incoming partons.
   */
  PBIPair thePartonBinInstances;

  /**
   * The pair of incoming particle instances.
   */
  PPair theLastParticles;

  /**
   * The pair of incoming parton instances.
   */
  PPair theLastPartons;

  /**
   * The last generated total energy squared of the incoming particles.
   */
  Energy2 theLastS;

  /**
   * The last generated total energy squared of the incoming prtons.
   */
  Energy2 theLastSHat;

  /**
   * The last rapidity of the sub process, log(x1/x2)/2.
   */
  double theLastY;

  /**
   * Log of one over the momentum fraction of the incoming particles.
   */
  DPair theLastP1P2;

  /**
   * Log of one over the incoming partons momentum fraction wrt. the
   * incoming particles.
   */
  DPair theLastL1L2;

  /**
   * The incoming partons momentum fraction wrt. the incoming
   * particles.
   */
  DPair theLastX1X2;

  /**
   * 1-lastX1() and 1-lastX2() to highest possible precision for
   * x\f$\rightarrow\f$ 1.
   */
  DPair theLastE1E2;

  /**
   * The last chosen scale of the hard scattering.
   */
  Energy2 theLastScale;

  /**
   * The last chosen central scale of the hard scattering.
   */
  Energy2 theLastCentralScale;

  /**
   * The last chosen shower scale.
   */
  Energy2 theLastShowerScale;

  /**
   * The \f$\alpha_S\f$ used in the hard scattering.
   */
  double theLastAlphaS;

  /**
   * The \f$\alpha_{EM}\f$ used in the hard scattering.
   */
  double theLastAlphaEM;

  /**
   * The maximum cm energy for this process.
   */
  Energy theMaxEnergy;

  /**
   * Information saved by the matrix element in the calculation of the
   * cross section to be used later when selecting diagrams and colour
   * flow.
   */
  DVector theMEInfo;

  /**
   * The SubProcess object corresponding to the last generated
   * sub-process.
   */
  SubProPtr theSub;

  /**
   * The meta information
   */
  map<int,AnyReference> theMeta;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<XComb> initXComb;
 
  /**
   * Private and non-existent assignment operator.
   */
  XComb & operator=(const XComb &);

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * XComb.
 */
template <>
struct BaseClassTrait<XComb,1>: public ClassTraitsType {
  /** Typedef of the base class of XComb. */
  typedef Base NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * XComb class.
 */
template <>
struct ClassTraits<XComb>:
    public ClassTraitsBase<XComb> {
  /** Return the class name. */
  static string className() { return "ThePEG::XComb"; }
};

/** @endcond */

}

#endif /* ThePEG_XComb_H */
