// -*- C++ -*-
//
// LastXCombInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LastXCombInfo_H
#define ThePEG_LastXCombInfo_H
// This is the declaration of the LastXCombInfo class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Cuts/Cuts.fh"
#include "XComb.h"

namespace ThePEG {

/**
 * LastXCombInfo is a templated class giving easy access to the
 * information in an XComb object. The default template argument is
 * the basic XComb class, but also subclasses of XComb can be
 * used. Classes which need to have easy access to the last selected
 * XComb object with information about the sub-process which is being
 * generated, should (possibly multiple) inherit from the
 * LastXCombInfo class. The LastXCombInfo is templated to enable
 * derived classes to only include dependencies necessary for the
 * access function which are actually used.
 * 
 */
template <typename XC = XComb>
class LastXCombInfo {

public:

  ThePEG_DECLARE_TEMPLATE_POINTERS(XC,XCPtr);

public:

  /** @name Acces to the actual XComb object. */
  //@{
  /**
   * Return a reference to the last selected XComb.
   */
  const XC & lastXComb() const { return *theLastXComb; }

  /**
   * Return a pointer to the last selected XComb.
   */
  tXCPtr lastXCombPtr() const { return theLastXComb; }

  /**
   * If the last selected XComb object belongs to a
   * group of XComb's return a reference to the head 
   * XComb object for this group.
   */
  const XC& lastHeadXComb() const { return *lastXComb().head(); }

  /**
   * If the last selected XComb object belongs to a
   * group of XComb's return a pointer to the head 
   * XComb object for this group.
   */
  tXCPtr lastHeadXCombPtr() const { return lastXComb().head(); }
  //@}

  /** @name Access the objects used by the XComb object. */
  //@{
  /**
   * Return a reference to the currently used EventHandler
   */
  const EventHandler & lastEventHandler() const { return lastXComb().eventHandler(); }

  /**
   * A pointer to the currently used parton extractor.
   */
  tPExtrPtr lastExtractor() const { return lastXComb().pExtractor(); }

  /**
   * Return the parton density used to extract the given parton. This
   * function is templated to avoid having to include the PDF.h and
   * all its dependencies in this header.
   */
  template <typename PDFT>
  PDFT pdf(tcPPtr parton) const {
    return PDFT(lastXComb().partonBinInstance(parton));
  }

  /**
   * A reference to the currently used kinematical cuts.
   */
  const Cuts & lastCuts() const { return *lastXComb().cuts(); }

  /**
   * A pointer to the currently used kinematical cuts.
   */
  tCutsPtr lastCutsPtr() const { return lastXComb().cuts(); }

  //@}

  /** @name Access information about the incoming particles and partons. */
  //@{
  /**
   * Return the pair of incoming parton instances.
   */
  const PPair & lastParticles() const { return lastXComb().lastParticles(); }

  /**
   * The last generated total energy squared of the incoming particles.
   */
  Energy2 lastS() const { return lastXComb().lastS(); }

  /**
   * Return the pair of incoming parton instances.
   */
  const PPair & lastPartons() const { return lastXComb().lastPartons(); }

  /**
   * The last used interval in total parton-parton energy squared
   */
  Energy2 lastSHat() const { return lastXComb().lastSHat(); }

  /**
   * Return lastSHat()/lastS().
   */
  double lastTau() const { return lastXComb().lastTau(); }

  /**
   * The generated rapidity of the hard scattering sub-system.
   */
  double lastY() const { return lastXComb().lastY(); }

  /**
   * Log of one over the momentum fraction of the first incoming
   * particle w.r.t. the maximum allowed energy.
   */
  double lastP1() const { return lastXComb().lastP1(); }

  /**
   * Log of one over the momentum fraction of the second incoming
   * particle w.r.t. the maximum allowed energy.
   */
  double lastP2() const { return lastXComb().lastP2(); }

  /**
   * Log of one over the first incoming parton momentum fraction w.r.t. the
   * first incoming particle.
   */
  double lastL1() const { return lastXComb().lastL1(); }

  /**
   * Log of one over the second incoming parton momentum fraction
   * w.r.t. the second incoming particle.
   */
  double lastL2() const { return lastXComb().lastL2(); }

  /**
   * The first incoming parton momentum fraction w.r.t. the
   * first incoming particle.
   */
  double lastX1() const { return lastXComb().lastX1(); }

  /**
   * The second incoming parton momentum fraction
   * w.r.t. the second incoming particle.
   */
  double lastX2() const { return lastXComb().lastX2(); }

  /**
   * Return 1-lastX1() to highest possible precision for
   * x \f$\rightarrow\f$ 1.
   */
  double lastE1() const { return lastXComb().lastE1(); }

  /**
   * Return 1-lastX2() to highest possible precision for
   * x\f$\rightarrow\f$ 1.
   */
  double lastE2() const { return lastXComb().lastE2(); }

  /**
   * The product of the parton density functions at the last generated
   * phase-space point.
   */
  double lastFL1L2() const { return lastXComb().lastFL1L2(); }
  //@}

  /** @name Access information of the hard sub-process. */
  //@{
  /**
   * The chosen scale of the hard scattering.
   */
  Energy2 lastScale() const { return lastXComb().lastScale(); }

  /**
   * Get the \f$\alpha_S\f$ used in the hard scattering. Is negative
   * if no value has been set.
   */
  double lastAlphaS() const { return lastXComb().lastAlphaS(); }

  /**
   * Get the \f$\alpha_{EM}\f$ used in the hard scattering. Is negative
   * if no value has been set.
   */
  double lastAlphaEM() const { return lastXComb().lastAlphaEM(); }

  /**
   * Return the momenta of the incoming and outgoing partons to be
   * used by the matrix element object, in the order specified by the
   * TreeDiagram objects given by the matrix element.
   */
  const vector<Lorentz5Momentum> & meMomenta() const { return lastXComb().meMomenta(); }

  /**
   * Return the matrix element squared as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  double lastME2() const { return lastXComb().lastME2(); }

  /**
   * Return the last preweight factor
   */
  double lastPreweight() const { return lastXComb().lastPreweight(); }

  /**
   * Get the last jacobian obtained when generating the kinematics
   * for the call to dSigHatDR.
   */
  double jacobian() const { return lastXComb().jacobian(); }

  /**
   * Return the partonic cross section as calculated
   * for the last phase space point. This may optionally
   * be used by a matrix element for caching.
   */
  CrossSection lastMECrossSection() const { return lastXComb().lastMECrossSection(); }

  /**
   * Return the PDF weight as calculated
   * for the last phase space point, if the matrix
   * element does supply PDF weights. This may optionally
   * be used by a matrix element for caching.
   */
  double lastMEPDFWeight() const { return lastXComb().lastMEPDFWeight(); }

  /**
   * Return the coupling weight as calculated
   * for the last phase space point.
   */
  double lastMECouplings() const { return lastXComb().lastMECouplings(); }

  /**
   * Return the SubProcess object corresponding to the last generated
   * sub-process.
   */
  tSubProPtr subProcess() const { return lastXComb().subProcess(); }

  /**
   * Return the incoming and outgoing parton types to be used by the
   * matrix element object, in the order specified by the TreeDiagram
   * objects given by the matrix element.
   */
  const cPDVector & mePartonData() const { return lastXComb().mePartonData(); }
  //@}

protected:

  /**
   * The pointer to the last selected XComb.
   */
  XCPtr theLastXComb;
};

}

#endif /* ThePEG_LastXCombInfo_H */
