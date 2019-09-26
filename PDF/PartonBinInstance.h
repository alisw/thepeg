// -*- C++ -*-
//
// PartonBinInstance.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_PartonBinInstance_H
#define THEPEG_PartonBinInstance_H
// This is the declaration of the PartonBinInstance class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDF/PartonBin.h"

namespace ThePEG {

ThePEG_DECLARE_CLASS_POINTERS(PartonBinInstance,PBIPtr);
/** A pair of pointers to PartonBinInstance objects. */
typedef pair<PBIPtr,PBIPtr> PBIPair;

ThePEG_DECLARE_CLASS_POINTERS(RemInfoBase,RemIPtr);

/**
 * PartonBinInstance is used to store information about the generation
 * of a given parton extraction for a corresponding PartonBin object.
 */
class PartonBinInstance: public PersistentBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PartonBinInstance();

  /**
   * Copy-constructor.
   */
  PartonBinInstance(const PartonBinInstance &);

  /**
   * Destructor.
   */
  virtual ~PartonBinInstance();

  /**
   * Constructor taking a PartonBin as argument. The second argument
   * should be used if the incoming bin is already known and exists.
   */
  PartonBinInstance(tcPBPtr, tPBIPtr = tPBIPtr());

  /**
   * Constructor using an already prepared extracted parton. This will
   * also initialize the x, and scale values. To calculate the
   * momentum fractions, a Direction<0> object must have been properly
   * initialized.
   *
   * @param parton the extracted parton which must have its first
   * parent set to define the particle extracted from.
   *
   * @param pb the PartonBin object corresponding to the extracted \a
   * parton. If the particle extracted from in turn has been
   * extracted, the incoming() member of the PartonBin must point to
   * the corresponding PartonBin.
   *
   * @param scale the resolution scale at which the \a parton was
   * extracted.
   */
  PartonBinInstance(tPPtr parton, tcPBPtr pb, Energy2 scale = ZERO);

  /**
   * Constructor using a parton which is to be extracted from the
   * given particle, but no mother-child relations exist, yet. This
   * will also initialize the x, and scale values. To calculate the
   * momentum fractions, a Direction<0> object must have been properly
   * initialized.
   */
  PartonBinInstance(tPPtr particle, tPPtr parton, tcPBPtr pb, 
		    Energy2 scale = ZERO);

  //@}

public:

  /** @name Access information about the corresponding PartonBin object. */
  //@{
  /**
   * Return a pointer to the PartonBin this instance refer to.
   */
  tcPBPtr bin() const { return theBin; }

  /**
   * Return pointers to the bins this instance refer to in case more
   * than one parton has been extracted.
   */
  const PartonVector & bins() const { return theBins; }

  /**
   * Return a pointer to the data object of the incoming particle.
   */
  tcPDPtr particleData() const { return bin()->particle(); }

  /**
   * Return a pointer to the data object of the extracted parton.
   */
  tcPDPtr partonData() const { return bin()->parton(); }

  /**
   * In the case the incoming particle in turn is extracted from
   * another particle, return the PartonBinInstance for that
   * extraction.
   */
  tPBIPtr incoming() const { return theIncoming; }

  /**
   * Return the parton bin instance corresponding to the first
   * incoming particle for this bin.
   */
  tPBIPtr getFirst();

  /**
   * The PDFBase object describing the momentum distribution of the
   * parton within the particle in this PartonBin.
   */
  tcPDFPtr pdf() const { return bin()->pdf(); }

  /**
   * The remnant handler associated with the pdf().
   */
  tcRemHPtr remnantHandler() const { return bin()->remnantHandler(); }

  /**
   * Return true if the corresponding PDFs has a pole at $x=1$ for the
   * current particle/parton combination.
   */
  bool hasPoleIn1() const;
  //@}

  /** @name Functions used for the generation. */
  //@{
  /**
   * Reset the current PartonBin, making room for a new event.
   */
  void reset(double lx = 0, Energy2 Q2 = ZERO);

  /**
   * Reset last generated l and Q2 values of this and parent bins.
   */
  void prepare();

  /**
   * Generate l and Q2 of this and parent bins.
   */
  void generate(const double * r);

  /**
   * Get the jacobian associated with the phase space point generated.
   */
  double jacobian() const { return theJacobian; }

  /**
   * Set the jacobian associated with the phase space point generated.
   */
  void jacobian(double j) { theJacobian = j; }
  //@}

  /** @name Access information about the generated extraction. */
  //@{
  /**
   * Get the current particle instance.
   */
  tPPtr particle() const { return theParticle; }

  /**
   * Set the current particle instance.
   */
  void particle(tPPtr p) { theParticle = p; }

  /**
   * Get the current parton instance.
   */
  tPPtr parton() const { return theParton; }

  /**
   * Set the current parton instance.
   */
  void parton(tPPtr p) { theParton = p; }

  /**
   * The currently extracted partons (in case of multiple
   * interactions.
   */
  const PVector & partons() const { return thePartons; }

  /**
   * Get the momentum fraction of this parton w.r.t. the incoming
   * particle in this bin.
   */
  double xi() const {
    if ( theXi < 0.0 ) theXi = exp(-li());
    return theXi;
  }


  /**
   * Get one minus the momentum fraction of this parton w.r.t. the
   * incoming particle in this bin.
   */
  double eps() const {
    if ( theEps < 0.0 ) theEps =  Math::exp1m(-li());
    return theEps;
  }

  /**
   * Get the logarithmic momentum fraction of this parton w.r.t. the
   * incoming particle in this bin.
   */
  double li() const { return theLi; }

  /**
   * Set the logarithmic momentum fraction of this parton w.r.t. the
   * incoming particle in this bin.
   */
  void li(double lx) {
    theLi = lx;
    theXi = theEps = -1.0;
  }


  /**
   * Get the momentum fraction of this parton w.r.t. the collidig
   * particles.
   */
  double x() const {
    if ( theX < 0.0 ) theX = exp(-l());
    return theX;
  }


  /**
   * Get the logarithmic momentum fraction of this parton w.r.t. the
   * collidig particles.
   */
  double l() const { return theL; }

  /**
   * Set the logarithmic momentum fraction of this parton w.r.t. the
   * collidig particles.
   */
  void l(double lx) {
    theL = lx;
    theX = -1.0;
  }


  /**
   * Get the scale at which the current parton was extracted.
   */
  Energy2 scale() const { return theScale; }
  

  /**
   * Set the scale at which the current parton was extracted.
   */
  void scale(Energy2 s) { theScale = s; }

  /**
   * Return the transverse momentum of the extracted parton.
   */
  const TransverseMomentum & kT() const { return theKT; }

  /**
   * Get the weight associated with the remnant generation.
   */
  double remnantWeight() const { return theRemnantWeight; }

  /**
   * Set the weight associated with the remnant generation.
   */
  void remnantWeight(double w) { theRemnantWeight = w; }

  /**
   * Get the current remnants.
   */
  const PVector & remnants() const { return theRemnants; }

  /**
   * Set the current remnants.
   */
  void remnants(const PVector & rems) { theRemnants = rems; }

  /**
   * Get information saved by the remnant handler from the generation,
   * to be used in the construction of the remnants. (In addition the
   * remnantWeight and remnants() may be used for this purpose.)
   */
  tRemIPtr remnantInfo() const { return theRemInfo; }

  /**
   * Set information saved by the remnant handler from the generation,
   * to be used in the construction of the remnants. (In addition the
   * remnantWeight and remnants() may be used for this purpose.)
   */
  void remnantInfo(tRemIPtr ri) { theRemInfo = ri; }
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

private:

  /**
   * Pointer to the main bin this instance refer to.
   */
  cPBPtr theBin;

  /**
   * Pointer to the main bin (and secondary in case several partons
   * have been extracted this instance refer to.
   */
  PartonVector theBins;

  /**
   * In the case the incoming particle in turn is extracted from
   * another particle, this is the PartonBinInstance for that
   * extraction.
   */
  PBIPtr theIncoming;

  /**
   * The jacobian associated with the phase space point generated.
   */
  double theJacobian;

  /**
   * The current particle instance.
   */
  PPtr theParticle;

  /**
   * The current parton instance.
   */
  PPtr theParton;

  /**
   * The currently extracted partons (in case of multiple
   * interactions.
   */
  PVector thePartons;

  /**
   * The momentum fraction (xi, li=log(xi), eps=1-xi), of this
   * parton w.r.t. the incoming particle in this
   * bin.
   */
  mutable double theXi;
  /**
   * The momentum fraction (xi, li=log(xi), eps=1-xi), of this
   * parton w.r.t. the incoming particle in this
   * bin.
   */
  mutable double theEps;
  /**
   * The momentum fraction (xi, li=log(xi), eps=1-xi), of this
   * parton w.r.t. the incoming particle in this
   * bin.
   */
  double theLi;

  /**
   * The momentum fraction (x, l=log(x)) of this parton
   * w.r.t. the collidig particles.
   */
  mutable double theX;
  /**
   * The momentum fraction (x, l=log(x)) of this parton
   * w.r.t. the collidig particles.
   */
  double theL;

  /**
   * The scale at which the current parton was extracted.
   */
  Energy2 theScale;

  /**
   * The transverse momentum of the extracted parton.
   */
  TransverseMomentum theKT;

  /**
   * The weight associated with the remnant generation.
   */
  double theRemnantWeight;

  /**
   * The current remnants.
   */
  PVector theRemnants;

  /**
   * The information saved by the remnant handler from the generation,
   * to be used in the construction of the remnants. (In addition the
   * remnantWeight and lastRemnants() may be used for this purpose.)
   */
  RemIPtr theRemInfo;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<PartonBinInstance> initPartonBinInstance;

  /**
   * Private and non-existent assignment operator.
   */
  PartonBinInstance & operator=(const PartonBinInstance &) = delete;

};

/** Empty base class. A RemnantHandler may use sub-classes to store
    information about the generation of remnants. */
class RemInfoBase: public Base {
public:
  /** The descructor. */
  virtual ~RemInfoBase() {}
};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of PartonBinInstance. */
template <>
struct BaseClassTrait<PartonBinInstance,1>: public ClassTraitsType {
  /** Typedef of the first base class of PartonBinInstance. */
  typedef Base NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  PartonBinInstance class. */
template <>
struct ClassTraits<PartonBinInstance>:
    public ClassTraitsBase<PartonBinInstance> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::PartonBinInstance"; }
};

/** @endcond */

}

#endif /* THEPEG_PartonBinInstance_H */
