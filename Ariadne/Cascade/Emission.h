// -*- C++ -*-
#ifndef Ariadne5_Emission_H
#define Ariadne5_Emission_H
//
// This is the declaration of the Emission class.
//

#include "Ariadne/Config/Ariadne5.h"
#include "Ariadne/Config/CloneBase.h"
#include "Emission.fh"
#include "EmitterBase.fh"
#include "DipoleBase.fh"
#include "Parton.fh"
#include "ThePEG/Utilities/Triplet.h"

namespace Ariadne5 {

/**
 * Emission is the base class to be used by all sub-classes of
 * EmissionModel to specify a generated emission for an Emitter.
 */
class Emission: public CloneBase {

public:

  /** Convenient typedef. */
  typedef Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum> Trip;

  /**
   * Enum giving the state of this Emission.
   */
  enum State {
    generated, /**< This Emission has been created but not yet performed. */
    performed, /**< This Emission has been performed. */
    reverted   /**< This Emission has been performed but then reverted. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor.
   */
  Emission(const EmitterBase & inmodel, const DipoleBase & indipole)
    : state(generated), emno(0), geno(0), model(&inmodel), cdipole(&indipole),
      rho(ZERO), y(0.0), ymax(-1.0),
      orderAlphaS(1), orderAlphaEW(0), prob(-1.0), weightPDF(-1.0),
      reversible(true), failsafe(false) {}

  /**
   * The deault constructor should not normally be used.
   */
  Emission():  state(generated), emno(0), geno(0), rho(ZERO), y(0.0), ymax(-1.0),
	       orderAlphaS(0), orderAlphaEW(0), prob(-1.0), weightPDF(-1.0),
	       reversible(false), failsafe(false) {}

  /**
   * The destructor.
   */
  virtual ~Emission() {}
  //@}

public:

  /** @name Functions relating to the DipoleState and CascadeHandler
   *  to which this belongs. */
  //@{
  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

public:

  /**
   * Perform this emission.
   */
  bool perform() const;

  /**
   * Revert this emission if it for some reason failed.
   */
  void revert() const;

  /**
   * The State of this Emission
   */
  mutable State state;

  /**
   * The emission number
   */
  mutable int emno;

  /**
   * The generation number
   */
  mutable int geno;

  /**
   * The emission model.
   */
  tcEmitterPtr model;

  /**
   * The unmutable dipole, available during the generation.
   */
  tcDBPtr cdipole;

  /**
   * The dipole. Is only set after the full emission has been generated.
   */
  tDBPtr dipole;

  /**
   * The momentum of the original partons.
   */
  mutable pair<Lorentz5Momentum,Lorentz5Momentum> pold;

  /**
   * The evolution scale of the generated emission.
   */
  Energy rho;

  /**
   *The invariant rapidity of the generated emission.
   */
  double y;

  /**
   * The generated mometa.
   */
  mutable Trip genmom;

  /**
   * The approximate maximum invariant rapidity possible for the
   * given scale (typically given by log(S/pt^2)/2. Is negative if
   * no such value is given.
   */
  double ymax;

  /**
   * The order in \f$\alpha_S\f$ of this emission. If this is a
   * pre-generated emission (e.g. a decay given already in the hard
   * sub-process) it should be set to the negative.
   */
  int orderAlphaS;

  /**
   * The order in \f$\alpha_{EW}\f$ of this emission. If this is
   * a pre-generated emission (e.g. a decay given already in the
   * hard sub-process) it shoule be set to the negative.
   */
  int orderAlphaEW;

  /**
   * The partons responsible for the emission (if any). Typically the
   * two partons spanning a dipole.
   */
  vector<tParPtr> radiators;

  /**
   * The parent which was the incoming (coloured) parton in the
   * original dipole.
   */
  tParPtr colourParent;
  
  /**
   * The parent which was the outgoing (anti-coloured) parton in the
   * original dipole.
   */
  tParPtr antiColourParent;

  /**
   * The parent from which the main particle was most likely emitted.
   */
  tParPtr mainParent;

  /**
   * The list of emitted partons
   */
  mutable vector<tParPtr> partons;

  /**
   * The list of partons which are affected by this emission.
   */
  mutable vector<tParPtr> affected;

  /**
   * The splitting probability for this (inverse) emission given in
   * units of 1/GeV.
   */
  double prob;

  /**
   * The PDF ratio weight for this (inverse) emission.
   */
  double weightPDF;

  /**
   * Can the emitter model revert this Emission so that the
   * DipoleSystem can be completely restored?
   */
  bool reversible;

  /**
   * Can the emitter model guarantee that this Emission is failsafe?
   */
  bool failsafe;

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

  /**
   * Standard debug function to be called from within a debugger.
   */
  void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Emission & operator=(const Emission &);

};

/**
 * Special container class to help selecting the emission with largest scale.
 */
class EmSel {

public:

  /**
   * The only constructor requires an Energy variable where the
   * largest scale so far is modified.
   */
  EmSel(Energy & rhom): rhomax(rhom) {}

  /**
   * Conditionally assing the given Emission if it is the largest
   * scale so far.
   */
  EmSel & operator=(EmPtr e) {
    if ( !sel || ( e && e->rho > sel->rho ) ) sel = e;
    if ( sel ) rhomax = max(sel->rho, rhomax);
    return *this;
  }

  operator EmPtr () {
    return sel;
  }

private:

  /**
   * The selected Emission.
   */
  EmPtr sel;

  /**
   * The Energy variable keeping track of the highest scale so far.
   */
  Energy & rhomax;

};

}

#endif /* Ariadne5_Emission_H */
