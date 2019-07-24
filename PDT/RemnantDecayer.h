// -*- C++ -*-
//
// RemnantDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_RemnantDecayer_H
#define THEPEG_RemnantDecayer_H
//
// This is the declaration of the RemnantDecayer class.
//

#include "ThePEG/PDT/Decayer.h"
#include "RemnantDecayer.fh"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/EventRecord/RemnantParticle.h"
#include "ThePEG/Handlers/PtGenerator.h"

namespace ThePEG {

/**
 * The RemnantDecayer class is the base class to be used for all
 * decayers capable of decaying a RemnantParticle object produced by a
 * SoftRemnantHandler object. A derived class must implement the
 * decay(const DecayMode &, const Particle &, Step &) function, while
 * the decay(const DecayMode &, const Particle &) function should
 * never be called.
 *
 * @see \ref RemnantDecayerInterfaces "The interfaces"
 * defined for RemnantDecayer.
 */
class RemnantDecayer: public Decayer {

public:

  /** A pointer to a PtGenerator object. */
  typedef Ptr<PtGenerator>::pointer PtGPtr;

public:

  /**
   * Enumerate the options for how to distribute recoils in the hard
   * subsystem when taking energy to produce remnants.
   */
  enum RecoilOption {
    boostAll,   /**< Boost all particles in the hard subsystem. */
    boostFinal, /**< Boost only final state particles in hard subsystem. */
    copyFinal   /**< Boost copies of final state particles in hard subsystem. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  RemnantDecayer() : respectDIS(2), theRecoilOption(copyFinal) {}

  /**
   * The destructor.
   */
  virtual ~RemnantDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Return true if this Decayer need to access the full current step
   * when a particle is decayed. If true is returned the standard
   * Decay Handler will call the decay(const DecayMode&,const
   * Particle&,Step&) function rather than the decay(const
   * DecayMode&,const Particle&) function.
   */
  virtual bool needsFullStep() const;

  /**
   * Perform a decay for a given DecayMode and a given Particle
   * instance. This version allows the decaying particle to borrow
   * energy/momentum from its sublings in the current step. This will
   * be called by the standard DecayHandler if the needsFullStep()
   * function returns true.
   *
   * @param dm   the DecayMode describing the decay.
   * @param p    the Particle instance to be decayed.
   * @param step the current step in which to find possible siblings to
   *             shuffle energy with.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p,
			       Step & step) const = 0;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}

  /**
   * Return true if this decayer can handle the extraction of the \a
   * extracted parton from the given \a particle.
   */
  virtual bool canHandle(tcPDPtr parent, tcPDPtr extracted) const;

  /**
   * Return true if this decayer can handle the extraction of the \a
   * extracted parton instance from the given \a particle instance. \a
   * pnew is the momentum of the resulting remnant. The default
   * version simply checks if the energy is positive.
   */
  virtual bool checkExtract(tcPPtr parent, tcPPtr extracted,
			    const LorentzMomentum & pnew) const;

  /**
   * Return true if this decayed can extract more than one parton from
   * a particle.
   */
  virtual bool multiCapable() const;

  /**
   * The option for how to distribute recoils in the hard subsystem
   * when taking energy to produce remnants.
   */
  RecoilOption recoilOption() const { return theRecoilOption; }

  /**
   * If true, do not boost a scattered lepton (and possible radiated
   * photons) in a DIS event, to ensure that \f$x\f$ and \f$Q^2\f$ is
   * unmodified.
   */
  int respectDISKinematics() const { return respectDIS; }

  /**
   * An object capable of generating an intrinsic transverse momentum
   * of the created remnants.
   */
  PtGPtr pTGenerator() const { return thePTGenerator; }

  /**
   * Static function to decay al remnants among the given \a
   * particles. The decay products are inserted in the \a step
   * provided.
   * @return a vector of the non-remnant particles together with the
   * remnant decay products.
   */
  static tPVector decayRemnants(const tPVector & particles, Step & step);

protected:

  /**
   * Access the RemnantData object of a \a remnant.
   */
  tRemPDPtr data(tcRemPPtr remnant) const { return remnant->remData; }

  /**
   * Access the parent  of a \a remnant.
   */
  tcPPtr parent(tcRemPPtr remnant) const { return remnant->parent; }

  /**
   * Access the vector of extracted particles of a \a remnant.
   */
  const PVector & extracted(tcRemPPtr remnant) const {
    return remnant->extracted();
  }

  /**
   * Recursively find all particles produced from an extracted parton.
   */
  virtual void fillSubSystem(tPPtr p, set<tPPtr> & sub) const;

  /**
   * Return the system of particles from the hard subsystem which may
   * be used to shuffle momenta to get the remnants on-shell. In this
   * version the particles are ordered in rapidity with the ones
   * closest to the remnant direction comes first. Other orderings can
   * be enforced by sub-classes.
   */
  virtual tPVector getSubSystem(tcPPtr parent, tPPtr parton) const;

  /**
   * Return a small boost along the z-axis. To cure rounding errors
   * when making large boosts it is sometimes necessary to correct the
   * plus (or minus) lightcone component with a small boost along the
   * z-axis. The resulting boost is constructed so that the momentum
   * \a p0 would be transformed to have the sam z-value as the
   * momentum \a p.
   */
  static LorentzRotation getZBoost(const LorentzMomentum & p0,
				   const LorentzMomentum & p);

public:

  /**
   * Exception used if getSubSystem fails.
   */
  struct SubSystemFail: public Exception {};

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
   * Return true if this object needs to be initialized before all
   * other objects because it needs to extract cuts from the event file.
   */
  virtual bool preInitialize() const;
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
  /**
   * If true, do not boost a scattered lepton (and possible radiated
   * photons) in a DIS event, to ensure that \f$x\f$ and \f$Q^2\f$ is
   * unmodified.
   */
  mutable int respectDIS;

private:

  /**
   * The option for how to distribute recoils in the hard subsystem
   * when taking energy to produce remnants.
   */
  RecoilOption theRecoilOption;

  /**
   * An object capable of generating an intrinsic transverse momentum
   * of the created remnants.
   */
  PtGPtr thePTGenerator;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractClassDescription<RemnantDecayer> initRemnantDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RemnantDecayer & operator=(const RemnantDecayer &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RemnantDecayer. */
template <>
struct BaseClassTrait<RemnantDecayer,1> {
  /** Typedef of the first base class of RemnantDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RemnantDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<RemnantDecayer>
  : public ClassTraitsBase<RemnantDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::RemnantDecayer"; }
};

/** @endcond */

}

#endif /* THEPEG_RemnantDecayer_H */
