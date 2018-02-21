// -*- C++ -*-
//
// Decayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Decayer_H
#define ThePEG_Decayer_H
// This is the declaration of the Decayer class.

#include "ThePEG/Config/ThePEG.h"
#include "Decayer.fh"
#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/MatrixElement/Amplitude.h"

namespace ThePEG {

/**
 * Decayer is an abstract base class to specify objects modelling the
 * decay of a particle.
 *
 * @see \ref DecayerInterfaces "The interfaces"
 * defined for Decayer.
 * @see ParticleData
 * @see DecayMode
 */
class Decayer: public HandlerBase {

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const = 0;

  /**
   * Return true if this Decayer need to access the full current step
   * when a particle is decayed. If true is returned the standard
   * Decay Handler will call the decay(const DecayMode&,const
   * Particle&,Step&) function rather than the decay(const
   * DecayMode&,const Particle&) function.
   */
  virtual bool needsFullStep() const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm,
			       const Particle & p) const = 0;

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
			       Step & step) const;

  /**
   * Calculate branching ratio. If this model has any oppinions on the
   * decay rate for a given decay mode \a dm, for a given particle
   * type \a pd, this method may be overriden to return this
   * oppinion. \a oldbrat is the branching ratio specified in the
   * DecayMode or by asking the WidthGenerator of \a pd.
   */
  virtual double brat(const DecayMode & dm, const ParticleData & pd,
		      double oldbrat) const;

  /**
   * Calculate branching ratio. If this model has any oppinions on the
   * decay rate for a given decay mode \a dm, for a given particle \a
   * p, this method may be overriden to return this oppinion. \a
   * oldbrat is the branching ratio specified in the DecayMode or by
   * asking the WidthGenerator of \a p.
   */
  virtual double brat(const DecayMode & dm, const Particle & p,
		      double oldbrat) const;

  /**
   * Produce the children. Can be used by sub-class decay() functions
   * to produce instances of the children. This default implementation
   * just calls the produceProducts() of the specified decay products
   * in DecayMode object, \a dm.
   */
  virtual ParticleVector getChildren(const DecayMode & dm,
				     const Particle & parent) const;

  /**
   * Boost the decay products. Can be used by sub-classes to perform
   * the final boost back from the parents cms. This default version
   * does just that.
   */
  virtual void finalBoost(const Particle & parent,
			  const ParticleVector & children) const;

  /**
   * Set the scales. Can be used by sub classes to set the production
   * scale of the children. This default version sets the scale to the
   * parents mass.
   */
  virtual void setScales(const Particle & parent,
			 const ParticleVector & children) const;
  //@}

  /**
   * Return an amplitude associated with this decay matrix
   * element. May return null.
   */
  Ptr<Amplitude>::pointer amplitude() const { return theAmplitude; }

  /**
   * Static function to administer the decay of a \a particle. The
   * children are properly added to the particle and to the given \a
   * step. Maximum \a maxtry attempts to decay the particle is
   * performed.
   * @return the produced children.
   */
  static ParticleVector
  DecayParticle(tPPtr parent, Step & step, long maxtry = 1000);

  /**
   * Exception class used if something goes wrong in DecayParticle().
   */
  struct DecayFailure: public Exception {};

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<Decayer> initDecayer;

  /**
   *  Private and non-existent assignment operator.
   */
  Decayer & operator=(const Decayer &);

  /**
   * A possible null pointer to an amplitude associated with this
   * matrix element.
   */
  Ptr<Amplitude>::pointer theAmplitude;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of Decayer. */
template <>
struct BaseClassTrait<Decayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of Decayer. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  Decayer class. */
template <>
struct ClassTraits<Decayer>: public ClassTraitsBase<Decayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Decayer"; }
};

/** @endcond */

}

#endif /* ThePEG_Decayer_H */
