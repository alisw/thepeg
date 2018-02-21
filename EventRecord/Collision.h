// -*- C++ -*-
//
// Collision.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Collision_H
#define ThePEG_Collision_H
// This is the decalaration of the Collision class. It

#include "EventConfig.h"
#include "Particle.h"
#include "StandardSelectors.h"
#include "ThePEG/Vectors/LorentzVector.h"
#include "ThePEG/Vectors/LorentzRotation.h"

namespace ThePEG {

/**
 * This is the decalaration of the Collision class. It contains all
 * <code>Particle</code>s produced in the generation of a collision
 * between two particles in an Event. The particles are divided into
 * <code>Step</code>s corresponding to the particles present after a
 * given step in the event generation. The collision also carries
 * information about the <code>SubProcesses</code> in the collision.
 *
 * @see Event
 * @see Step
 * @see SubProcess
 * @see Particle
 */
class Collision: public EventRecordBase {

public:

  /**
   * EventHandler is a friend of most Event classes.
   */
  friend class EventHandler;
  /** Most of the Event classes are friends with each other. */
  friend class Event;
  /** Most of the Event classes are friends with each other. */
  friend class Step;

public:

  /**
   * The standard constructor takes a pair of incoming particles as
   * argument. Optionally can be given a pointer to the Event which
   * this Collision belongs, and a pointer to the EventHandler
   * which produced this collision.
   * @param newIncoming a pair of incoming particles.
   * @param newEvent the Event to which this Collision belongs.
   * @param newHandler the handler object in charge of the generation
   * of this Collision.
   */
  Collision(const PPair & newIncoming, tEventPtr newEvent = tEventPtr(),
	    tcEventBasePtr newHandler = tcEventBasePtr()) 
    : theIncoming(newIncoming), theEvent(newEvent), theHandler(newHandler) {
    addParticle(incoming().first);
    addParticle(incoming().second);
  }

  /**
   * The destructor
   */
  ~Collision();

public:

  /**
   * Create a new step in this collision, which is a copy of
   * the last step (if any) and return a pointer to it.
   * @param newHandler the handler object in charge of generating the
   * new step.
   */
  tStepPtr newStep(tcEventBasePtr newHandler = tcEventBasePtr());

  /**
   * Add a new Step to this Collision.
   */
  void addStep(tStepPtr s);

  /**
   * Return a pointer to the EventHandler which produced this
   * Collision. May be the null pointer.
   */
  tcEventBasePtr handler() const { return theHandler; }

  /**
   * Return a pointer to the Event to which this Collision
   * belongs. May be the null pointer.
   */
  tEventPtr event() const { return theEvent; }

  /** @name Functions for accessing particles etc. */
  //@{
  /**
   * Extract particles from this Collision which satisfies the
   * requirements given by an object of the SelectorBase class.
   * @param r an output iterator specifying where the extracted
   * (pointers to) particles will be appended.
   * @param s SelectorBase object defining which particles should be
   * extracted.
   */
  template <class OutputIterator>
  void select(OutputIterator r, const SelectorBase & s) const;

  /**
   * Extract all final state particles in this Collision.
   * @param r an output iterator specifying where the extracted
   * (pointers to) particles will be appended.
   */
  template <class OutputIterator>
  void selectFinalState(OutputIterator r) const {
    select(r, SelectFinalState());
  }

  /**
   * Extract all final state particles in this Collision.
   * @return a vector of pointers to the extracted particles.
   */
  tPVector getFinalState() const {
    tPVector ret;
    selectFinalState(back_inserter(ret));
    return ret;
  }

  /**
   * Return a pointer to the primary SubProcess in this Collision. May
   * be the null pointer.
   */
  tSubProPtr primarySubProcess() const {
    return subProcesses().empty()? SubProPtr(): subProcesses().front();
  }

  /**
   * Return the possibly empty list of sub processes in this Collision.
   */
  const SubProcessVector & subProcesses() const {
    return theSubProcesses;
  }

  /**
   * Return a const pointer to the last step in this Collission.
   */
  tcStepPtr finalStep() const {
    return steps().empty()? tcStepPtr(): tcStepPtr(steps().back());
  }

  /**
   * Return a pointer to the last step in this Collission.
   */
  tStepPtr finalStep() {
    return steps().empty()? StepPtr(): steps().back();
  }

  /**
   * Return the vector of steps in this Collision. 
   */
  const StepVector & steps() const { return theSteps; }

  /**
   * Return a pointer to a given Step in this Collision.
   */
  tcStepPtr step(unsigned int i) const {
    return i < steps().size()? tcStepPtr(theSteps[i]): tcStepPtr();
  }

  /**
   * Return a reference to the pair of colliding particles in this
   * Collision.
   */

  const PPair & incoming() const { return theIncoming; }

  /**
   * Return the set of remnants in this collision. Remnants are
   * defined as the daughters of the incoming particles which are not
   * incoming particles to any SubProcess or children thereof which
   * are present in the final state.
   */
  tParticleSet getRemnants() const;

  /**
   * Return true if the given particle is a remnant of the colliding
   * particles. Calls the getRemnants method, so to check several
   * particles it is better to call getRemnants directly and check if
   * the particles are members of the resulting set by hand.
   */
  bool isRemnant(tPPtr p) const { return member(getRemnants(), p); }

  //@}

  /**
   * Return the vertex position of this Collision.
   */
  const LorentzPoint & vertex() const { return theVertex; }

  /**
   * Set the vertex position of this Collision.
   */
  void vertex(const LorentzPoint & p) { theVertex = p; }

  /**
   * Transform all particles in this Collision.
   */
  void transform(const LorentzRotation &);

  /**
   * Return the total invariant mass squared of the final-state
   * particles in this Collision.
   */
  Energy2 m2() const {
    return ( incoming().first->momentum() + incoming().second->momentum() ).m2();
  }

  /** @name Functions for removing entires from a Collision. */
  //@{
  /**
   * Remove (recursively) the decay products from a given Particle and
   * add the particle to the list of final state particles.
   */
  void removeDecay(tPPtr);

  /**
   * Remove (recursively) the given Particle from the Collision. If
   * this was the last daughter of the mother Particle, the latter is
   * added to the list of final state particles.
   */
  void removeParticle(tPPtr);

  /**
   * Remove all steps which have no new particles introduced in them.
   */
  void cleanSteps();

  /**
   * Remove the last Step in this Collision.
   */
  void popStep();

  //@}

public:

  /**
   * Standard function for writing to a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard functions for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function. @see Base::Init().
   */
  static void Init();

protected:

  /** @name Internal functions for adding and removing entires. */
  //@{
  /**
   * Add a new SubProcess to this Collision.
   */
  void addSubProcess(tSubProPtr p);

  /**
   * Remove a SubProcess from this Collision.
   */
  void removeSubProcess(tSubProPtr p);

  /**
   * Add a range of particles to this Collision.
   */
  template <class Iterator>
  void addParticles(Iterator first, Iterator last);

  /**
   * Add a particle to this Collision.
   */
  void addParticle(tPPtr p);

  /**
   * Remove a given Particle entry.
   */
  void removeEntry(tPPtr p);

  //@}

  /**
   * Return a reference to the list of all particles in this Collision.
   */
  const ParticleSet & all() const { return allParticles; }

  /**
   * Clone this Collision. This also makes clones of all steps, sub
   * processes and particles in this Collision.
   */
  CollPtr clone() const;

  /**
   * Rebind to cloned objects. When a Collision is cloned, a shallow
   * copy is done first, then all <code>Particle</code>s etc, are
   * cloned, and finally this method is used to see to that the
   * pointers in the cloned Collision points to the cloned
   * <code>Particle</code>s etc.
   */
  void rebind(const EventTranslationMap & trans);

private:

  /**
   * The pair of colliding particles.
   */
  PPair theIncoming;

  /**
   * A vector of all steps in this Collision.
   */
  StepVector theSteps;

  /**
   * A vector of all sub-processes in this Collision. The front
   * element points to the primary sub-process.
   */
  SubProcessVector theSubProcesses;

  /**
   * A set of all particles in this Collision.
   */
  ParticleSet allParticles;

  /**
   * A pointer to the Event to which this Collision belongs.
   */
  tEventPtr theEvent;

  /**
   * A pointer to the EventHandler which performed the generation
   * of this Collision.
   */
  tcEventBasePtr theHandler;

  /**
   * The vertex position of this Collision
   */
  LorentzPoint theVertex;

private:

  /**
   * Describe concrete class with persistent data.
   */
  static ClassDescription<Collision> initCollision;

  /**
   * Private default constructor must only be used by the
   * PersistentIStream class via the ClassTraits<Collision> class .
   */
  Collision() {}

  /**
   * The ClassTraits<Collision> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<Collision>;

  /**
   * The assignment operator is private and not implemented.
   */
  Collision & operator=(const Collision &);

  /** Output to a standard ostream. */
  friend ostream & operator<<(ostream & os, const Collision & c);

};

/** Output a Collision to a standard ostream. */
ostream & operator<<(ostream &, const Collision &);


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base class of Collision. */
template <>
struct BaseClassTrait<Collision,1>: public ClassTraitsType {
  /** Typedef of the first base class of Collision. */
  typedef EventRecordBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Collision class and how to create it. */
template <>
struct ClassTraits<Collision>: public ClassTraitsBase<Collision> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Collision"; }
  /** Create a Collision object. */
  static TPtr create() { return TPtr::Create(Collision()); }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Collision.tcc"
#endif

#endif /* ThePEG_Collision_H */
