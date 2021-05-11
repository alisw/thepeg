// -*- C++ -*-
//
// Event.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Event_H
#define ThePEG_Event_H
// This is the decalaration of the Event class.

#include "Particle.h"
#include "StandardSelectors.h"
#include "SubProcess.h"
#include "ThePEG/Utilities/Named.h"
#include "ThePEG/Utilities/AnyReference.h"

namespace ThePEG {

/**
 * The Event class contains all Particles produced in the generation
 * of an event. The particles are divided into Collisions
 * corresponding to the actiual collisions between incoming particles
 * in a bunch crossing.
 *
 * Event inherits from the Named which holds the name of an event.
 *
 * @see Collision
 * @see Step
 * @see SubProcess
 * @see Particle
 * @see SelectorBase
 * @see Named
 *
 */
class Event : public EventRecordBase, public Named {

public:

  /**
   * EventHandler is a friend of most Event classes.
   */
  friend class EventHandler;
  /** Most of the Event classes are friends with each other. */
  friend class Collision;

  /** Map colour lines to indices. */
  typedef map<tcColinePtr, int> ColourLineMap;

public:

  /**
   * The standard constructor for an Event takes as arguments a pair
   * of colliding particles (corresponding to the primary collision in
   * case of multiple collisions in an event). Optionally a pointer to
   * the EventHandler which performed the generation, an event name
   * and event number can be given.
   * @param newIncoming a pair of incoming particles to the prinary Collision.
   * @param newHandler the handler object in charge of the generation
   * of this Event.
   * @param newName the name of this event.
   * @param newNumber the number of this event.
   * @param weight the weight of this event
   */
  Event(const PPair & newIncoming, tcEventBasePtr newHandler = tcEventBasePtr(),
	string newName = "", long newNumber = -1, double weight = 1.0);

  /**
   * The copy constructor.
   */
  Event(const Event&);

  /**
   * The destructor.
   */
  ~Event();

  /**
   * Returns a full clone of this Event. All collisions,
   * <code>Particle</code>s etc. in this Event are cloned.
   */
  EventPtr clone() const;

public:

  /**
   * Return a pointer to the EventHandler which produced this
   * Event. May be the null pointer.
   */
  tcEventBasePtr handler() const { return theHandler; }

  /** @name Functions for accessing particles etc. */
  //@{
  /**
   * Extract particles from this event which satisfies the
   * requirements given by an object of the SelectorBase class.
   * @param r an output iterator specifying where the extracted
   * (pointers to) particles will be appended.
   * @param s SelectorBase object defining which particles should be
   * extracted.
    */
  template <class OutputIterator>
  void select(OutputIterator r, const SelectorBase & s) const;

  /**
   * Extract all final state particles in this Event.
   * @param r an output iterator specifying where the extracted
   * (pointers to) particles will be appended.
   */
  template <class OutputIterator>
  void selectFinalState(OutputIterator r) const {
    select(r, SelectFinalState());
  }

  /**
   * Extract all final state particles in this Event.
   * @param c a container where the extracted (pointers to) particles
   * will be appended.
   */
  template <class Container>
  void getFinalState(Container & c) const {
    selectFinalState(inserter(c));
  }

  /**
   * Extract all final state particles in this Event.
   * @return a vector of pointers to the extracted particles.
   */
  tPVector getFinalState() const {
    tPVector ret;
    selectFinalState(back_inserter(ret));
    return ret;
  }

  /**
   * Return a pointer to the primary Collision in this Event. May
   * be the null pointer.
   */
  tCollPtr primaryCollision() const {
    return collisions().empty() ? tCollPtr() : tCollPtr(collisions()[0]);
  }

  /**
   * Return a possibly empty list of collisions in this Event.
   */
  const CollisionVector & collisions() const { return theCollisions; }

  /**
   * Return a pointer to the primary SubProcess in the prinmary
   * Collision in this Event. May be the null pointer.
   */
  tSubProPtr primarySubProcess() const;

  /**
   * Return a reference to the pair of colliding particles in the
   * primary Collision of this Event.
   */
  const PPair & incoming() const { return theIncoming; }

  //@}

  /**
   * Create a new Collision in this event and return a pointer to it.
   */
  tCollPtr newCollision();

  /**
   * Create a new Step in the current Collision, which is a copy of
   * the last Step (if any) and return a pointer to it. If no
   * collision exists, one will be added.
   */
  tStepPtr newStep();

  /**
   * Transform all particles in this Event.
   */
  void transform(const LorentzRotation &);

  /**
   * Return the number assigned to this Event. The name is accessed
   * with the name() method of the Named base class.
   */
  long number() const { return theNumber; }

  /**
   * Return the index of the given colour line.
   */
  int colourLineIndex(tcColinePtr) const;

  /** @name Functions for removing entires from an Event. */
  //@{
  /**
   * Remove (recursively) the decay products from a given Particle and
   * add the particle to the list of final state particles.
   */
  void removeDecay(tPPtr);

  /**
   * Remove the given Particle from the Collision. If this was the
   * last daughter of the mother Particle, the latter is added to the
   * list of final state particles.
   */
  void removeParticle(tPPtr);

  /**
   * Remove all steps which have no new particles introduced in them.
   */
  void cleanSteps();

  //@}

  /**
   * Return the weight associated with this event.
   */
  double weight() const { return theWeight; }

  /**
   * Return an optional named weight associated to this event. Returns
   * 0, if no weight identified by this name is present.
   */
  double optionalWeight(const string& name) const;

  /**
   * Return the optional named weights associated to this event.
   */
  const map<string,double>& optionalWeights() const { return theOptionalWeights; }

  /**
   * Print this Event in Graphviz format on the standard output.
   */
  void printGraphviz() const;

  /**
   * Set the weight associated with this event.
   */
  void weight(double w) { theWeight = w; }

  /**
   * Set an optional named weight associated to this event.
   */
  void optionalWeight(const string& name, double value);

  /**
   * Access the optional named weights associated to this event.
   */
  map<string,double>& optionalWeights() { return theOptionalWeights; }

  /**
   * Set event info.
   */
  void setInfo(tcEventBasePtr newHandler, string newName,
	       long newNumber, double weight);

  /**
   * Add a collision to this Event.
   */
  void addCollision(tCollPtr c);

  /**
   * Set the primary collision in this Event.
   */
  void primaryCollision(tCollPtr c);

public:

  /**
   * Check for meta information
   */
  bool hasMeta(const string& id) const {
    return theMeta.find(id) != theMeta.end();
  }

  /**
   * Set meta information.
   */
  template<class T>
  void meta(const string& id, T& ref) {
    theMeta[id] = AnyReference(ref);
  }

  /**
   * Erase meta information.
   */
  void eraseMeta(const string& id) {
    theMeta.erase(id);
  }

  /**
   * Retrieve meta information.
   */
  template<class T>
  T& meta(const string& id) const {
    return theMeta.find(id)->second.cast<T>();
  }

protected:

  /**
   * Add a range of particles to this Collision.
   */
  template <class Iterator>
  void addParticles(Iterator first, Iterator last) {
    while ( first != last ) addParticle(*first++);
  }

  /**
   * Add a particle to this Collision.
   */
  void addParticle(tPPtr p);

  /**
   * Add a new SubProcess to this Event. For book keeping purposes
   * only. The sub-processes are accessed from the different
   * Collisions in this Event.
   */
  void addSubProcess(tSubProPtr p) {
    if ( p ) allSubProcesses.insert(p);
  }

  /**
   * Remove a SubProcess from this Event.
   */
  void removeSubProcess(tSubProPtr p) { allSubProcesses.erase(p); }

  /**
   * Add a new Step to this Collision. For book keeping purposes
   * only. The steps are accessed from the different Collisions in
   * this Event.
   */
  void addStep(tStepPtr s) { 
    if ( s ) allSteps.insert(s); 
  }

  /**
   * Remove a given Particle entry.
   */
  void removeEntry(tPPtr p);

  /**
   * Rebind to cloned objects. When an Event is cloned, a shallow
   * copy is done first, then all <code>Particle</code>s etc, are
   * cloned, and finally this method is used to see to that the
   * pointers in the cloned Event points to the cloned
   * <code>Particle</code>s etc.
   */
  void rebind(const EventTranslationMap & trans);

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

private:

  /**
   * The pair of colliding particles.
   */
  PPair theIncoming;

  /**
   * A vector of collisions in this Event.
   */
  CollisionVector theCollisions;

  /**
   * A set of all particles in this Event.
   */
  StepSet allSteps;

  /**
   * A set of all sub-processes in this Event.
   */
  SubProcessSet allSubProcesses;

  /**
   * A set of all particles in this Event.
   */
  ParticleSet allParticles;

  /**
   * A pointer to the EventHandler which performed the generation
   * of this Event.
   */
  tcEventBasePtr theHandler;

  /**
   * Map of all registered colour lines to their index numbers.
   */
  mutable ColourLineMap theColourLines;

  /**
   * The number assigned to this Event.
   */
  long theNumber;

  /**
   * The weight associated with this event.
   */
  double theWeight;

  /**
   * Optional named weights
   */
  map<string,double> theOptionalWeights;

  /**
   * Counter to keep track of particle numbering.
   */
  long theParticleNumber;

  /**
   * The meta information
   */
  map<string,AnyReference> theMeta;

public:

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger via the debug() function.
   */
  virtual void debugme() const;

private:

  /**
   * Describe concrete class with persistent data.
   */
  static ClassDescription<Event> initEvent;

  /**
   * Private default constructor must only be used by the
   * PersistentIStream class via the ClassTraits<Event> class .
   */
  Event() : theNumber(-1), theWeight(1.0), theParticleNumber(0) {}

  /**
   * The ClassTraits<Event> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<Event>;

  /**
   * The assignment operator is private and not implemented.
   */
  Event & operator=(const Event&) = delete;

};

/** Output a Event to a standard ostream. */
ostream & operator<<(ostream &, const Event &);

/** Print event tree in Graphviz format, ready for plotting. */
void printGraphviz(ostream &, tcEventPtr);

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base class of Event. */
template <>
struct BaseClassTrait<Event,1>: public ClassTraitsType {
  /** Typedef of the first base class of Collision. */
  typedef EventRecordBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Event class and how to create it. */
template <>
struct ClassTraits<Event>: public ClassTraitsBase<Event> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Event"; }
  /** Create a Event object. */
  static TPtr create() { return TPtr::Create(Event()); }
};

/** @endcond */

}

#include "Collision.h"

inline ThePEG::tSubProPtr ThePEG::Event::primarySubProcess() const {
  return collisions().empty() ? ThePEG::tSubProPtr() :
    ThePEG::tSubProPtr(primaryCollision()->primarySubProcess());
}

namespace ThePEG {
  template <class OutputIterator>
  void Event::select(OutputIterator r, const SelectorBase & s) const {
    if ( s.allCollisions() ) {
      for ( CollisionVector::const_iterator it = theCollisions.begin();
	    it != theCollisions.end(); ++it ) (**it).select(r, s);
    } else {
      primaryCollision()->select(r, s);
    }
  }
}

#endif /* ThePEG_Event_H */
