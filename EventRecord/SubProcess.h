// -*- C++ -*-
//
// SubProcess.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SubProcess_H
#define ThePEG_SubProcess_H
// This is the declaration of the SubProcess class.


#include <vector>
#include "ThePEG/EventRecord/Particle.h"

namespace ThePEG {

  class SubProcessGroup;

/**
 * A SubProcess object represents a hard \f$2\rightarrow n\f$
 * sub-process in a collision. It carries information about the
 * incoming and outgoing particles, as well as possible intermediate
 * ones. It also has a pointer to the MEBase object which generated
 * the sub-process.
 *
 * @see Event
 * @see Particle
 * @see SubProcessGroup
 */
class SubProcess: public EventRecordBase {

public:

  /** Most of the Event classes are friends with each other. */
  friend class Step;
  /** Most of the Event classes are friends with each other. */
  friend class Collision;
  /** Most of the Event classes are friends with each other. */
  friend class SubProcessGroup;

public:

  /**
   * Standard constructor.
   * @param newIncoming the two incoming partons.
   * @param newCollision the Collision to which this SubProcess belongs.
   * @param newHandler the MEBase object which generated this SubProcess.
   */
  SubProcess(const PPair & newIncoming,
	     tCollPtr newCollision = tCollPtr(),
	     tcEventBasePtr newHandler = tcEventBasePtr(),
	     tSubProPtr newHead = tSubProPtr(),
	     double newGroupWeight = 1.0);

  /**
   * Destructor.
   */
  virtual ~SubProcess();

  /**
   * A pointer to the MEBase object which generated this SubProcess.
   */
  tcEventBasePtr handler() const { return theHandler; }

  /**
   * A pointer to the collision to which this sub-process belongs.
   */
  tCollPtr collision() const { return theCollision; }

  /**
   * The pair of incoming partons.
   */
  const PPair & incoming() const { return theIncoming; }

  /**
   * A reference to the vector of intermediate partons.
   */
  const ParticleVector & intermediates() const { return theIntermediates; }

  /**
   * A reference to the vector of outgoing particles.
   */
  const ParticleVector & outgoing() const { return theOutgoing; }

  /**
   * Set the vector of outgoing particles.
   */
  template <class InputIterator>
  void setOutgoing(InputIterator, InputIterator);

  /**
   * Add a particle to the list of outgoing ones. If \a fixrelations
   * is true the mother daughter pointers will be set to/from the
   * incoming partons.
   */
  void addOutgoing(tPPtr p, bool fixrelations = true);

  /**
   * Change the incoming parton
   */
  void changeIncoming(tPPtr pnew, tPPtr pold); 

  /**
   * Set the vector of intermediate particles.
   */
  template <class InputIterator>
  void setIntermediates(InputIterator, InputIterator);

  /**
   * Add a particle to the list of intermediate ones. If \a fixrelations
   * is true the mother daughter pointers will be set to/from the
   * incoming partons.
   */
  void addIntermediate(tPPtr p, bool fixrelations = true);

  /**
   * Remove a particle entry from this sub-process.
   */
  void removeEntry(tPPtr p);

  /**
   * Return a clone of this sub process.
   */
  virtual SubProPtr clone() const;

  /**
   * True if a perturbative cascade has been applied to this sub
   * process.
   */
  bool decayed() const { return isDecayed; }

  /**
   * Set to true if a perturbative cascade has been applied to this
   * sub process.
   */
  void decayed(bool x) { isDecayed = x; }

  /**
   * Return the head SubProcess, if this SubProcess
   * object belongs to a SubProcessGroup. Return NULL
   * if head of a SubProcessGroup or not member of
   * a SubProcessGroup at all.
   */
  tSubProPtr head() const { return theHead; }

  /**
   * Set the head SubProcess
   */
  void head(tSubProPtr newHead) { theHead = newHead; }

  /**
   * If this SubProcess belongs to a SubProcessGroup,
   * return its relative weight w.r.t. the head's
   * weight.
   */
  double groupWeight() const { return theGroupWeight; }

  /**
   * If this SubProcess belongs to a SubProcessGroup,
   * set its relative weight w.r.t. the head's
   * weight.
   */
  void groupWeight(double w) { theGroupWeight = w; }

protected:

  /**
   * Rebind to cloned objects. When a SubProcess is cloned, a shallow
   * copy is done first, then all <code>Particle</code>s etc, are
   * cloned, and finally this method is used to see to that the
   * pointers in the cloned SubProcess points to the cloned
   * <code>Particle</code>s etc.
   */
  virtual void rebind(const EventTranslationMap & trans);


public:

  /**
   * Perform a LorentzTransformation of all particles in the sub
   * process.
   */
  virtual void transform(const LorentzRotation &);

  /**
   * Return the value of the Mandelstam variable \f$\hat{s}\f$ in this
   * SubProcess. It is calculated using the incoming particles.
   */
  Energy2 shat() const {
    return (incoming().first->momentum() + incoming().second->momentum()).m2();
  }

  /**
   * Return the value of the Mandelstam variable \f$\hat{t}\f$ in this
   * SubProcess. It is calculated using the first incoming and first outgoing
   * particle.
   */
  Energy2 that() const {
    return (incoming().first->momentum() - outgoing()[0]->momentum()).m2();
  }

  /**
   * Return the value of the Mandelstam variable \f$\hat{u}\f$ in this
   * SubProcess. It is calculated using the first incoming and last outgoing
   * particle.
   */
  Energy2 uhat() const {
    return (incoming().second->momentum() - outgoing()[0]->momentum()).m2();
  }

public:

  /**
   * Standard function for writing to a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function. @see Base::Init().
   */
  static void Init();

private:

  /**
   * A pointer to the MEBase object which generated this sub-process.
   */
  tcEventBasePtr theHandler;

  /**
   * A pointer to the collision to which this sub-process belongs.
   */
  tCollPtr theCollision;

  /**
   * The pair of incoming particles.
   */
  PPair theIncoming;

  /**
   * The vector of intermediate particles,
   */
  ParticleVector theIntermediates;

  /**
   * The vector of outgoing particles.
   */
  ParticleVector theOutgoing;

  /**
   * True if a perturbative cascade has been applied to this sub process.
   */
  bool isDecayed;

  /**
   * The head SubProcess, if this SubProcess
   * object belongs to a SubProcessGroup. NULL
   * if head of a SubProcessGroup or not member of
   * a SubProcessGroup at all.
   */
  tSubProPtr theHead;

  /**
   * If this SubProcess belongs to a SubProcessGroup,
   * this gives its relative weight w.r.t. the head's
   * weight.
   */
  double theGroupWeight;

public:

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger via the debug() function.
   */
  virtual void debugme() const;

  /**
   * Put to ostream
   */
  virtual void printMe(ostream&) const;

private:

  /**
   * Default constructor
   */
  SubProcess() : isDecayed(false), theGroupWeight(1.0) {}

  /**
   * Describe concrete class with persistent data.
   */
  static ClassDescription<SubProcess> initSubProcess;

  /**
   * The ClassTraits<SubProcess> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<SubProcess>;

  /**
   * Assignment is forbidden.
   */
  SubProcess & operator=(const SubProcess &) = delete;

};

/** Output a SubProcess to an ostream. */
ostream & operator<<(ostream &, const SubProcess &);


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base class of Collision. */
template <>
struct BaseClassTrait<SubProcess,1>: public ClassTraitsType {
  /** Typedef of the first base class of SubProcess. */
  typedef EventRecordBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SubProcess class and how to create it. */
template <>
struct ClassTraits<SubProcess>: public ClassTraitsBase<SubProcess> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::SubProcess"; }
  /** Create a SubProcess object. */
  static TPtr create() { return TPtr::Create(SubProcess()); }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "SubProcess.tcc"
#endif

#endif /* ThePEG_SubProcess_H */
