// -*- C++ -*-
//
// SelectorBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SelectorBase_H
#define ThePEG_SelectorBase_H
// This is the declaration of the SelectorBase class.


#include "EventConfig.h"

namespace ThePEG {

/**
 * Classes derived from the <code>SelectorBase</code> class are used
 * to extract particles from an Event with
 * <code>Event::select()</code> method. There are five different kinds
 * of checks done by a selector object in the
 * <code>Event::select</code> method. If the
 * <code>allCollisions()</code> method returns false, only particles
 * which belongs to the primary collision in an event will be
 * considered for extraction. Furthermore if the
 * <code>allSteps()</code> method returns false, only particles
 * present in the final step of each collision will be considered. If
 * <code>finalState()</code> returns false, final state particles will
 * not be considered and if <code>intermediate()</code> returns false,
 * intermediate particles will not be considered. Finally among all
 * considered particles, only the ones for which the <code>check(const
 * Particle &)</code> returns true will be extracted.
 *
 *
 *
 *
 * @see StandardSelectors
 * @see Event
 * @see Collision
 * @see Step
 * @see Particle
 * 
 */
class SelectorBase {

public:

  /**
   * Virtual destructor.
   */
  virtual ~SelectorBase() {}

  /**
   * Static method corresponding to the virtual check() method.
   */
  static bool Check(const Particle &) { return true; }

  /**
   * Static method corresponding to the virtual intermediate() method.
   */
  static bool Intermediate() { return true; }

  /**
   * Static method corresponding to the virtual finalState() method.
   */
  static bool FinalState() { return true; }

  /**
   * Static method corresponding to the virtual allSteps() method.
   */
  static bool AllSteps() { return true; }

  /**
   * Static method corresponding to the virtual allCollisions() method.
   */
  static bool AllCollisions() { return true; }

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const { return Check(p); }

  /**
   * Return true if final state particles are to be considered.
   */
  virtual bool finalState() const { return FinalState(); }

  /**
   * Return true if intermediate particles should be considered.
   */
  virtual bool intermediate() const { return Intermediate(); }

  /**
   * Return true if all steps should be considered. Otherwise only the
   * last step in each collision is considered.
   */
  virtual bool allSteps() const { return AllSteps(); }

  /**
   * Return ture if all collisions should be considered. Otherwise
   * only the primary collision will be considered.
   */
  virtual bool allCollisions() const { return AllCollisions(); }

};

/**
 * The templated <code>ParticleSelector</code> class may be used to
 * implement derived classes from the <code>SelectorBase</code>
 * class. The requirement on the template class is that it implements
 * the static <code>AllCollisions()</code>, <code>AllSteps()</code>,
 * <code>FinalState()</code>, <code>Intermediate()</code> and
 * <code>Check(const Particle &)</code> (corresponding to the virtual
 * ones in <code>ParticleSelector</code>).
 */
template <class T>
struct ParticleSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual check() method.
   */
  static bool Check(const Particle & p) { return T::Check(p); }

  /**
   * Static method corresponding to the virtual intermediate() method.
   */
  static bool Intermediate() { return T::Intermediate(); }

  /**
   * Static method corresponding to the virtual finalState() method.
   */
  static bool FinalState() { return T::FinalState(); }

  /**
   * Static method corresponding to the virtual allSteps() method.
   */
  static bool AllSteps() { return T::AllSteps(); }

  /**
   * Static method corresponding to the virtual allCollisions() method.
   */
  static bool AllCollisions() { return T::AllCollisions(); }

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const { return Check(p); }

  /**
   * Return true if final state particles are to be considered.
   */
  virtual bool finalState() const { return FinalState(); }

  /**
   * Return true if intermediate particles should be considered.
   */
  virtual bool intermediate() const { return Intermediate(); }

  /**
   * Return true if all steps should be considered. Otherwise only the
   * last step in each collision is considered.
   */
  virtual bool allSteps() const { return AllSteps(); }

   /**
   * Return ture if all collisions should be considered. Otherwise
   * only the primary collision will be considered.
   */
  virtual bool allCollisions() const { return AllCollisions(); }

};  

/**
 * The SelectIfNot classes can be used to negate the meaning of
 * another SelectorBase object.
 */
class SelectIfNot: public SelectorBase {

public:

  /** Constructor taking the SelectorBase object to be negated. */
  explicit SelectIfNot(const SelectorBase & S) : s(S) {}

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const { return !s.check(p); }

  /**
   * Return true if final state particles are to be considered.
   */
  virtual bool finalState() const { return !s.finalState(); }

  /**
   * Return true if intermediate particles should be considered.
   */
  virtual bool intermediate() const { return !s.intermediate(); }

  /**
   * Return true if all steps should be considered. Otherwise only the
   * last step in each collision is considered.
   */
  virtual bool allSteps() const { return !s.allSteps(); }

  /**
   * Return ture if all collisions should be considered. Otherwise
   * only the primary collision will be considered.
   */
  virtual bool allCollisions() const { return !s.allCollisions(); }

private:

  /**
   * The selector to be negated.
   */
  const SelectorBase & s;
};

/**
 * The SelectIfBoth class can be used to combine other selector
 * objects. Particles which would be extracted with either selectors
 * will be extractor.
 */
class SelectIfBoth: public SelectorBase {

public:

  /**
   * Constructor taking two SelectorBase object to be combiden.
   */
  SelectIfBoth(const SelectorBase & S1, const SelectorBase & S2)
    : s1(S1), s2(S2) {}

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const {
    return s1.check(p) && s2.check(p);
  }

  /**
   * Return true if final state particles are to be considered.
   */
  virtual bool finalState() const {
    return s1.finalState() && s2.finalState();
  }

  /**
   * Return true if intermediate particles should be considered.
   */
  virtual bool intermediate() const {
    return s1.intermediate() && s2.intermediate();
  }

  /**
   * Return true if all steps should be considered. Otherwise only the
   * last step in each collision is considered.
   */
  virtual bool allSteps() const {
    return s1.allSteps() && s2.allSteps();
  }

  /**
   * Return ture if all collisions should be considered. Otherwise
   * only the primary collision will be considered.
   */
  virtual bool allCollisions() const {
    return s1.allCollisions() && s2.allCollisions();
  }

private:

  /**
   * One selector to be combined.
   */
  const SelectorBase & s1;

  /**
   * The other selector to be combined.
   */
  const SelectorBase & s2;

};

/**
 * The SelectIfEither class can be used to combine other selector
 * objects. Only particles which would be extracted with both selectors
 * will be extractor.
 */
class SelectIfEither: public SelectorBase {

public:

  /**
   * Constructor taking two SelectorBase object to be combiden.
   */
  SelectIfEither(const SelectorBase & S1, const SelectorBase & S2)
    : s1(S1), s2(S2) {}

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const {
    return s1.check(p) || s2.check(p);
  }

  /**
   * Return true if final state particles are to be considered.
   */
  virtual bool finalState() const {
    return s1.finalState() || s2.finalState();
  }

  /**
   * Return true if intermediate particles should be considered.
   */
  virtual bool intermediate() const {
    return s1.intermediate() || s2.intermediate();
  }

  /**
   * Return true if all steps should be considered. Otherwise only the
   * last step in each collision is considered.
   */
  virtual bool allSteps() const {
    return s1.allSteps() || s2.allSteps();
  }

  /**
   * Return ture if all collisions should be considered. Otherwise
   * only the primary collision will be considered.
   */
  virtual bool allCollisions() const {
    return s1.allCollisions() || s2.allCollisions();
  }

private:

  /**
   * One selector to be combined.
   */
  const SelectorBase & s1;

  /**
   * The other selector to be combined.
   */
  const SelectorBase & s2;

};

/** Helper function to be used together with SelectorBase objects. */
template <typename OutputIterator, typename Container>
inline void copyIfCheck(OutputIterator r, const Container & c,
			const SelectorBase & s) {
  for ( typename Container::const_iterator it = c.begin();
	it != c.end(); ++it )
    if ( s.check(**it) ) *r++ = *it;
}

}

#endif /* ThePEG_SelectorBase_H */
