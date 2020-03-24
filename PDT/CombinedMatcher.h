// -*- C++ -*-
//
// CombinedMatcher.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_CombinedMatcher_H
#define ThePEG_CombinedMatcher_H
// This is the declaration of the AndMatcher, OrMatcher and NotMatcher.


#include "Matcher.h"

namespace ThePEG {

/**
 * The AndMatcher class represents the boolean <i>and</i> operation between
 * its two template argument classes of base type MatcherBase.
 *
 * @see MatcherBase
 */
template <class T1, class T2>
struct AndMatcher: public MatcherType {

  /**
   * Typedef for the class representing the matcher for the
   * charge-gonjugate of particles matched by this class.
   */
  typedef AndMatcher<typename T1::CC, typename T2::CC> CC;

  /**
   * Check match. Return true if the particle type \a pd is matched by
   * this class (ie. by both of the template argument classes).
   */
  static bool Check(const ParticleData & pd) {
    return T1::Check(pd) && T2::Check(pd);
  }

};

/**
 * The OrMatcher class represents the boolean <i>or</i> operation
 * between its two template argument classes of base type MatcherBase.
 *
 * @see MatcherBase
 */
template <class T1, class T2>
struct OrMatcher: public MatcherType {

  /**
   * Typedef for the class representing the matcher for the
   * charge-gonjugate of particles matched by this class.
   */
  typedef OrMatcher<typename T1::CC, typename T2::CC> CC;

  /**
   * Check match. Return true if the particle type \a pd is matched by
   * this class (ie. by either of the template argument classes).
   */
  static bool Check(const ParticleData & pd) {
    return T1::Check(pd) || T2::Check(pd);
  }

};

/**
 * The NotMatcher class represents the boolean <i>not</i> operation
 * for its template argument class of base type MatcherBase.
 *
 * @see MatcherBase
 */
template <class T>
struct NotMatcher: public MatcherType {

  /**
   * Typedef for the class representing the matcher for the
   * charge-gonjugate of particles matched by this class.
   */
  typedef NotMatcher<typename T::CC> CC;

  /**
   * Check match. Return true if the particle type \a pd is matched by
   * this class (ie. if it is not matched by the template argument
   * class).
   */
  static bool Check(const ParticleData & pd) {
    return !T::Check(pd);
  }

};

}

#endif /* ThePEG_CombinedMatcher_H */
