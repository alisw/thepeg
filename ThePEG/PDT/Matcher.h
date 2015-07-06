// -*- C++ -*-
//
// Matcher.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Matcher_H
#define ThePEG_Matcher_H
// This is the declaration of the Matcher class.

#include "MatcherBase.h"

namespace ThePEG {

/**
 * Matcher is a templated class inheriting from MatcherBase. It is
 * used to conveniently create interfaced classes inheriting from
 * MatcherBase giving a class with a static T::Check() method as
 * template argument.
 *
 * @see MatcherBase
 * 
 */
template <class T>
class Matcher: public MatcherBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Destructor.
   */
  virtual ~Matcher();
  //@}

  /** @name Special clone and create functions used by the Repository. */
  //@{
  /**
   * Create and clone Matcher objects.
   */
  static PMPtr Create(const string & newName, string antiName);
  /**
   * Create and clone Matcher objects.
   */
  virtual PMPtr pmclone() const;
  //@}


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


  /** @name Virtual and static versions of the check functions. */
  //@{
  /**
   * Virtual function overriding the one in MatcherBase. Simply calls
   * Check().
   */
  virtual bool check(const ParticleData & pd) const { return T::Check(pd); }

  /**
   * Static check function. Return true if a given particle type, \a
   * pd, is matched by this Matcher, ie. if the T::Check() function of
   * the template argument returns true.
   */
  static bool Check(const ParticleData & pd) { return T::Check(pd); }
  //@}

protected:

  /**
   * Set antipartner.
   */
  static void setCC(tPMPtr pm, tPMPtr apm) { MatcherBase::setCC(pm,apm); }

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription< Matcher<T> > initMatcher;

};

/**
 * MatcherType is an empty non-polymorphic base class for all matcher
 * classes to be used as template argument to Matcher.
 */
struct MatcherType {};

/** @cond TRAITSPECIALIZATIONS */

/** This partial template specialization informs ThePEG about the base
 *  classes of Matcher<T>. */
template <typename T>
struct BaseClassTrait<Matcher<T>,1>: public ClassTraitsType {
  /** Typedef of the first base class of Matcher<T>. */
  typedef MatcherBase NthBase;
};

/** This partial template specialization informs ThePEG about the name
 *  of the Matcher<T> class. Note that the template argument class is
 *  assumed to have a specialization of ClassTraits of its own.*/
template <typename T>
struct ClassTraits< Matcher<T> >: public ClassTraitsBase< Matcher<T> > {
  /** Return a platform-independent class name */
  static string className() {
    return "ThePEG::Matcher<" + T::className() + ">";
  }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Matcher.tcc"
#endif

#endif /* ThePEG_Matcher_H */
