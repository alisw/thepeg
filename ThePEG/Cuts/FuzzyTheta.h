// -*- C++ -*-
//
// FuzzyTheta.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FuzzyTheta_H
#define ThePEG_FuzzyTheta_H
//
// This is the declaration of the FuzzyTheta class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Repository/EventGenerator.h"
#include <cassert>

namespace ThePEG {

  namespace CutTypes {

    /**
     * Identify an energy-type cut
     */
    struct Energy {};

    /**
     * Identify a momentum-type cut
     */
    struct Momentum {};

    /**
     * Identify a rapidity-type cut
     */
    struct Rapidity {};

    /**
     * Identify an azimuth-type cut
     */
    struct Azimuth {};

    /**
     * Identify an polar-angle-type cut
     */
    struct Polar {};

  }

/**
 * FuzzyTheta implements fuzzy cut prescriptions
 *
 * @see \ref FuzzyThetaInterfaces "The interfaces"
 * defined for FuzzyTheta.
 */
class FuzzyTheta: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FuzzyTheta();

  /**
   * The destructor.
   */
  virtual ~FuzzyTheta();
  //@}

public:

  /**
   * Return the (compact) support of the delta approximation
   * considered, given its center value. This default version assumes
   * a box approximation. All values are assumed to be in units of the
   * width considered.
   */
  virtual pair<double,double> support(double x) const {
    return make_pair(x-0.5,x+0.5);
  }

  /**
   * Return the overlap integral of the delta approximation with the
   * given box and center. This default version assumes
   * a box approximation. All values are assumed to be in units of the
   * width considered.
   */
  virtual double overlap(double x, const pair<double,double>& box) const {
    if ( x - 0.5 >= box.first && x + 0.5 <= box.second )
      return 1.;
    if ( x - 0.5 > box.second || x + 0.5 < box.first )
      return 0.;
    return min(box.second,x+0.5) - max(box.first,x-0.5);
  }

  /**
   * Return the overlap, optionally considering absolute lower and
   * upper bounds.
   */
  double overlap(double x,
		 pair<double,double> box,
		 const pair<double,double>& support) const {
    box.first = max(box.first,support.first);
    box.second = min(box.second,support.second);
    assert(x >= support.first && x <= support.second);
    assert(support.second - support.first >= 1.);
    if ( x - 0.5 < support.first )
      x = support.first + 0.5;
    if ( x + 0.5 > support.second )
      x = support.second - 0.5;
    return overlap(x,box);
  }

  /**
   * Return the bounds for an energy-type cut
   */
  pair<double,double> bounds(const CutTypes::Energy&) const {
    return make_pair(0.,generator()->maximumCMEnergy()/theEnergyWidth);
  }

  /**
   * Return the width for an energy-type cut
   */
  Energy width(const CutTypes::Energy&) const {
    return theEnergyWidth;
  }

  /**
   * Return the bounds for a momentum-type cut
   */
  pair<double,double> bounds(const CutTypes::Momentum&) const {
    return make_pair(0.,0.5*generator()->maximumCMEnergy()/theEnergyWidth);
  }

  /**
   * Return the width for a momentum-type cut
   */
  Energy width(const CutTypes::Momentum&) const {
    return theEnergyWidth;
  }

  /**
   * Return the bounds for a rapidity-type cut
   */
  pair<double,double> bounds(const CutTypes::Rapidity&) const {
    return make_pair(-Constants::MaxDouble,Constants::MaxDouble);
  }

  /**
   * Return the width for a rapidity-type cut
   */
  double width(const CutTypes::Rapidity&) const {
    return theRapidityWidth;
  }

  /**
   * Return the bounds for a azimuth-type cut
   */
  pair<double,double> bounds(const CutTypes::Azimuth&) const {
    return make_pair(0.0,2.*Constants::pi/theAngularWidth);
  }

  /**
   * Return the width for a azimuth-type cut
   */
  double width(const CutTypes::Azimuth&) const {
    return theAngularWidth;
  }

  /**
   * Return the bounds for a polar-angle-type cut
   */
  pair<double,double> bounds(const CutTypes::Polar&) const {
    return make_pair(0.0,Constants::pi/theAngularWidth);
  }

  /**
   * Return the width for a polar-type cut
   */
  double width(const CutTypes::Polar&) const {
    return theAngularWidth;
  }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isInside(const Value& v, const Value& lower, const Value& upper, double& weight) const {
    CutType type;
    Value w = width(type);
    weight *=
      overlap(v/w,pair<double,double>(lower/w,upper/w),bounds(type));
    if ( weight == 0.0 )
      return false;
    return true;
  }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isLessThan(const Value& v, const Value& upper, double& weight) const {
    CutType type;
    Value w = width(type);
    pair<double,double> b = bounds(type);
    weight *=
      overlap(v/w,pair<double,double>(b.first,upper/w),b);
    if ( weight == 0.0 )
      return false;
    return true;
  }

  /**
   * Check for value inside the given bounds and update the weight
   */
  template<class CutType, class Value>
  bool isLargerThan(const Value& v, const Value& lower, double& weight) const {
    CutType type;
    Value w = width(type);
    pair<double,double> b = bounds(type);
    weight *=
      overlap(v/w,pair<double,double>(lower/w,b.first),b);
    if ( weight == 0.0 )
      return false;
    return true;
  }

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The width to be considered for momenta
   */
  Energy theEnergyWidth;

  /**
   * The width to be considered for rapidity quantities
   */
  double theRapidityWidth;

  /**
   * The width to be considered for angular quantities
   */
  double theAngularWidth;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FuzzyTheta> initFuzzyTheta;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FuzzyTheta & operator=(const FuzzyTheta &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FuzzyTheta. */
template <>
struct BaseClassTrait<FuzzyTheta,1> {
  /** Typedef of the first base class of FuzzyTheta. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FuzzyTheta class and the shared object where it is defined. */
template <>
struct ClassTraits<FuzzyTheta>
  : public ClassTraitsBase<FuzzyTheta> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::FuzzyTheta"; }
};

/** @endcond */

}

#endif /* ThePEG_FuzzyTheta_H */
