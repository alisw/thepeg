// -*- C++ -*-
//
// CompSelector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_CompSelector_H
#define THEPEG_CompSelector_H
//
// This is the declaration of the CompSelector class.
//

#include "ThePEG/Utilities/Selector.h"

namespace ThePEG {

/**
 * The CompSelector class works like the Selector class in that it can
 * be used to randomly select objects according to associated
 * probabilities. In addition, the CompSelector class is able to
 * handle the case where the associated probabilities are
 * overestimates and the selected object will be discarded according
 * to some weight. If then a weight above one is encountered, this
 * means that the overestimated probability for the selected object
 * was wrong and it should in fact have been higher. If this happens,
 * the CompSelecteor will go into compensation mode, which means that
 * the selected object will be oversampled a period after the
 * violation to compensate for having been undersampled before. Also
 * the associated probability is adjusted to reflect the new
 * overestimate.
 *
 * The available functions are not as many as in Selector, and some of
 * the works somewhat differently. Before starting sampling the
 * objects should be added to a CompSelector object with the insert()
 * function. To selct an object the select() function should be
 * used. After that the weight with which the object should be
 * accepted should be presented with the reweight() function which
 * normally returns zero. If, however, the weight is larger than unity
 * the new overestimated probability is returned and the CompSelector
 * enters the compensating mode. Note that the weight is passed as a
 * reference and may be changed in by the reweight function if in the
 * compensating mode.
 */
template <typename T, typename WeightType = double>
class CompSelector {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. The optional argument gives the margin
   * used to get a new overestimated probability for an object when
   * entering compensation mode.
   */
  CompSelector(double newMargin = 1.1, double newTolerance = 1.0e-6)
    : N(0), last(), theMargin(newMargin), theTolerance(newTolerance) {}
  //@}

public:

  /** @name The main function controlling the selection and compensation. */
  //@{
  /**
   * Insert an object given a probability for this object. If the
   * probability is zero or negative, the object will not be inserted
   * and the probability itself is returned. Otherwise the sum of
   * probabilities so far is returned. Note that if selection has
   * already started and this CompSelector is in compensating mode, it
   * will immediately leave this mode and the selection procedure will
   * start from scratch.
   */
  WeightType insert(WeightType d, const T & t) {
    reset();
    return selector.insert(d, t);
  }

  /**
   * Selct an object randomly. Given a random number generator which
   * generates flat random numbers in the interval ]0,1[ with the
   * <code>operator()()</code> function, select an object according to
   * the individual probabilities specified when they were
   * inserted. If the generated number is outside the allowed range or
   * the Selector is empty, a range_error will be thrown. The
   * generator should have a push_back function which will be used
   * push back a uniform random number in the interval ]0,1[
   * calculated from the fraction of rnd which was in the range of the
   * selected object.
   */
  template <typename RNDGEN>
  T & select(RNDGEN & rnd) {
    ++N;
    if ( !compensating() ) last = selector.select(rnd);
    return last;
  }

  /**
   * Report the weight associated with the last selected
   * object. Returns the zero if weight was below unity, otherwise the
   * compensation mode will be entered and the new overestimated
   * probabilty for the last selected object will be returned.
   */
  WeightType reweight(double & weight) {
    if ( abs(weight) > 1.0 + tolerance() ) {
      // Retrieve the old overestimate of the object by seing how much
      // the summed weights are decreased when removing the object.
      WeightType oldtot = selector.sum();
      WeightType oldmax = oldtot - selector.erase(last);
      WeightType newmax = oldmax*abs(weight)*margin();
      WeightType newtot = selector.insert(newmax, last);
      double rat = newmax/oldmax;
      
      // Setup the new compensation level.
      Level level;
      level.weight = 1.0/rat;
      level.lastN = long(N*newtot/oldtot);
      
      // If we are already compensating, reweight the previous
      // compensation levels.
      for ( int i = 0, M = levels.size(); i < M; ++i ) {
	levels[i].lastN = long(levels[i].lastN*newtot/oldtot);
	levels[i].weight /= rat;
      }
      levels.push_back(level);
      weight /= rat;
      return newmax;
    }
    
    // If we are compensating we should only accept the selection if the
    // weight is above the previous overestimate.
    if ( compensating() ) if ( abs(weight) < levels.back().weight ) weight = 0.0;
    
    return WeightType();
  }

  /**
   * Exit compensation mode and start selection procedure from
   * scratch.
   */
  void reset() {
    N = 0;
    levels.clear();
    last = T();
  }

  /**
   * Erases all objects.
   */
  void clear() {
    selector.clear();
    reset();
  }

  /**
   * Set the margin used to get a new overestimated probability for an
   * object when entering compensation mode.
   */
  void margin(double m) { theMargin = m; }

  /**
   * Set the tolerance for how much a weight is allowed to be
   * larger than unity before starting the compensation.
   */
  void tolerance(double t) { theTolerance = t; }
  //@}


  /** @name Simple access functions. */
  //@{
  /**
   * Return true if this CompSelector is in a compensating state.
   */
  bool compensating() {
    // Leave all levels which has reached there 'expiry date'.
    while ( levels.size() && levels.back().lastN < N ) levels.pop_back();
    return !levels.empty();
  }

  /**
   * If in a compensating mode, return the number of selection needed
   * before exiting this mode.
   */
  long compleft() const { return levels.empty()? 0: levels.back().lastN - N; }

  /**
   * Return the sum of probabilities of the objects inserted. Note
   * that probabilities specified when objects are inserted are
   * rescaled with this number to give unit probability for
   * 'select()'.
   */
  WeightType sum() const { return selector.sum(); }

  /**
   * Return the margin used to get a new overestimated probability for an
   * object when entering compensation mode.
   */
  double margin() const { return theMargin; }

  /**
   * Return the tolerance for how much a weight is allowed to be
   * larger than unity before starting the compensation.
   */
  double tolerance() const { return theTolerance; }
  //@}

  /** @name I/O functions. */
  //@{
  /**
   * Output to a stream.
   */
  template <typename OStream>
  void output(OStream & os) const {
    os << selector << N << last << theMargin << theTolerance << levels.size();
    for ( int i = 0, M = levels.size(); i < M; ++i )
      os << levels[i].lastN << levels[i].weight;
  }

  /**
   * Input from a stream.
   */
  template <typename IStream>
  void input(IStream & is) {
    long M;
    is >> selector >> N >> last >> theMargin >> theTolerance >> M;
    levels.resize(M);
    for ( int i = 0; i < M; ++i ) is >> levels[i].lastN >> levels[i].weight;
  }
  //@}

private:

  /**
   * Internal struct used for bookkeeping when compensating.
   */
  struct Level {

    /**
     * The selection number at which point this level of compensation
     * is ended.
     */
    long lastN;

    /**
     * The minimum weight allowed when compensating on this level.
     */
    double weight;

  };

private:

  /**
   * The underlying selector
   */
  Selector<T,WeightType> selector;

  /**
   * The number of selections so far.
   */
  long N;

  /**
   * The last selected object.
   */
  T last;

  /**
   * The margin used to get a new overestimated probability for an
   * object when entering compensation mode.
   */
  double theMargin;

  /**
   * Set the tolerance for how much a weight is allowed to be
   * larger than unity before starting the compensation.
   */
  double theTolerance;

  /**
   * The currently active compensation levels.
   */
  vector<Level> levels;

};

/**
 * Output a Selector to a stream.
 */
template <typename OStream, typename T, typename WeightType>
inline OStream & operator<<(OStream & os,
			    const CompSelector<T,WeightType> & s) {
  s.output(os);
  return os;
}

/**
 * Input a Selector from a stream.
 */
template <typename IStream, typename T, typename WeightType>
inline IStream & operator>>(IStream & is,
			    CompSelector<T,WeightType> & s) {
  s.input(is);
  return is;
}

}

#endif /* THEPEG_CompSelector_H */
