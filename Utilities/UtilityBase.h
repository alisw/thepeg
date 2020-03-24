// -*- C++ -*-
//
// UtilityBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_UtilityBase_H
#define ThePEG_UtilityBase_H

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Utilities/Triplet.h"

namespace ThePEG {

/**
 * UtilityBase is a base class implementing a number of static utility
 * functions. It should be used as a base class to give acces to these
 * functions to a class. A class can safely multiply inherit from this
 * class as it only contains static functions.
 */
struct UtilityBase {

  /**
   * Sums the four-momentum of given container. The class
   * <code>Cont::value_type</code> must be of a type <code>T</code>
   * for which <code>ParticleTraits<T>::momentum(const T&)</code> is
   * implemented correctly.
   */
  template <typename Cont>
  static LorentzMomentum sumMomentum(const Cont & c) {
    return sumMomentum(c.begin(), c.end());
  }

  /**
   * Sums the four-momentum of the entries between first and last. The
   * class <code>Iterator::value_type</code> must be of a type
   * <code>T</code> for which <code>ParticleTraits<T>::momentum(const
   * T&)</code> is implemented correctly.
   */
  template <typename Iterator>
  static LorentzMomentum sumMomentum(Iterator first, Iterator last) {
    LorentzMomentum sum;
    typedef typename std::iterator_traits<Iterator>::value_type PType;
    typedef ParticleTraits<PType> Traits;
    
    while ( first != last ) sum += Traits::momentum(*first++);
    return sum;
  }

  /**
   * Transform the entries between \a first and \a last. The class
   * <code>Iterator::value_type</code> must be of a type
   * <code>T</code> for which <code>ParticleTraits<T>::momentum(const
   * T&)</code> is implemented correctly.
   */
  template <typename Iterator>
  static void transform(Iterator first, Iterator last,
			const LorentzRotation & boost) {
    typedef typename std::iterator_traits<Iterator>::value_type PType;
    typedef ParticleTraits<PType> Traits;
    
    while ( first != last ) Traits::transform(*first++, boost);
  }

  /**
   * Transform the entries in a container \a cont. The class
   * <code>Cont::value_type</code> must be of a type <code>T</code>
   * for which <code>ParticleTraits<T>::momentum(const T&)</code> is
   * implemented correctly.
   */
  template <typename Cont>
  static void transform(Cont & cont, const LorentzRotation & boost) {
    transform(cont.begin(), cont.end(), boost);
  }

  /**
   * Boost the two objects in the pair to their CM system. Also rotate
   * so that the first is along the z-axis. The class
   * <code>PType</code> must have
   * <code>ParticleTraits<PType>::momentum(const PType&)</code> and
   * <code>ParticleTraits<PType>::transform(PType&, const
   * LorentzRotation&)</code> implemented correctly.
   */
  template <typename PType>
  static LorentzRotation boostToCM(const pair<PType,PType> & pp);

  /**
   * Boost the three objects in the Triplet to their CM system. Also
   * rotate so that the first is along the z-axis and the second is in
   * the x-z plane with positive x. The class <code>PType</code> must
   * have <code>ParticleTraits<PType>::momentum(const PType&)</code>
   * and <code>ParticleTraits<PType>::transform(PType&, const
   * LorentzRotation&)</code> implemented correctly.
   */
  template <typename PType>
  static LorentzRotation boostToCM(const Triplet<PType,PType,PType> & pt);

  /**
   * Obtain the LorentzRotation needed to boost the two objects in the
   * pair to their CM system. Also rotate the LorentzRotation so that
   * the first is along the z-axis. The class <code>PType</code> must
   * have <code>ParticleTraits<PType>::momentum(const PType&)</code>
   * implemented correctly.
   */
  template <typename PType>
  static LorentzRotation getBoostToCM(const pair<PType,PType> & pp);

  /**
   * Obtain the LorentzRotation needed to boost the three objects in
   * the Triplet to their CM system. Also rotate the LorentzRotation
   * so that the first is along the z-axis and the secons i in the x-z
   * plane with positive x. The class <code>PType</code>
   * must have <code>ParticleTraits<PType>::momentum(const
   * PType&)</code> implemented correctly.
   */
  template <typename PType>
  static LorentzRotation getBoostToCM(const Triplet<PType,PType,PType> & pt);

  /**
   * Get the inverse boost as compared to getBoostToCM.
   */
  template <typename PType>
  static LorentzRotation getBoostFromCM(const pair<PType,PType> & pp);

  /**
   * Get the inverse boost as compared to getBoostToCM.
   */
  template <typename PType>
  static LorentzRotation getBoostFromCM(const Triplet<PType,PType,PType> & pt);

  /**
   * Boost the entries between fisrt and last into their CM system.
   * The class <code>Iterator::value_type</code> must be of a type
   * <code>T</code> for which <code>ParticleTraits<T>::momentum(const
   * T&)</code> and <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   */
  template <typename Iterator>
  static LorentzRotation boostToCM(Iterator first, Iterator last) {
    return boostToCM(first, last, last, last);
  }

  /**
   * Boost the entries between fisrt and last into their CM system. If
   * zAxis != last, also rotate the entries so that zAxis becomes
   * paralell to the z-axis. The class
   * <code>Iterator::value_type</code> must be of a type
   * <code>T</code> for which <code>ParticleTraits<T>::momentum(const
   * T&)</code> and <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   */
  template <typename Iterator>
  static LorentzRotation boostToCM(Iterator first, Iterator last, Iterator zAxis) {
    return boostToCM(first, last, zAxis, last);
  }

  /**
   * Boost the entries between fisrt and last into their CM system. If
   * zAxis != last, also rotate the entries so that zAxis becomes
   * paralell to the z-axis. Also, if xzPlane != last, rotate the
   * entries so that xzPlane is placed in the xz-plane.  The class
   * <code>Iterator::value_type</code> must be of a type
   * <code>T</code> for which <code>ParticleTraits<T>::momentum(const
   * T&)</code> and <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   */
  template <typename Iterator>
  static LorentzRotation boostToCM(Iterator first, Iterator last,
				   Iterator zAxis, Iterator xzPlane);

  /**
   * Rotate p to the z-axis and boost it to its CMS, then boost it
   * along the z-axis and rotate it so that it ends up with momentum
   * q. If p is massless - simply set its momentum. The class
   * <code>PType</code> must have
   * <code>ParticleTraits<PType>::momentum(const PType&)</code>
   * implemented correctly.
   */
  template <typename PType>
  static void setMomentum(PType & p, const Momentum3 & q);

  /**
   * Boost p along the z-axis and rotate it so that, if it was
   * previously at rest, it ends up with momentum q. If p is massless
   * - simply set its momentum to q. The class
   * <code>PType</code> must have
   * <code>ParticleTraits<PType>::momentum(const PType&)</code>
   * implemented correctly.
   */
  template <typename PType>
  static void setMomentumFromCMS(PType & p, const Momentum3 & q);

  /**
   * Rotate the range of particles so their sum is along z-axis and
   * boost them to their CMS, then boost them along the z-axis and
   * rotate them so that they end up with total momentum q. The class
   * <code>Iter::value_type</code> must be of a type <code>T</code>
   * for which <code>ParticleTraits<T>::momentum(const T&)</code> and
   * <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   */
  template <typename Iter>
  static void setMomentum(Iter first, Iter last, const Momentum3 & q);

  /**
   * Rotate the range of particles so their sum is along z-axis then
   * boost them along the z-axis and rotate them so that they end up
   * with total momentum q. If a single boost does not succeed to
   * obtain the required precision within eps times the total energy,
   * the boost is redone. The class <code>Iter::value_type</code> must
   * be of a type <code>T</code> for which
   * <code>ParticleTraits<T>::momentum(const T&)</code> and
   * <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   */
  template <typename Iter>
  static void setMomentum(Iter first, Iter last,
			  const Momentum3 & q, double eps);

  /**
   * Boost the range of particles along the z-axis and rotate them so
   * that, if they were previously in their rest frame, they end up
   * with total momentum q. The class <code>Iter::value_type</code> must
   * be of a type <code>T</code> for which
   * <code>ParticleTraits<T>::momentum(const T&)</code> and
   * <code>ParticleTraits<T>::transform(T&, const
   * LorentzRotation&)</code> are implemented correctly.
   * @param first iterator pointing to the first particle in the range.
   * @param last iterator indicating the end of the range.
   * @param m2 the invariant mass squared of the particles.
   * @param q final summed momentum of the particles.
   */
  template <typename Iter>
  static void setMomentumFromCMS(Iter first, Iter last,
				 Energy2 m2, const Momentum3 & q);

  /**
   * Return the transformation needed to rotate \a p to the z-axis and
   * boost it to its CMS, then boost it along the z-axis and rotate it
   * so that it ends up with momentum \a q. The class
   * <code>PType</code> must have
   * <code>ParticleTraits<PType>::momentum(const PType&)</code>
   * implemented correctly. <b>Warning</b> This function only works
   * properly if \a p has a well defined direction in both polar and
   * azimuth angles.
   * \deprecated{Use getTransformToMomentum() instead.}
   */
  template <typename PType>
  static LorentzRotation transformToMomentum(const PType & p, 
					     const Momentum3 & q) {
    typedef ParticleTraits<PType> Traits;
    LorentzMomentum q4(q, sqrt(q.mag2() + Traits::momentum(p).m2()));
    return transformToMomentum(p, q4);
  }

  /**
   * Return the transformation needed to rotate \a p to the z-axis and
   * boost it to its CMS, then boost it along the z-axis and rotate it
   * so that it ends up with momentum \a q. The class <code>PType</code>
   * must have <code>ParticleTraits<PType>::momentum(const
   * PType&)</code> implemented correctly. <b>Warning</b> This
   * function only works properly if \a p has a well defined direction
   * in both polar and azimuth angles.
   * \deprecated{Use getTransformToMomentum() instead.}
   */
  template <typename PType>
  static LorentzRotation transformToMomentum(const PType & p, 
					     const LorentzMomentum & q) {
    return transformFromCMS(q)*transformToCMS(p);
  }

  /**
   * Return a transformation appropriate for transforming \a p to have
   * the momentum \a q. The transformation is done so that the
   * auxiliary vector \a k is left unchanged.
   */
  template <typename PType>
  static LorentzRotation getTransformToMomentum(const PType & p, 
						const LorentzMomentum & q,
						const LorentzMomentum & k) {
    typedef ParticleTraits<PType> Traits;
    LorentzMomentum k0 = Traits::momentum(p) - k;
    LorentzMomentum k1 = Traits::momentum(q) - k;
    return getBoostFromCM(make_pair(k1, k))*getBoostToCM(make_pair(k0, k));
  }

  /**
   * Return a transformation appropriate for transforming \a p to have
   * the momentum \a q. The transformation is done so that the
   * auxiliary vector \a k is left unchanged.
   */
  template <typename PType>
  static LorentzRotation getTransformToMomentum(const PType & p, 
						const Momentum3 & q,
						const LorentzMomentum & k) {
    typedef ParticleTraits<PType> Traits;
    LorentzMomentum q4(q, sqrt(q.mag2() + Traits::momentum(p).m2()));
    return getTransformToMomentum(p, q4, k);
  }

  /**
   * Create a rotation corresponding to transforming p to its current
   * value from its CMS by first boosting along the z-axis and then
   * rotating. The class <code>LV</code> must have methods
   * <code>rho()</code> and <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformFromCMS(const LV & p);

  /**
   * Create a rotation corresponding to transforming sum to its
   * current value from its CMS, with zAxis along the z-axis in that
   * CMS frame. The class <code>LV</code> must have methods
   * <code>rho()</code>, <code>phi()</code> <code>theta()</code> and
   * <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformFromCMS(const LV & sum, LV zAxis);

  /**
   * Create a rotation corresponding to transforming sum to its
   * current value from its CMS, with zAxis along the z-axis and
   * xyPlane in the x-y plane in that CMS frame. The class
   * <code>LV</code> must have methods <code>rho()</code>,
   * <code>phi()</code> <code>theta()</code> and <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformFromCMS(const LV & sum,
					  const LV & zAxis, LV xyPlane);

  /**
   * Create a rotation which would transform sum to its CMS frame with
   * zAxis along the z-axis in that frame. The class <code>LV</code>
   * must have methods <code>rho()</code>, <code>phi()</code>
   * <code>theta()</code> and <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformToCMS(const LV & sum, LV zAxis);

  /**
   * Create a rotation which would transform sum to its CMS frame
   * first rotating it to the z-axis and then boost it along the
   * z-axis. The class <code>LV</code> must have methods
   * <code>rho()</code>, <code>phi()</code> <code>theta()</code> and
   * <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformToCMS(const LV & p);

  /**
   * Create a rotation which would transform sum to its CMS frame with
   * zAxis along the z-axis and xyPlane in the x-y plane in that
   * frame. The class <code>LV</code> must have methods
   * <code>rho()</code>, <code>phi()</code> <code>theta()</code> and
   * <code>e()</code>.
   */
  template <typename LV>
  static LorentzRotation transformToCMS(const LV & sum,
					const LV & zAxis, LV xyPlane);

  /**
   * Add the elements in Cont2 to Cont1, appending them to the end if
   * possible.
   */
  template <typename Cont1, typename Cont2>
  static void add(Cont1 & c1, const Cont2 & c2);

};

/** Concrete class with UtilityBase as base class. */
struct Utilities: public UtilityBase {};

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "UtilityBase.tcc"
#endif

#endif /* ThePEG_UtilityBase_H */
