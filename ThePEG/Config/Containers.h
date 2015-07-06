// -*- C++ -*-
//
// Containers.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Containers_H
#define ThePEG_Containers_H

/** \file
 * This file defines a number of containers. Some are just typedefs of
 * std containers, while others are wrappers around
 * std containers introduced in the hope of reducing the
 * amount of debugging and code duplication.
 *
 * Do not make changes in this file. If you need to modify any of the
 * standard containers used in ThePEG, edit a copy of this file and
 * include it in an alternative config file which can be included in
 * the main ThePEG.h config file
 * using the macro ThePEG_ALTERNATE_CONFIG.
 */

#include "ThePEG/Utilities/UnitIO.h"

namespace ThePEG {

/** A set of pointers to ParticleData objects. */
ThePEG_DECLARE_SET(PDPtr,ParticleDataSet);

/** A vector of pointers to ParticleData objects. */
typedef vector<PDPtr> PDVector;

/** A vector of pointers to const ParticleData objects. */
typedef vector<cPDPtr> cPDVector;

/** A vector of transient pointers to ParticleData objects. */
typedef vector<tPDPtr> tPDVector;

/** A vector of transient pointers to const ParticleData objects. */
typedef vector<tcPDPtr> tcPDVector;

/** A set of pointers to MatcherBase objects. */
ThePEG_DECLARE_SET(PMPtr,MatcherSet);

/** A set of pointers to DecayMode objects */
ThePEG_DECLARE_SET(DMPtr,DecayModeSet);

/** A set of pointers to InterfacedBase objects */
ThePEG_DECLARE_SET(IBPtr,ObjectSet);

/** A set of pointers to InterfacedBase objects */
ThePEG_DECLARE_SET(IBPtr,DependencySet);

/** A map relating integers to ParticleData objects */
ThePEG_DECLARE_MAP(long,PDPtr,ParticleMap);

/** A map relating character strings to InterfacedBase objects */
ThePEG_DECLARE_MAP(string,IBPtr,ObjectMap);

/** A map relating InterfacedBase objects to corresponding
 * DependencySet containers */
ThePEG_DECLARE_MAP(IBPtr,DependencySet,DependencyMap);

/** A vector of pointers to InterfacedBase objects. */
typedef vector<IBPtr> IVector;

/** A vector of pointers to const InterfacedBase objects. */
typedef vector<cIBPtr> CIVector;

/** A vector of pointers to Particle objects. */
typedef vector<PPtr> ParticleVector;

/** A vector of pointers to Particle objects. */
typedef vector<PPtr> PVector;

/** A vector of pointers to const Particle objects. */
typedef vector<cPPtr> cPVector;

/** A vector of transient pointers to Particle objects. */
typedef vector<tPPtr> tPVector;

/** A vector of transient pointers to const Particle objects. */
typedef vector<tcPPtr> tcPVector;

/** A list of pointers to Particle objects. */
typedef list<PPtr> ParticleList;

/** A list of pointers to Particle objects. */
typedef list<PPtr> PList;

/** A list of pointers to const Particle objects. */
typedef list<cPPtr> cPList;

/** A list of transient pointers to Particle objects. */
typedef list<tPPtr> tPList;

/** A list of transient pointers to const Particle objects. */
typedef list<tcPPtr> tcPList;

/** A map relating character strings to bare pointers to InterfaceBase objects */
ThePEG_DECLARE_MAP(string,const InterfaceBase *,InterfaceMap);

/** A rebinder for InterfacedBase objects. */
typedef Rebinder<InterfacedBase> TranslationMap;

/** A map relating character strings to EventGenerator objects */
ThePEG_DECLARE_MAP(string,EGPtr,GeneratorMap);

/** A vector of pointers to AnalysisHandler objects. */
typedef vector<AnaPtr> AnalysisVector;

/** A pair of pointers to ParticleData objects. */
typedef pair<PDPtr, PDPtr> PDPair;

/** A pair of pointers to const ParticleData objects. */
typedef pair<cPDPtr, cPDPtr> cPDPair;

/** A pair of transient pointers to ParticleData objects. */
typedef pair<tPDPtr, tPDPtr> tPDPair;

/** A pair of transient pointers to const ParticleData objects. */
typedef pair<tcPDPtr, tcPDPtr> tcPDPair;

/** A pair of pointers to Particle objects. */
typedef pair<PPtr, PPtr> PPair;

/** A pair of pointers to const Particle objects. */
typedef pair<cPPtr, cPPtr> cPPair;

/** A pair of transient pointers to const Particle objects. */
typedef pair<tPPtr, tPPtr> tPPair;

/** A pair of transient pointers to const Particle objects. */
typedef pair<tcPPtr, tcPPtr> tcPPair;

/** An Interval in scale. */
typedef Interval<Energy2> SInterval;

/** A vector of intervals of scales. */
typedef vector<SInterval> SIntervalVector;

/** A vector of pairs of transient pointers to PartonBins. */
typedef vector<tPDPair> tPartonPairVec;

/** A pair of transient pointers to ColourLine objects. */
typedef pair<tColinePtr,tColinePtr> tColinePair;

/** A set of pointers to DecayMode objects. */
ThePEG_DECLARE_SET(tDMPtr,DecaySet);

/** A set oc character strings. */
ThePEG_DECLARE_SET(string,StringSet);

/** A vector of energies. */
typedef vector<Energy> EnergyVector;

/** A vector of pointers to EventInfoBase objects. */
typedef vector<EIPtr> EIVector;

/** A vector of doubles. */
typedef vector<double> DVector;

/** A pair of doubles. */
typedef pair<double,double> DPair;

/** @name Global shift operators to simplify adding and removing
 * objects to containers. */
//@{
/**
 * Overload the left shift operator for vector to push_back objects to
 * a vector.
 * @param tv the vector being filled by push_back.
 * @param u the object being pushed back.
 * @return a referens to the vector.
 */
template <typename T, typename U>
vector<T> & operator<<(vector<T> & tv, const U & u) {
  tv.push_back(u);
  return tv;
}

/**
 * Overload the right shift operator for vector to pop objects from
 * a vector.
 * @param tv the vector being popped by pop_back.
 * @param u the object at the back of the vector before popping.
 * @return a referens to the vector.
 */
template <typename T, typename U>
vector<T> & operator>>(vector<T> & tv, U & u) {
  u = tv.back();
  tv.pop_back();
  return tv;
}

/**
 * Overload the left shift operator for stack to push objects to
 * a vector.
 * @param ts the stack being filled by push.
 * @param u the object being pushed.
 * @return a referens to the stack.
 */
template <typename T, typename U>
stack<T> & operator<<(stack<T> & ts, const U & u) {
  ts.push(u);
  return ts;
}

/**
 * Overload the right shift operator for stack to pop objects from
 * a stack.
 * @param ts the stack being popped.
 * @param u the object at the top of the stack before popping.
 * @return a referens to the stack.
 */
template <typename T, typename U>
stack<T> & operator>>(stack<T> & ts, U & u) {
  u = ts.top();
  ts.pop();
  return ts;
}

/**
 * Overload the left shift operator for deque to push_back objects to
 * a deque.
 * @param td the deque being filled by push_back.
 * @param u the object being pushed back.
 * @return a referens to the deque.
 */
template <typename T, typename U>
deque<T> & operator<<(deque<T> & td, const U & u) {
  td.push_back(u);
  return td;
}

/**
 * Overload the right shift operator for vector to pop objects from
 * a deque.
 * @param td the deque being popped by pop_front.
 * @param u the object at the front of the deque before popping.
 * @return a referens to the deque.
 */
template <typename T, typename U>
deque<T> & operator>>(deque<T> & td, U & u) {
  u = td.front();
  td.pop_front();
  return td;
}

/**
 * Overload the left shift operator for std::set to insert objects in
 * a set.
 * @param ts the set being filled by insert.
 * @param u the object being inserted.
 * @return a referens to the set.
 */
template <typename T, typename U>
set<T> & operator<<(set<T> & ts, const U & u) {
  ts.insert(u);
  return ts;
}
//@}

/** @name Functions for I/O of containers of objects with unit. */
//@{
/**
 * Ouput a vector of objects with the specified unit.
 * @param os the stream used for output.
 * @param v the vector to be output.
 * @param u the unit to be used.
 */
template <typename OStream, typename T, typename Alloc, typename UT>
void ounitstream(OStream & os, const vector<T,Alloc> & v, UT & u) {
  os << v.size();
  for ( typename vector<T,Alloc>::const_iterator i = v.begin();
	i != v.end(); ++i )
    os << ounit(*i, u);
}

/**
 * Input a vector of objects with the specified unit.
 * @param is the stream used for input.
 * @param v the vector to be input.
 * @param u the unit to be used.
 */
template <typename IStream, typename T, typename Alloc, typename UT>
void iunitstream(IStream & is, vector<T,Alloc> & v, UT & u) {
  typename vector<T,Alloc>::size_type l;
  is >> l;
  v.resize(l);
  for ( typename vector<T,Alloc>::iterator i = v.begin(); i != v.end(); ++i )
    is >> iunit(*i, u);
}

/**
 * Ouput a set of objects with the specified unit.
 * @param os the stream used for output.
 * @param s the set to be output.
 * @param u the unit to be used.
 */
template <typename OStream, typename T, typename CMP, typename A, typename UT>
void ounitstream(OStream & os, const set<T,CMP,A> & s, UT & u) {
  os << s.size();
  for ( typename set<T,CMP,A>::const_iterator i = s.begin(); i != s.end(); ++i )
    os << ounit(*i, u);
}

/**
 * Input a set of objects with the specified unit.
 * @param is the stream used for input.
 * @param s the set to be input.
 * @param u the unit to be used.
 */
template <typename IStream, typename T, typename CMP, typename A, typename UT>
void iunitstream(IStream & is, set<T,CMP,A> & s, UT & u) {
  s.clear();
  typename set<T,CMP,A>::size_type l;
  is >> l;
  T t;
  while ( l-- ) {
    is >> iunit(t, u);
    s.insert(t);
  }
}

/**
 * Ouput a map of keys and objects where the objects are output with
 * the specified unit.
 * @param os the stream used for output.
 * @param m the map to be output.
 * @param u the unit to be used for the mapped objects.
 */
template <typename OStream, typename K, typename T,
          typename CMP, typename A, typename UT>
void ounitstream(OStream & os, const map<K,T,CMP,A> & m, UT & u) {
  os << m.size();
  for ( typename map<K,T,CMP,A>::const_iterator i = m.begin();
	i != m.end(); ++i )
    os << i->first << ounit(i->second, u);
}

/**
 * Input a map of keys and objects where the objects are input with
 * the specified unit.
 * @param is the stream used for input.
 * @param m the map to be input.
 * @param u the unit to be used for the mapped objects.
 */
template <typename IStream, typename K, typename T,
          typename CMP, typename A, typename UT>
void iunitstream(IStream & is, map<K,T,CMP,A> & m, UT & u) {
  m.clear();
  typename map<K,T,CMP,A>::size_type l;
  is >> l;
  T t;
  K k;
  while ( l-- ) {
    is >> k >> iunit(t, u);
    m[k] = t;
  }
}
//@}

}

#endif /* ThePEG_Containers_H */
